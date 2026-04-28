#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc


log = logging.getLogger(__name__)


def setup_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )


DEFAULT_SAMPLES = [
    "Smp1", "Smp2", "Smp3", "Smp4",
    "Smp5", "Smp6", "Smp7",
    "SmpA", "SmpB", "SmpC", "SmpD", "SmpE",
]

_DEFAULT_MOUSE_MAP: dict[str, str] = {
    "Smp1": "Mouse#1", "Smp2": "Mouse#1", "Smp3": "Mouse#1", "Smp4": "Mouse#1",
    "Smp5": "Mouse#2", "Smp6": "Mouse#2", "Smp7": "Mouse#2",
    "SmpA": "Mouse#3", "SmpB": "Mouse#3", "SmpC": "Mouse#3", "SmpD": "Mouse#3",
    "SmpE": "Mouse#4",
}

_STALE_COLS = [
    "barcode", "pGen", "MinCut", "polylox",
    "polylox_new", "Polylox_BC", "MouseBC", "cb",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assign Polylox barcodes to cells using sample-specific, mouse-consistent matching.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--adata", required=True, help="Input .h5ad file.")
    p.add_argument("--output", required=True, help="Output .h5ad file.")
    p.add_argument(
        "--base-path",
        default="/home/jxfeng/DataProcessing/PolyloxData/2025_Processed_Polylox/CSV_files_Snakemake/",
        help="Directory containing <sample>_seg_assemble.csv files.",
    )
    p.add_argument(
        "--rare-barcode-list",
        default="/home/jxfeng/DataProcessing/PolyloxData/2025_Processed_Polylox/merged_polylox_data_pgen.txt",
    )
    p.add_argument(
        "--barcode-lib",
        default="/home/jxfeng/DataProcessing/PolyloxData/polylox_barcodelib.txt",
    )
    p.add_argument(
        "--min-recomb",
        default="/home/jxfeng/DataProcessing/PolyloxData/PolyloxMatlab_2024/min_reconb_list.txt",
    )
    p.add_argument(
        "--sample-mouse-map",
        default=None,
        help="Optional TSV without header mapping sample to mouse ID. Columns: sample<TAB>mouse.",
    )
    p.add_argument("--samples", nargs="*", default=DEFAULT_SAMPLES)
    p.add_argument("--cb-length", type=int, default=16)
    p.add_argument("--missing-barcode-pgen", type=float, default=1e-8)
    p.add_argument("--default-mincut", type=int, default=2)
    p.add_argument("--no-reverse-complement", action="store_true")
    return p.parse_args()


_RC_TABLE = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(seq: str) -> str:
    if pd.isna(seq):
        return seq
    seq = str(seq).upper()
    return seq.translate(_RC_TABLE)[::-1]


def load_sample_mouse_map(path: str | None) -> dict[str, str]:
    if path is None:
        log.info("No --sample-mouse-map provided; using hardcoded default.")
        return _DEFAULT_MOUSE_MAP.copy()

    df = pd.read_csv(path, sep="\t", header=None, names=["sample", "mouse"])
    df["sample"] = df["sample"].astype(str)
    df["mouse"] = df["mouse"].astype(str)
    mapping = df.dropna().drop_duplicates(subset=["sample"]).set_index("sample")["mouse"].to_dict()
    log.info("Loaded %d sample-to-mouse mappings from %s", len(mapping), path)
    return mapping


def validate_paths(*paths: Path) -> None:
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(f"Required file not found: {p}")


def strip_sample_suffix(obs_names: pd.Index, cb_length: int) -> pd.Series:
    raw = obs_names.astype(str).str.split("-").str[0]
    return raw.str[:cb_length].str.strip().str.upper()


def read_barcode_lib(barcode_lib_path: Path, min_recomb_path: Path, default_mincut: int) -> pd.DataFrame:
    barcode_lib = pd.read_csv(barcode_lib_path, sep="\t", header=None)
    min_recomb = pd.read_csv(min_recomb_path, sep="\t", header=None)

    n = min(len(barcode_lib), len(min_recomb))
    if len(barcode_lib) != len(min_recomb):
        log.warning(
            "Barcode library (%d rows) and MinCut list (%d rows) differ in length; using first %d rows.",
            len(barcode_lib), len(min_recomb), n,
        )

    df = barcode_lib.iloc[:n, [0]].copy()
    df.columns = ["barcode"]
    df["MinCut"] = pd.to_numeric(min_recomb.iloc[:n, 0], errors="coerce").fillna(default_mincut).astype(int)
    df["barcode"] = df["barcode"].astype(str).str.strip().str.upper()
    df = df.drop_duplicates(subset=["barcode"])
    log.info("Barcode library loaded: %d entries", len(df))
    return df


_KNOWN_BC_COLS = ["PolyloxBC", "Barcode", "barcode", "polylox"]


def detect_rare_key_column(df: pd.DataFrame) -> str:
    for col in _KNOWN_BC_COLS:
        if col in df.columns and df[col].notna().any():
            return col
    raise ValueError(
        f"No usable barcode key column found in rare barcode list. "
        f"Looked for: {_KNOWN_BC_COLS}. Available columns: {list(df.columns)}"
    )


def load_sample_barcodes(
    samples: list[str],
    base_path: Path,
    rare_barcode_path: Path,
    barcode_lib_df: pd.DataFrame,
    sample_mouse_map: dict[str, str],
    missing_barcode_pgen: float,
    default_mincut: int,
) -> dict[str, pd.DataFrame]:
    rare_df = pd.read_csv(rare_barcode_path, sep="\t", header=0)
    rare_key = detect_rare_key_column(rare_df)
    rare_df = rare_df.copy()
    rare_df[rare_key] = rare_df[rare_key].astype(str).str.strip().str.upper()

    barcodes: dict[str, pd.DataFrame] = {}

    for sample in samples:
        if sample not in sample_mouse_map:
            log.warning("Sample %s is not present in the sample-to-mouse map; skipping.", sample)
            continue

        sample_file = base_path / f"{sample}_seg_assemble.csv"
        if not sample_file.exists():
            log.warning("Sample file not found, skipping: %s", sample_file)
            continue

        df = pd.read_csv(sample_file, header=0)

        missing_cols = {"cb", "polylox"} - set(df.columns)
        if missing_cols:
            log.warning("Sample %s missing required columns %s; skipping.", sample, missing_cols)
            continue

        df = df[["cb", "polylox"]].drop_duplicates().copy()
        df["sample"] = sample
        df["mouse"] = sample_mouse_map[sample]
        df["cb"] = df["cb"].astype(str).str.strip().str.upper()
        df["polylox"] = df["polylox"].astype(str).str.strip().str.upper()
        df["pGen"] = 1.0

        df = df.merge(
            barcode_lib_df[["barcode", "MinCut"]],
            left_on="polylox",
            right_on="barcode",
            how="left",
        )

        absent_mask = df["barcode"].isna()
        df.loc[absent_mask, "pGen"] = missing_barcode_pgen

        if sample in rare_df.columns:
            rare_map = (
                rare_df[[rare_key, sample]]
                .dropna()
                .drop_duplicates(subset=[rare_key])
                .set_index(rare_key)[sample]
                .to_dict()
            )
            df["pGen"] = df["polylox"].map(rare_map).fillna(df["pGen"])

        df["MinCut"] = pd.to_numeric(df["MinCut"], errors="coerce").fillna(default_mincut).astype(int)
        df["pGen"] = pd.to_numeric(df["pGen"], errors="coerce")
        df = df[["sample", "mouse", "cb", "polylox", "pGen", "MinCut"]].drop_duplicates()
        barcodes[sample] = df
        log.info("%s: %d barcode assignments loaded", sample, len(df))

    return barcodes


def build_lookup_df(barcodes: dict[str, pd.DataFrame]) -> pd.DataFrame:
    combined = pd.concat(barcodes.values(), ignore_index=True)
    combined["sample"] = combined["sample"].astype(str)
    combined["mouse"] = combined["mouse"].astype(str)
    combined["cb"] = combined["cb"].astype(str).str.strip().str.upper()
    combined["polylox"] = combined["polylox"].astype(str).str.strip().str.upper()
    combined["pGen"] = pd.to_numeric(combined["pGen"], errors="coerce")
    combined["MinCut"] = pd.to_numeric(combined["MinCut"], errors="coerce")

    combined["_order"] = range(len(combined))
    combined = combined.sort_values(
        by=["sample", "mouse", "cb", "pGen", "MinCut", "_order"],
        ascending=[True, True, True, False, True, True],
        na_position="last",
    )

    lookup = (
        combined
        .drop_duplicates(subset=["sample", "mouse", "cb"], keep="first")
        .drop(columns="_order")
        .reset_index(drop=True)
    )

    log.info("Lookup table built: %d unique sample + mouse + cell barcode assignments", len(lookup))
    return lookup


def main() -> None:
    setup_logging()
    args = parse_args()

    adata_path = Path(args.adata).resolve()
    output_path = Path(args.output).resolve()
    base_path = Path(args.base_path).resolve()
    rare_barcode_path = Path(args.rare_barcode_list).resolve()
    barcode_lib_path = Path(args.barcode_lib).resolve()
    min_recomb_path = Path(args.min_recomb).resolve()

    validate_paths(adata_path, rare_barcode_path, barcode_lib_path, min_recomb_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sample_mouse_map = load_sample_mouse_map(args.sample_mouse_map)

    log.info("Loading AnnData: %s", adata_path)
    adata = sc.read_h5ad(adata_path)
    log.info("AnnData: %d cells x %d genes", adata.n_obs, adata.n_vars)

    if "sample" not in adata.obs.columns:
        log.error("adata.obs must contain a 'sample' column for sample-specific Polylox barcode assignment.")
        sys.exit(1)

    stale = [c for c in _STALE_COLS if c in adata.obs.columns]
    if stale:
        log.info("Dropping stale obs columns from previous run: %s", stale)
        adata.obs.drop(columns=stale, inplace=True)

    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata.obs["mouse"] = adata.obs["sample"].map(sample_mouse_map)

    unassigned_mask = adata.obs["mouse"].isna()
    if unassigned_mask.any():
        log.warning(
            "%d cells have samples absent from the sample-to-mouse map and will not receive Polylox assignments. Missing samples: %s",
            int(unassigned_mask.sum()),
            sorted(adata.obs.loc[unassigned_mask, "sample"].unique().tolist()),
        )
        adata.obs.loc[unassigned_mask, "mouse"] = "NoAssign"

    adata.obs["cb"] = strip_sample_suffix(adata.obs_names, args.cb_length).values

    log.info("Loading barcode library and per-sample barcode data.")
    barcode_lib_df = read_barcode_lib(barcode_lib_path, min_recomb_path, args.default_mincut)
    barcodes = load_sample_barcodes(
        samples=args.samples,
        base_path=base_path,
        rare_barcode_path=rare_barcode_path,
        barcode_lib_df=barcode_lib_df,
        sample_mouse_map=sample_mouse_map,
        missing_barcode_pgen=args.missing_barcode_pgen,
        default_mincut=args.default_mincut,
    )

    if not barcodes:
        log.error("No sample data loaded; cannot proceed.")
        sys.exit(1)

    lookup = build_lookup_df(barcodes)
    lookup_indexed = lookup.set_index(["sample", "mouse", "cb"])

    cell_index = pd.MultiIndex.from_frame(adata.obs[["sample", "mouse", "cb"]])
    adata.obs["polylox"] = lookup_indexed["polylox"].reindex(cell_index).to_numpy()
    adata.obs["pGen"] = lookup_indexed["pGen"].reindex(cell_index).to_numpy()
    adata.obs["MinCut"] = lookup_indexed["MinCut"].reindex(cell_index).to_numpy()

    if args.no_reverse_complement:
        adata.obs["barcode"] = adata.obs["polylox"]
    else:
        adata.obs["barcode"] = adata.obs["polylox"].apply(reverse_complement)

    n_total = adata.n_obs
    n_matched = adata.obs["polylox"].notna().sum()
    n_missing = n_total - n_matched
    pct = 100 * n_matched / n_total if n_total else 0

    log.info(
        "Barcode assignment: %d / %d cells matched (%.1f%%); %d unmatched.",
        n_matched, n_total, pct, n_missing,
    )
    log.info("Matching key: sample + mouse + cb. Duplicate resolution: highest pGen, then lowest MinCut.")

    log.info("Writing output: %s", output_path)
    adata.write_h5ad(output_path)
    log.info("Done.")


if __name__ == "__main__":
    main()
