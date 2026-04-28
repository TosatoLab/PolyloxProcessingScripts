import argparse
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


def run_doublet_detection(
    input_path,
    output_path,
    figure_path,
    batch_key="sample",
    threshold_line=0.2,
):
    adata = sc.read_h5ad(input_path)

    if adata.raw is not None:
        adata_raw = sc.AnnData(
            X=adata.raw.X,
            obs=adata.obs.copy(),
            var=adata.raw.var.copy(),
        )
    else:
        adata_raw = adata.copy()

    if batch_key not in adata_raw.obs.columns:
        batch_key = None

    sc.pp.scrublet(adata_raw, batch_key=batch_key)

    adata.obs["doublet_score"] = adata_raw.obs["doublet_score"]
    adata.obs["predicted_doublet"] = adata_raw.obs["predicted_doublet"]

    plt.figure(figsize=(6, 4))
    sns.histplot(adata.obs["doublet_score"], bins=50, kde=True)
    plt.axvline(threshold_line, color="red", linestyle="--", label=f"Threshold = {threshold_line}")
    plt.legend()
    plt.title("Doublet score distribution")
    plt.xlabel("Doublet score")
    plt.ylabel("Cell count")
    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()

    adata.write_h5ad(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-f", "--figure", default="doublet_score_distribution.png")
    parser.add_argument("--batch-key", default="sample")
    parser.add_argument("--threshold-line", type=float, default=0.2)

    args = parser.parse_args()

    run_doublet_detection(
        input_path=args.input,
        output_path=args.output,
        figure_path=args.figure,
        batch_key=args.batch_key,
        threshold_line=args.threshold_line,
    )
