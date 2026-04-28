import argparse
import numpy as np
import scanpy as sc
import decoupler as dc


def run_cell_type_annotation(
    adata_path,
    output_path,
    groupby="leiden",
    species_column="mouse",
    min_n=3,
    n_ctypes=10,
    annotation_key="ora_cell_type",
    top_celltypes_key="ora_top_celltypes",
):
    adata = sc.read_h5ad(adata_path)

    markers = dc.get_resource("PanglaoDB")
    markers = markers[markers[species_column] == True].copy()
    markers = markers[~markers.duplicated(["cell_type", "genesymbol"])].copy()
    markers["genesymbol"] = markers["genesymbol"].str.capitalize()

    dc.run_ora(
        mat=adata,
        net=markers,
        source="cell_type",
        target="genesymbol",
        min_n=min_n,
        verbose=False,
    )

    acts = dc.get_acts(adata, obsm_key="ora_estimate")

    acts_v = acts.X.ravel()
    max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
    acts.X[~np.isfinite(acts.X)] = max_e

    df = dc.rank_sources_groups(
        acts,
        groupby=groupby,
        reference="rest",
        method="t-test_overestim_var",
    )

    ctypes_dict = (
        df.groupby("group")
        .head(n_ctypes)
        .groupby("group")["names"]
        .apply(list)
        .to_dict()
    )

    top_annotation = df.groupby("group").head(1).set_index("group")["names"].to_dict()

    adata.obs[annotation_key] = adata.obs[groupby].astype(str).map(top_annotation)
    adata.obs[top_celltypes_key] = (
        adata.obs[groupby]
        .astype(str)
        .map(lambda x: "; ".join(ctypes_dict.get(x, [])))
    )

    adata.uns["ora_cell_type_rankings"] = df
    adata.uns["ora_top_celltypes_by_cluster"] = ctypes_dict

    adata.write_h5ad(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--groupby", default="leiden")
    parser.add_argument("--species-column", default="mouse")
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--n-celltypes", type=int, default=10)
    parser.add_argument("--annotation-key", default="ora_cell_type")
    parser.add_argument("--top-celltypes-key", default="ora_top_celltypes")

    args = parser.parse_args()

    run_cell_type_annotation(
        adata_path=args.input,
        output_path=args.output,
        groupby=args.groupby,
        species_column=args.species_column,
        min_n=args.min_n,
        n_ctypes=args.n_celltypes,
        annotation_key=args.annotation_key,
        top_celltypes_key=args.top_celltypes_key,
    )
