import argparse
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


def run_true_barcode_pgen_heatmap(
    input_path,
    figure_path,
    table_path,
    pgen_cutoff=1e-6,
):
    adata = sc.read_h5ad(input_path)

    adata.obs["pGen"] = adata.obs["pGen"].astype(float)
    adata_subset = adata[adata.obs["pGen"] < pgen_cutoff].copy()

    if "cell_type" in adata_subset.obs.columns:
        adata_subset = adata_subset[adata_subset.obs["cell_type"] != "Doublets"].copy()

    celltype_barcode = adata_subset.obs[["barcode", "cell_type"]].copy()
    barcode_celltype_table = pd.crosstab(
        index=celltype_barcode["barcode"],
        columns=celltype_barcode["cell_type"]
    )

    column_order = (
        (barcode_celltype_table > 0)
        .sum(axis=0)
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    barcode_celltype_table = barcode_celltype_table[column_order]

    row_order = (
        (barcode_celltype_table > 0)
        .sum(axis=1)
        .sort_values(ascending=False)
        .index
        .tolist()
    )
    barcode_celltype_table = barcode_celltype_table.loc[row_order]

    barcode_celltype_table.to_csv(table_path)

    plt.figure(figsize=(6, 12))
    sns.heatmap(
        barcode_celltype_table,
        cmap="coolwarm",
        linewidths=0.2,
        cbar_kws={"shrink": 0.5}
    )
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(rotation=0, ha="right", fontsize=7)
    plt.title(f"Heatmap of Cell Types and Barcode Kinds (pGen < {pgen_cutoff})", pad=20)
    plt.xlabel("cell_type")
    plt.ylabel("barcode")
    plt.grid(visible=False, which="major")
    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-f", "--figure", default="true_barcode_pgen_heatmap.png")
    parser.add_argument("-t", "--table", default="true_barcode_pgen_table.csv")
    parser.add_argument("--pgen-cutoff", type=float, default=1e-6)

    args = parser.parse_args()

    run_true_barcode_pgen_heatmap(
        input_path=args.input,
        figure_path=args.figure,
        table_path=args.table,
        pgen_cutoff=args.pgen_cutoff,
    )
