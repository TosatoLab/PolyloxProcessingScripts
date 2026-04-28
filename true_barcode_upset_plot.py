import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from upsetplot import UpSet


def run_true_barcode_upset_plot(
    input_path,
    figure_path,
    pgen_cutoff=1e-6,
):
    adata = sc.read_h5ad(input_path)

    adata.obs["pGen"] = pd.to_numeric(adata.obs["pGen"], errors="coerce")

    filtered_obs = adata.obs.loc[
        (adata.obs["pGen"] < pgen_cutoff)
        & (~adata.obs["barcode"].isna())
        & (adata.obs["barcode"].astype(str) != "nan"),
        ["barcode", "cell_type"]
    ].copy()

    barcode_celltype_table = pd.crosstab(
        index=filtered_obs["barcode"],
        columns=filtered_obs["cell_type"]
    )

    bool_df = barcode_celltype_table > 0

    if "Doublets" in bool_df.columns:
        bool_df = bool_df.drop(columns=["Doublets"])

    priority_cols = ["Endothelial cells"]
    priority_cols = [c for c in priority_cols if c in bool_df.columns]
    other_cols = [c for c in bool_df.columns if c not in priority_cols]
    bool_df = bool_df[priority_cols + other_cols]

    grouped_data = bool_df.groupby(list(bool_df.columns)).size()

    upset = UpSet(
        grouped_data,
        show_counts=True,
        sort_by="degree",
        orientation="horizontal",
        intersection_plot_elements=10,
    )

    if "Mesenchymal type" in bool_df.columns and "Endothelial cells" in bool_df.columns:
        upset.style_subsets(
            present=["Mesenchymal type", "Endothelial cells"],
            facecolor="blue",
        )
        upset.style_subsets(
            present=["Mesenchymal type"],
            absent=["Endothelial cells"],
            facecolor="orange",
        )
        upset.style_subsets(
            present=["Endothelial cells"],
            absent=["Mesenchymal type"],
            facecolor="red",
        )
    elif "Endothelial cells" in bool_df.columns:
        upset.style_subsets(
            present=["Endothelial cells"],
            facecolor="red",
        )
    elif "Mesenchymal type" in bool_df.columns:
        upset.style_subsets(
            present=["Mesenchymal type"],
            facecolor="orange",
        )

    fig = upset.plot()
    fig["totals"].get_yaxis().set_visible(False)
    plt.subplots_adjust(left=0.2, top=0.95)
    plt.suptitle("")
    plt.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-f", "--figure", default="true_barcode_upset_plot.png")
    parser.add_argument("--pgen-cutoff", type=float, default=1e-6)

    args = parser.parse_args()

    run_true_barcode_upset_plot(
        input_path=args.input,
        figure_path=args.figure,
        pgen_cutoff=args.pgen_cutoff,
    )
