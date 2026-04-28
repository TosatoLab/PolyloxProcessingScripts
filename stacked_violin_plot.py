import argparse
import scanpy as sc
import matplotlib.pyplot as plt


def run_stacked_violin_plot(
    input_path,
    output_path,
    figure_path,
):
    adata = sc.read_h5ad(input_path)

    genes = ["Cdh5", "Pecam1", "Kdr", "Flt1", "Cxcl12", "Pdgfrb", "Lepr", "Runx1", "Col1a2"]

    group_names = {
        "0": "Endothelial cells",
        "1": "Endothelial cells",
        "13": "Endothelial cells",
        "22": "Endothelial cells",
        "14": "Mesenchymal type",
    }

    adata.obs["group_labels"] = adata.obs["leiden"].astype(str).map(group_names)

    target_clusters = list(group_names.keys())
    adata_subset = adata[adata.obs["leiden"].astype(str).isin(target_clusters)].copy()

    plt.figure(dpi=300)
    sc.pl.stacked_violin(
        adata_subset,
        var_names=genes,
        groupby="group_labels",
        use_raw=True,
        figsize=(6, 2),
        dendrogram=False,
        swap_axes=False,
        stripplot=False,
        show=False,
    )
    plt.subplots_adjust(right=1.5)
    plt.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close()

    adata.write_h5ad(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-f", "--figure", default="stacked_violin_plot.png")

    args = parser.parse_args()

    run_stacked_violin_plot(
        input_path=args.input,
        output_path=args.output,
        figure_path=args.figure,
    )
