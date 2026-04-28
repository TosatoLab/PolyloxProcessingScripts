import argparse
import math
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def plot_barcode_list_5cols_single_page(
    adata_in,
    barcode_col="barcode",
    group_col="True BC",
    palette=None,
    sort_barcodes=True,
    ncols=5,
    figsize=(6, 10),
    fontsize=8,
    box_size=0.010,
    xpad=0.004,
    top=0.95,
    bottom=0.05,
    left=0.03,
    right=0.98,
):
    if palette is None:
        palette = {"Not True barcodes": "tab:blue", "True barcodes": "tab:orange"}

    df = adata_in.obs[[barcode_col, group_col]].copy()
    df = df.dropna(subset=[barcode_col, group_col])
    df[barcode_col] = df[barcode_col].astype(str)

    bc_to_group = df.groupby(barcode_col)[group_col].apply(
        lambda s: "True barcodes" if (s == "True barcodes").any() else "Not True barcodes"
    )

    items = list(bc_to_group.items())
    if sort_barcodes:
        items = sorted(items, key=lambda x: x[0])

    if len(items) == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No barcodes found", ha="center", va="center", transform=ax.transAxes)
        plt.tight_layout()
        return fig, ax

    rows_per_col = int(math.ceil(len(items) / ncols))
    col_width = (right - left) / ncols
    col_x = [left + i * col_width for i in range(ncols)]
    usable_h = top - bottom
    row_step = usable_h / (rows_per_col + 1)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_axis_off()

    for i, (bc, g) in enumerate(items):
        col = i // rows_per_col
        row = i % rows_per_col
        if col >= ncols:
            break

        x0 = col_x[col]
        y0 = top - (row + 1) * row_step

        box_w = box_size * 1.2
        box_h = box_size * 0.6

        rect = mpatches.Rectangle(
            (x0, y0 - box_h * 0.6),
            box_w,
            box_h,
            transform=ax.transAxes,
            facecolor=palette.get(g, "gray"),
            edgecolor="black",
            linewidth=0.3,
        )
        ax.add_patch(rect)

        ax.text(
            x0 + box_size + xpad,
            y0,
            bc,
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=fontsize,
        )

    handles = [
        mpatches.Patch(facecolor=palette["True barcodes"], edgecolor="black", label="True barcodes"),
        mpatches.Patch(facecolor=palette["Not True barcodes"], edgecolor="black", label="Not True barcodes"),
    ]
    ax.legend(handles=handles, loc="lower center", ncol=2, frameon=False)
    plt.tight_layout()
    return fig, ax


def run_true_barcode_umap(
    input_path,
    output_path,
    umap_figure_path,
    barcode_list_figure_path,
    pgen_cutoff=1e-6,
):
    adata = sc.read_h5ad(input_path)

    adata.obs["pGen"] = pd.to_numeric(adata.obs["pGen"], errors="coerce")
    adata.obs["True BC"] = "Not True barcodes"
    adata.obs.loc[adata.obs["pGen"] < pgen_cutoff, "True BC"] = "True barcodes"

    adata_subset = adata[
        (adata.obs["cell_type"] != "Doublets")
        & (adata.obs["barcode"].astype(str) != "nan")
        & (~adata.obs["barcode"].isna())
    ].copy()

    fig, ax = plt.subplots(figsize=(5, 5))
    sc.pl.umap(
        adata_subset,
        color=["True BC"],
        ax=ax,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(umap_figure_path, dpi=300)
    plt.close()

    palette = {"Not True barcodes": "tab:blue", "True barcodes": "tab:orange"}
    fig, ax = plot_barcode_list_5cols_single_page(
        adata_subset,
        barcode_col="barcode",
        group_col="True BC",
        palette=palette,
        ncols=5,
        figsize=(6, 10),
        fontsize=8,
        box_size=0.010,
    )
    fig.savefig(barcode_list_figure_path, dpi=300)
    plt.close(fig)

    adata.write_h5ad(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--umap-figure", default="true_barcode_umap.png")
    parser.add_argument("--barcode-list-figure", default="true_barcode_list.png")
    parser.add_argument("--pgen-cutoff", type=float, default=1e-6)

    args = parser.parse_args()

    run_true_barcode_umap(
        input_path=args.input,
        output_path=args.output,
        umap_figure_path=args.umap_figure,
        barcode_list_figure_path=args.barcode_list_figure,
        pgen_cutoff=args.pgen_cutoff,
    )
