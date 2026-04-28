import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import gridspec
from matplotlib.colorbar import ColorbarBase


def make_table(adata, pgen_cutoff):
    obs = adata.obs.copy()
    obs["pGen"] = pd.to_numeric(obs["pGen"], errors="coerce")
    obs = obs.loc[obs["pGen"] < pgen_cutoff, ["barcode", "cell_type"]].dropna()
    return pd.crosstab(index=obs["barcode"], columns=obs["cell_type"])


def order_heatmap_data(table, target_cell_type, require_other_link):
    if target_cell_type not in table.columns:
        return pd.DataFrame()

    linked_table = table.loc[table[target_cell_type] > 0].copy()

    if require_other_link:
        other_columns = [col for col in linked_table.columns if col != target_cell_type]
        if other_columns:
            linked_table = linked_table.loc[(linked_table[other_columns] > 0).sum(axis=1) > 0]

    if linked_table.empty:
        return pd.DataFrame()

    heatmap_data = linked_table.T
    row_totals = heatmap_data.sum(axis=1).sort_values(ascending=False)
    row_order = [target_cell_type] + [idx for idx in row_totals.index if idx != target_cell_type]
    heatmap_data = heatmap_data.loc[row_order]
    col_order = heatmap_data.sum(axis=0).sort_values(ascending=False).index
    heatmap_data = heatmap_data.loc[:, col_order]

    return heatmap_data


def make_cmap():
    colors = [
        "grey",
        "lightblue",
        "blue",
        "magenta",
        "lightgreen",
        "green",
        "orange",
        "red",
    ]
    bounds = [0, 1, 2, 3, 4, 5, 6, 7, 100]
    return ListedColormap(colors), BoundaryNorm(bounds, len(colors)), bounds


def plot_simple_heatmap(heatmap_data, figure_path, title, figsize):
    cmap, norm, bounds = make_cmap()

    if heatmap_data.empty:
        plt.figure(figsize=figsize, dpi=300)
        plt.title(title)
        plt.text(0.5, 0.5, "No matching barcodes", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(figure_path, dpi=300)
        plt.close()
        return

    plt.figure(figsize=figsize, dpi=300)
    sns.heatmap(
        heatmap_data,
        cmap=cmap,
        norm=norm,
        linewidths=0.2,
        cbar_kws={"shrink": 0.5, "ticks": list(range(9))},
        annot=True,
        fmt="d",
        annot_kws={"size": 7},
    )
    plt.title(title, pad=20)
    plt.xticks(rotation=90, ha="center", fontsize=8)
    plt.yticks(rotation=0, ha="right", fontsize=10)
    plt.xlabel("barcode")
    plt.ylabel("cell_type")
    plt.grid(visible=False, which="major")
    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()


def plot_heatmap_with_totals(heatmap_data, figure_path, title, figsize):
    cmap, norm, bounds = make_cmap()

    if heatmap_data.empty:
        plt.figure(figsize=figsize, dpi=300)
        plt.title(title)
        plt.text(0.5, 0.5, "No matching barcodes", ha="center", va="center")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(figure_path, dpi=300)
        plt.close()
        return

    row_totals = heatmap_data.sum(axis=1).astype(int)

    fig = plt.figure(figsize=figsize, dpi=300)
    gs = gridspec.GridSpec(1, 3, width_ratios=[20, 1.5, 0.6], wspace=0.05)

    ax0 = plt.subplot(gs[0])
    sns.heatmap(
        heatmap_data,
        ax=ax0,
        cmap=cmap,
        norm=norm,
        linewidths=0.2,
        cbar=False,
        annot=True,
        fmt="d",
        annot_kws={"size": 8},
    )
    ax0.set_ylabel("cell_type", fontsize=10)
    ax0.set_xlabel("barcode", fontsize=10)
    ax0.set_title(title, pad=20)
    ax0.set_yticklabels(heatmap_data.index, rotation=0, fontsize=8)
    ax0.set_xticklabels(ax0.get_xticklabels(), rotation=90, fontsize=8)

    ax1 = plt.subplot(gs[1], sharey=ax0)
    ax1.axis("off")
    for i, value in enumerate(row_totals):
        ax1.text(0, i + 0.5, str(value), va="center", ha="left", fontsize=8)

    ax2 = plt.subplot(gs[2])
    ColorbarBase(
        ax2,
        cmap=cmap,
        norm=norm,
        boundaries=bounds,
        ticks=range(9),
        orientation="vertical",
    )
    ax2.set_title("Count", fontsize=10)

    plt.tight_layout()
    plt.savefig(figure_path, dpi=300)
    plt.close()


def run_linked_barcode_heatmaps(input_path, output_dir, hspc_pgen_cutoff, linked_pgen_cutoff):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(input_path)

    hspc_table = make_table(adata, hspc_pgen_cutoff)
    linked_table = make_table(adata, linked_pgen_cutoff)

    hspc_heatmap = order_heatmap_data(hspc_table, "HSPC", False)
    endothelial_heatmap = order_heatmap_data(linked_table, "Endothelial cells", True)
    mesenchymal_heatmap = order_heatmap_data(linked_table, "Mesenchymal type", True)

    hspc_heatmap.to_csv(output_dir / "hspc_linked_barcode_table.csv")
    endothelial_heatmap.to_csv(output_dir / "endothelial_linked_barcode_table.csv")
    mesenchymal_heatmap.to_csv(output_dir / "mesenchymal_linked_barcode_table.csv")

    plot_simple_heatmap(
        hspc_heatmap,
        output_dir / "hspc_linked_barcode_heatmap.png",
        "HSPC-Linked Cell Types per Barcode",
        (12, 6),
    )

    plot_heatmap_with_totals(
        endothelial_heatmap,
        output_dir / "endothelial_linked_barcode_heatmap.png",
        "Endothelial cells-Linked Cell Types per Barcode",
        (12, 3.5),
    )

    plot_heatmap_with_totals(
        mesenchymal_heatmap,
        output_dir / "mesenchymal_linked_barcode_heatmap.png",
        "Mesenchymal-Linked Cell Types per Barcode",
        (4, 4),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output-dir", default="linked_barcode_heatmaps")
    parser.add_argument("--hspc-pgen-cutoff", type=float, default=1e-6)
    parser.add_argument("--linked-pgen-cutoff", type=float, default=1e-6)

    args = parser.parse_args()

    run_linked_barcode_heatmaps(
        input_path=args.input,
        output_dir=args.output_dir,
        hspc_pgen_cutoff=args.hspc_pgen_cutoff,
        linked_pgen_cutoff=args.linked_pgen_cutoff,
    )
