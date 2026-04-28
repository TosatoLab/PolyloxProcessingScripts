import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import mygene


def convert_human_to_mouse_genes(human_genes):
    mg = mygene.MyGeneInfo()
    query_result = mg.querymany(human_genes, scopes="symbol", fields="symbol", species="mouse")
    mouse_genes = []
    for entry in query_result:
        if "symbol" in entry:
            mouse_genes.append(entry["symbol"])
    return list(dict.fromkeys(mouse_genes))


def plot_phase_heatmap(phase_props, figure_path):
    fig, ax = plt.subplots(figsize=(6, 12), dpi=300)
    im = ax.imshow(phase_props.values, aspect="auto")

    ax.set_xticks(range(len(phase_props.columns)))
    ax.set_xticklabels(phase_props.columns)
    ax.set_yticks(range(len(phase_props.index)))
    ax.set_yticklabels(phase_props.index)
    ax.set_xlabel("Cell Cycle Phase")
    ax.set_ylabel("Leiden Cluster")
    ax.set_title("Cell Cycle Phase Proportions per Leiden Cluster")

    for i in range(phase_props.shape[0]):
        for j in range(phase_props.shape[1]):
            ax.text(j, i, f"{phase_props.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Fraction")

    plt.tight_layout()
    plt.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close()


def run_cell_cycle_scoring(
    input_path,
    output_path,
    figure_path,
):
    adata = sc.read_h5ad(input_path)

    human_s_genes = [
        "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
        "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "RPA2", "NASP",
        "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
        "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
        "USP1", "CLSPN", "POLA1", "CHAF1B", "MRPL36", "E2F8"
    ]
    human_g2m_genes = [
        "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
        "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
        "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
        "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
        "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
        "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"
    ]

    mouse_s_genes = convert_human_to_mouse_genes(human_s_genes)
    mouse_g2m_genes = convert_human_to_mouse_genes(human_g2m_genes)

    if adata.raw is not None:
        adata_cc = sc.AnnData(X=adata.raw.X, obs=adata.obs.copy(), var=adata.raw.var.copy())
    else:
        adata_cc = adata.copy()

    s_genes = [g for g in mouse_s_genes if g in adata_cc.var_names]
    g2m_genes = [g for g in mouse_g2m_genes if g in adata_cc.var_names]

    sc.tl.score_genes_cell_cycle(adata_cc, s_genes=s_genes, g2m_genes=g2m_genes)

    adata.obs["S_score"] = adata_cc.obs["S_score"]
    adata.obs["G2M_score"] = adata_cc.obs["G2M_score"]
    adata.obs["phase"] = adata_cc.obs["phase"]

    phase_counts = pd.crosstab(adata.obs["leiden"], adata.obs["phase"])
    phase_props = phase_counts.div(phase_counts.sum(axis=1), axis=0)

    plot_phase_heatmap(phase_props, figure_path)

    adata.write_h5ad(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-f", "--figure", default="cell_cycle_phase_heatmap.png")

    args = parser.parse_args()

    run_cell_cycle_scoring(
        input_path=args.input,
        output_path=args.output,
        figure_path=args.figure,
    )
