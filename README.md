# PolyloxProcessingScripts

This repository contains Python scripts used to process and analyze single-cell RNA-seq and Polylox lineage-tracing data associated with the manuscript:

**Single-cell lineage tracing identifies hemogenic endothelial cells in the adult mouse bone marrow**

The scripts support 10x Genomics Cell Ranger processing, Scanpy-based single-cell RNA-seq analysis, Polylox barcode assignment, automated cell-type annotation, doublet detection, cell-cycle analysis, and barcode-sharing visualization.

## Overview

The workflow is organized into four major stages:

1. Cell Ranger processing from FASTQ files.
2. Scanpy preprocessing, quality control, dimensionality reduction, and clustering.
3. Polylox barcode assignment using per-sample Polylox segmentation/assembly files.
4. Downstream annotation and visualization, including true-barcode plots, barcode-sharing heatmaps, UpSet plots, cell-cycle scoring, and marker-expression plots.

These scripts are provided for reproducibility and transparency. They are not intended as a general-purpose software package.

## Repository contents

| File | Purpose |
|---|---|
| `Biowulf_Cellranger.py` | Generates per-sample `cellranger count` shell scripts and a Biowulf/NIH swarm submission file from 10x-style FASTQ files. |
| `ScanpyAnanlysis.py` | Runs a Scanpy workflow for multiple Cell Ranger outputs, including sample loading, QC filtering, normalization, HVG selection, PCA, UMAP, Leiden clustering, marker plotting, and output export. |
| `PolyloxBarcodeAssignment.py` | Assigns Polylox barcodes to cells using sample-specific and mouse-consistent matching. Matching is performed by `sample + mouse + cb`; duplicate calls within the same sample/mouse/cell barcode are resolved by highest `pGen`, then lowest `MinCut`. |
| `doublet_detection.py` | Runs Scrublet-based doublet detection on an AnnData object and saves doublet scores, predicted doublet calls, and a doublet-score distribution plot. |
| `cell_type_annotation.py` | Performs automated cell-type enrichment using `decoupler` and PanglaoDB markers, then stores top cell-type calls in `adata.obs` and rankings in `adata.uns`. |
| `cell_cycle_phase_heatmap.py` | Scores cell-cycle phase using S and G2/M gene sets, stores `S_score`, `G2M_score`, and `phase`, and generates a Leiden-cluster cell-cycle phase heatmap. |
| `true_barcode_umap.py` | Labels cells as true or non-true Polylox barcode cells based on a `pGen` cutoff and generates a UMAP plus a barcode list plot. |
| `true_barcode_pgen_heatmap.py` | Generates a barcode-by-cell-type heatmap for true barcodes and exports the barcode-by-cell-type table. |
| `true_barcode_upset_plot.py` | Generates an UpSet plot showing cell-type sharing among true Polylox barcodes. |
| `linked_barcode_heatmaps.py` | Generates HSPC-linked, endothelial-linked, and mesenchymal-linked barcode-sharing heatmaps and corresponding tables. |
| `stacked_violin_plot.py` | Generates a stacked violin plot for selected endothelial and mesenchymal marker genes across predefined Leiden-cluster groups. |

## Requirements

### System requirements

- Linux or macOS
- Python 3.9 or later recommended
- Cell Ranger installed for `Biowulf_Cellranger.py`
- NIH Biowulf/SLURM environment for swarm submission, if using the generated swarm file

### Python packages

Install the required Python packages with:

```bash
pip install scanpy anndata pandas numpy matplotlib seaborn decoupler upsetplot mygene scrublet
```

Additional notes:

- `cell_type_annotation.py` uses `decoupler` and the PanglaoDB marker resource.
- `cell_cycle_phase_heatmap.py` uses `mygene` to convert human cell-cycle genes to mouse orthologs and therefore requires internet access unless the converted gene lists are already available locally.
- `doublet_detection.py` and the optional `--scrublet` mode in `ScanpyAnanlysis.py` require Scrublet functionality.

## Expected inputs

### Cell Ranger FASTQ input

`Biowulf_Cellranger.py` expects 10x-style FASTQ filenames such as:

```text
SampleName_S1_L001_R1_001.fastq.gz
SampleName_S1_L001_R2_001.fastq.gz
```

A typical project layout is:

```text
Project/
├── FASTQ files or FASTQ subdirectories
├── Ref/
│   └── cellranger_reference/
├── scripts/
├── swarm_logs/
└── cellranger_out/
```

The script scans the project root recursively and skips the following top-level directories:

```text
cellranger_out
scripts
swarm_logs
Ref
```

### Scanpy input

`ScanpyAnanlysis.py` expects Cell Ranger output folders in this structure:

```text
CellRangerCounts/
├── Smp1/
│   └── filtered_feature_bc_matrix/
├── Smp2/
│   └── filtered_feature_bc_matrix/
└── Smp3/
    └── filtered_feature_bc_matrix/
```

### Polylox barcode-assignment input

`PolyloxBarcodeAssignment.py` expects:

1. an AnnData `.h5ad` file with `adata.obs["sample"]`;
2. per-sample Polylox segmentation/assembly CSV files named `<sample>_seg_assemble.csv`;
3. a rare barcode/pGen table;
4. a Polylox barcode library;
5. a MinCut/recombination-support file;
6. optionally, a sample-to-mouse map.

Each segmentation/assembly CSV must contain:

```text
cb
polylox
```

The optional sample-to-mouse map should be a tab-separated file with no header:

```text
Smp1    Mouse#1
Smp2    Mouse#1
SmpA    Mouse#3
```

## Typical workflow

### 1. Generate Cell Ranger scripts on Biowulf

```bash
python Biowulf_Cellranger.py \
  --root /path/to/project \
  --ref /path/to/cellranger/reference
```

This creates:

```text
cellranger_count.swarm
scripts/cellranger_count_<sample>.sh
cellranger_out/
swarm_logs/
```

Submit the swarm job with:

```bash
swarm -f cellranger_count.swarm -t 16 -g 110 --time 48:00:00 --gres=lscratch:200 --logdir swarm_logs
```

Dry run:

```bash
python Biowulf_Cellranger.py \
  --root /path/to/project \
  --ref /path/to/cellranger/reference \
  --dry-run
```

### 2. Run Scanpy preprocessing and clustering

```bash
python ScanpyAnanlysis.py \
  --input-root /path/to/CellRangerCounts \
  --samples Smp1 Smp2 Smp3 Smp4 Smp5 Smp6 \
  --outdir /path/to/scanpy_out
```

Optional Scrublet filtering during preprocessing:

```bash
python ScanpyAnanlysis.py \
  --input-root /path/to/CellRangerCounts \
  --samples Smp1 Smp2 Smp3 Smp4 Smp5 Smp6 \
  --outdir /path/to/scanpy_out \
  --scrublet
```

Main outputs include:

```text
adata_merged_rawcounts.h5ad
adata_processed.h5ad
adata_obs_processed.tsv
adata_var_processed.tsv
cluster_cell_counts.tsv
run_parameters.json
pipeline.log
marker_umaps/
```

### 3. Assign Polylox barcodes

```bash
python PolyloxBarcodeAssignment.py \
  --adata /path/to/scanpy_out/adata_processed.h5ad \
  --output /path/to/scanpy_out/adata_polylox_annotated.h5ad \
  --base-path /path/to/CSV_files_Snakemake/ \
  --rare-barcode-list /path/to/merged_polylox_data_pgen.txt \
  --barcode-lib /path/to/polylox_barcodelib.txt \
  --min-recomb /path/to/min_reconb_list.txt \
  --sample-mouse-map /path/to/sample_mouse.tsv
```

The script adds the following fields to `adata.obs`:

```text
mouse
cb
polylox
pGen
MinCut
barcode
```

Important matching behavior:

- Matching is performed using `sample + mouse + cb`.
- Barcodes are not matched across different samples.
- Cells from samples missing in the sample-to-mouse map are marked as `NoAssign` and do not receive Polylox assignments.
- Duplicate assignments within the same `sample + mouse + cb` are resolved by highest `pGen`, then lowest `MinCut`.
- The final `barcode` field is the reverse complement of `polylox` unless `--no-reverse-complement` is used.

### 4. Optional standalone doublet detection

```bash
python doublet_detection.py \
  -i /path/to/adata_polylox_annotated.h5ad \
  -o /path/to/adata_doublet_checked.h5ad \
  -f doublet_score_distribution.png
```

This adds:

```text
doublet_score
predicted_doublet
```

to `adata.obs`.

### 5. Automated cell-type annotation

```bash
python cell_type_annotation.py \
  -i /path/to/adata_doublet_checked.h5ad \
  -o /path/to/adata_celltype_annotated.h5ad \
  --groupby leiden
```

This adds:

```text
ora_cell_type
ora_top_celltypes
```

to `adata.obs`, and stores ranking outputs in:

```text
adata.uns["ora_cell_type_rankings"]
adata.uns["ora_top_celltypes_by_cluster"]
```

### 6. Cell-cycle scoring

```bash
python cell_cycle_phase_heatmap.py \
  -i /path/to/adata_celltype_annotated.h5ad \
  -o /path/to/adata_cell_cycle.h5ad \
  -f cell_cycle_phase_heatmap.png
```

This adds:

```text
S_score
G2M_score
phase
```

to `adata.obs`.

### 7. True-barcode UMAP and barcode list

```bash
python true_barcode_umap.py \
  -i /path/to/adata_cell_cycle.h5ad \
  -o /path/to/adata_true_bc.h5ad \
  --umap-figure true_barcode_umap.png \
  --barcode-list-figure true_barcode_list.png
```

Default true-barcode threshold:

```text
pGen < 1e-6
```

This adds:

```text
True BC
```

to `adata.obs`.

### 8. True-barcode pGen heatmap

```bash
python true_barcode_pgen_heatmap.py \
  -i /path/to/adata_true_bc.h5ad \
  -f true_barcode_pgen_heatmap.png \
  -t true_barcode_pgen_table.csv
```

Default threshold:

```text
pGen < 1e-6
```

### 9. True-barcode UpSet plot

```bash
python true_barcode_upset_plot.py \
  -i /path/to/adata_true_bc.h5ad \
  -f true_barcode_upset_plot.png
```

Default threshold:

```text
pGen < 1e-6
```

### 10. Linked-barcode heatmaps

```bash
python linked_barcode_heatmaps.py \
  -i /path/to/adata_true_bc.h5ad \
  -o linked_barcode_heatmaps
```

This creates:

```text
linked_barcode_heatmaps/
├── hspc_linked_barcode_heatmap.png
├── endothelial_linked_barcode_heatmap.png
├── mesenchymal_linked_barcode_heatmap.png
├── hspc_linked_barcode_table.csv
├── endothelial_linked_barcode_table.csv
└── mesenchymal_linked_barcode_table.csv
```

Default thresholds:

```text
HSPC-linked heatmap: pGen < 1e-6
Endothelial-linked heatmap: pGen < 1e-6
Mesenchymal-linked heatmap: pGen < 1e-6
```

### 11. Stacked violin marker plot

```bash
python stacked_violin_plot.py \
  -i /path/to/adata_true_bc.h5ad \
  -o /path/to/adata_group_labels.h5ad \
  -f stacked_violin_plot.png
```

The default marker genes are:

```text
Cdh5
Pecam1
Kdr
Flt1
Cxcl12
Pdgfrb
Lepr
Runx1
Col1a2
```

The script groups the following Leiden clusters:

```text
0, 1, 13, 22 -> Endothelial cells
14 -> Mesenchymal type
```

These cluster labels are dataset-specific and should be updated if the clustering result changes.

## Important parameters and defaults

### Scanpy preprocessing

| Parameter | Default | Meaning |
|---|---:|---|
| `--min-genes` | `200` | Minimum genes detected per cell. |
| `--min-counts` | `500` | Minimum UMI counts per cell. |
| `--min-cells` | `3` | Minimum cells expressing a gene. |
| `--max-genes-by-counts` | `8000` | Upper gene-count filter for cells. |
| `--max-pct-mt` | `10.0` | Maximum mitochondrial percentage. |
| `--target-sum` | `10000` | Library-size normalization target. |
| `--n-neighbors` | `10` | Neighbor graph parameter. |
| `--n-pcs` | `40` | Number of PCs used for neighbor graph. |
| `--resolution` | `0.8` | Leiden clustering resolution. |

### Polylox barcode assignment

| Parameter | Default | Meaning |
|---|---:|---|
| `--cb-length` | `16` | Number of nucleotides used to extract the 10x cell barcode from `adata.obs_names`. |
| `--missing-barcode-pgen` | `1e-8` | pGen assigned to barcodes absent from the barcode library. |
| `--default-mincut` | `2` | MinCut value used when MinCut is missing. |
| `--no-reverse-complement` | off | Stores `polylox` directly as `barcode` instead of reverse complementing it. |

### Barcode visualization

| Script | Default threshold |
|---|---:|
| `true_barcode_umap.py` | `pGen < 1e-6` |
| `true_barcode_pgen_heatmap.py` | `pGen < 1e-6` |
| `true_barcode_upset_plot.py` | `pGen < 1e-6` |
| `linked_barcode_heatmaps.py`, HSPC-linked plot | `pGen < 1e-6` |
| `linked_barcode_heatmaps.py`, endothelial/mesenchymal-linked plots | `pGen < 1e-6` |

## Reproducibility notes

- Scripts are designed to run non-interactively from the command line.
- Most thresholds are exposed as command-line parameters and should be reported when used in a manuscript or data descriptor.
- The main Scanpy pipeline writes `run_parameters.json` to preserve run settings.
- Cell-type annotations generated by `cell_type_annotation.py` are automated marker-enrichment calls and should be biologically reviewed before manuscript-level interpretation.
- The stacked violin cluster grouping in `stacked_violin_plot.py` is dataset-specific and may need revision if Leiden cluster numbering changes.
- Polylox barcode assignment is sample- and mouse-aware. The script intentionally avoids assigning a Polylox barcode from one sample to a cell from another sample.

## Suggested complete command sequence

```bash
python Biowulf_Cellranger.py \
  --root /path/to/project \
  --ref /path/to/cellranger/reference

swarm -f cellranger_count.swarm -t 16 -g 110 --time 48:00:00 --gres=lscratch:200 --logdir swarm_logs

python ScanpyAnanlysis.py \
  --input-root /path/to/CellRangerCounts \
  --samples Smp1 Smp2 Smp3 Smp4 Smp5 Smp6 \
  --outdir /path/to/scanpy_out

python PolyloxBarcodeAssignment.py \
  --adata /path/to/scanpy_out/adata_processed.h5ad \
  --output /path/to/scanpy_out/adata_polylox_annotated.h5ad \
  --base-path /path/to/CSV_files_Snakemake/ \
  --rare-barcode-list /path/to/merged_polylox_data_pgen.txt \
  --barcode-lib /path/to/polylox_barcodelib.txt \
  --min-recomb /path/to/min_reconb_list.txt \
  --sample-mouse-map /path/to/sample_mouse.tsv

python doublet_detection.py \
  -i /path/to/scanpy_out/adata_polylox_annotated.h5ad \
  -o /path/to/scanpy_out/adata_doublet_checked.h5ad \
  -f doublet_score_distribution.png

python cell_type_annotation.py \
  -i /path/to/scanpy_out/adata_doublet_checked.h5ad \
  -o /path/to/scanpy_out/adata_celltype_annotated.h5ad \
  --groupby leiden

python cell_cycle_phase_heatmap.py \
  -i /path/to/scanpy_out/adata_celltype_annotated.h5ad \
  -o /path/to/scanpy_out/adata_cell_cycle.h5ad \
  -f cell_cycle_phase_heatmap.png

python true_barcode_umap.py \
  -i /path/to/scanpy_out/adata_cell_cycle.h5ad \
  -o /path/to/scanpy_out/adata_true_bc.h5ad \
  --umap-figure true_barcode_umap.png \
  --barcode-list-figure true_barcode_list.png

python true_barcode_pgen_heatmap.py \
  -i /path/to/scanpy_out/adata_true_bc.h5ad \
  -f true_barcode_pgen_heatmap.png \
  -t true_barcode_pgen_table.csv

python true_barcode_upset_plot.py \
  -i /path/to/scanpy_out/adata_true_bc.h5ad \
  -f true_barcode_upset_plot.png

python linked_barcode_heatmaps.py \
  -i /path/to/scanpy_out/adata_true_bc.h5ad \
  -o linked_barcode_heatmaps

python stacked_violin_plot.py \
  -i /path/to/scanpy_out/adata_true_bc.h5ad \
  -o /path/to/scanpy_out/adata_group_labels.h5ad \
  -f stacked_violin_plot.png
```

## Citation

If these scripts are used, please cite the associated manuscript:

**Single-cell lineage tracing identifies hemogenic endothelial cells in the adult mouse bone marrow**

## License

This repository is released under the MIT License. See `LICENSE.md` for details.
