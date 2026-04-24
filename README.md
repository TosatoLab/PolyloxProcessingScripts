# PolyloxProcessingScripts

This repository contains Python scripts used for processing and analyzing single-cell Polylox lineage-tracing data associated with the manuscript:

**Single-cell lineage tracing identifies hemogenic endothelial cells in the adult mouse bone marrow**

The scripts were used to support preprocessing of 10x Genomics single-cell RNA-seq data, Polylox barcode assignment, and downstream Scanpy-based analysis.

## Repository contents

| File | Description |
|---|---|
| `Biowulf_Cellranger.py` | Generates per-sample `cellranger count` shell scripts and a Biowulf/NIH swarm submission file from 10x-style FASTQ files. |
| `ScanpyAnanlysis.py` | Runs a Scanpy workflow for multiple Cell Ranger outputs, including sample loading, QC filtering, normalization, highly variable gene selection, PCA, UMAP, Leiden clustering, marker visualization, and output export. |
| `PolyloxBarcodeAssignment.py` | Assigns Polylox barcodes to cells in an AnnData object by matching cell barcodes to per-sample Polylox segmentation/assembly CSV files. |

## Purpose

These scripts are provided to document the computational procedures used for Polylox single-cell lineage-tracing analysis. They are intended primarily for reproducibility and transparency rather than as a general-purpose software package.

## Requirements

The scripts require Python 3 and the following Python packages:

```bash
pandas
scanpy
anndata
matplotlib
