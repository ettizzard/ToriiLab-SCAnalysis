## The Problem

Given a set of Illumina RNA-Seq data from an *A. thaliana* sample:

1) Filter the low quality cells (low UMI transcripts per cell).

2) Perform UMAP analysis and clustering.

3) Report any interesting findings.

## Prerequisites

The dataset can be found [here](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR34735789&display=metadata) at the NCBI site.

In order to download the data locally, the [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki) was used due to the dataset's large size.

After successful download, the dataset was stored locally with the *"SRR34735789.sra"* filename.

Filtering and UMAP analysis will be carried out with [Seurat](https://satijalab.org/seurat/), so the file was converted to .fastq with SRA Toolkit before proceeding.

Once that was complete, 