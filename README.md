# Deliverables

Given a set of Illumina snRNA-Seq data from an *A. thaliana* sample:

1) Filter the low quality cells (low UMI transcripts per cell).

2) Perform UMAP analysis and clustering.

3) Report any interesting findings.


## Prerequisites

The provided dataset can be found [here](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR34735789&display=metadata) at the NCBI site. More information about the origin and generation of this data can be found [here](https://www.ncbi.nlm.nih.gov/sra/SRX29870010). Notably, the snRNA-Seq library was developed with the [10X Genomics Single Cell 3' Library & Gel Bead Kit v3.1 protocol](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf) and was subsequently sequenced on an Illumina NovaSeq 6000. 

In order to download the data locally, the [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki) was used due to the dataset's large size.

After successful download, the dataset was stored locally with the *"SRR34735789.sra"* filename. The file was then converted to `.fastq` and split into Read_1 and Read_2 files with SRA Toolkit before proceeding.

`$ fastq-dump --split-files /Users/evan/Desktop/sra_cache/sra/SRR34735789.sra`


## Step 1: Filtering

Before any filtering can be accomplished, the reads must be mapped to a reference genome, then to genes, and finally to individual cells/nuclei.

A toplevel, unmasked DNA reference genome assembly was acquired from EnsemblPlants [here](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz). According to a [10X Genomics Article](https://kb.10xgenomics.com/s/article/360060307872-Does-my-genome-sequence-needs-be-masked-or-unmasked-for-custom-reference-generation), the primary/toplelvel DNA assembly is preferable for 10X 3' assay analysis.

A corresponding genetic feature format file (`.gff3`) was downloaded from EnsemblPlants [here](https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz).

Although CellRanger is a first-party and perfectly suitable option for mapping and filtering these reads, I opted for [STAR](https://github.com/alexdobin/STAR), specifically utilizing STARsolo, as I could run it on my local machine. Additionally, STARsolo touts faster speed and equivalent quality to CellRanger.


#### Genome Index Generation:
`$ STAR \`
`--runThreadN 8 \`
`--runMode genomeGenerate \`
`--genomeDir /Users/evan/bioinfo/ToriiLab-SCAnalysis/genome_indices/ \`
`--genomeFastaFiles /Users/evan/bioinfo/ToriiLab-SCAnalysis/reference_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \`
`--sjdbGTFfile /Users/evan/bioinfo/ToriiLab-SCAnalysis/gene_annotation/Arabidopsis_thaliana.TAIR10.62.gff3 \`
`--sjdbGTFtagExonParentTranscript Parent \`
`--genomeSAindexNbases 12`

Notably, the `--sjdbGTFtagExonParentTranscript Parent` is necessary here due to the gene annotation file's `.gff3` format. Additionally, `--genomeSAindexNbases` had to be scaled down from the default value of 14 to 12 due to the relatively small genome size. After successful genome index generation, we can proceed with mapping/alignment.


#### Mapping and Filtering:

While trying to map, I kept encountering a fatal error with the bioconda version of STAR. I had to recompile a newer [alpha version](https://github.com/dobinlab/STAR_pre_releases/releases/tag/2.7.11b_alpha_2024-02-09) and then modify the source makefile to get this working, [as suggested by a fellow researcher on GitHub](https://github.com/alexdobin/STAR/issues/2142).

`$ STAR \`  
`--runThreadN 8 \`  
`--runMode alignReads \`  
`--genomeDir /Users/evan/bioinfo/ToriiLab-SCAnalysis/genome_indices/ \`  
`--readFilesIn /Users/evan/bioinfo/ToriiLab-SCAnalysis/snRNA-Seq_reads/SRR34735789_2.fastq /Users/evan/bioinfo/ToriiLab-SCAnalysis/snRNA-Seq_reads/SRR34735789_1.fastq \`  
`--sjdbGTFfile /Users/evan/bioinfo/ToriiLab-SCAnalysis/gene_annotation/Arabidopsis_thaliana.TAIR10.62.gff3 \`  
`--sjdbGTFtagExonParentTranscript Parent \`  
`--soloUMIlen 12 \`  
`--soloType CB_UMI_Simple \`  
`--soloCBwhitelist /Users/evan/bioinfo/ToriiLab-SCAnalysis/cellbarcode_whitelist/3M-february-2018_TRU.txt \`  
`--soloFeatures GeneFull`  










Notably, Read_1 contains barcode sequences and Read_2 contains the cDNA sequences. The additional parameters `--soloUMIlen 12`, `--soloType CB_UMI_Simple`, `--soloCBwhitelist /Users/evan/bioinfo/ToriiLab-SCAnalysis/cellbarcode_whitelist/3M-february-2018_TRU.txt`, and `--soloFeatures GeneFull` were chosen based on [this guide](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) for running STARsolo on snRNA-Seq 10X Chromium V3 data.

##### Current Issue

The "features.tsv" file is persistently blank after succesful mapping runs, likely indicating an issue with the `.gff3` annotation file. Presently working on a conversion/formatting script before proceeding with downstream visualization via Seurat.

<!-- Additional filtering parameters were chosen from [this STARsolo publication](https://www.biorxiv.org/content/10.1101/2021.05.05.442755v1.full.pdf) in an effort to most closely mirror CellRanger's filtering behavior within STAR. Namely, `` -->







<!-- Filtering and UMAP analysis will be carried out with [Seurat](https://satijalab.org/seurat/), so a fresh version of Seurat and its additional packages were installed within R Studio [as per the developer recommendations](https://satijalab.org/seurat/articles/install_v5).

Before diving into Seurat, however, unique molecular identified (UMI) count matrices are needed. Since this dataset was generated via 10X Genomics materials, CellRanger would be the ideal software to generate matrices. However, I chose to complete this analysis locally, so I opted for STARsolo instead. -->

