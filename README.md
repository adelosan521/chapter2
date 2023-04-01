This Github repository contains the scripts used in Chapter 2 of Alejandro De Los Angeles's DPhil thesis. The three bioinformatics tools that were used in the Chapter were DESEQ2, DEXSEQ, and MAJIQ. DESEQ2 is used for gene expression analysis; DEXSEQ is used to assess differential exon usage (DEU); and MAJIQ is used to quantify differential splice junction utilization.

The programming code for the computational pipeline used to pre-process raw sequencing data for DESEQ2 is the "DESEQ2_initial_processing.sh" script. The code for running DESEQ2 is contained in the "DESEQ2_analyses.R" script.

The programming code for the computational pipeline used to pre-process raw sequencing data for DEXSEQ is the "DEXSEQ_initial_processing.sh" script. The code for running DEXSEQ is contained in the "DEXSEQ_analyses.R" script.

Finally, the programming code for running MAJIQ and VOILA is contained in the "majiq.sh" script.

Contact: Alejandro De Los Angeles (adelosan@gmail.com)

<img src="https://raw.githubusercontent.com/PKief/vscode-material-icon-theme/ec559a9f6bfd399b82bb44393651661b08aaf7ba/icons/folder-github-open.svg" width="80" />

## ⚙️ Project Structure

```bash
.
├── DESEQ2_initial_processing.sh
├── DESEQ_analyses.R
├── DEXSEQ_analyses.R
├── DEXSEQ_initial_processing.sh
├── README.md
└── majiq.sh

1 directory, 6 files
```
---
<details closed><summary>.</summary>

| File                         | Summary                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|:-----------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| DESEQ_analyses_R.R           | This code uses the DESeq2 library to analyze gene expression data from two conditions, iPSC and iPSC-neuron. It performs a Principal Component Analysis (PCA) and generates a heatmap to visualize the data.                                                                                                                                                                                                                                                                |
| DESEQ2_initial_processing.sh | This code is a shell script that takes in an argument from the command line and loads modules for trim_galore, fastqc, rna-star, and python3-cbrg.                                                                                                                                                                                                                                                                                                                          |
| DEXSEQ_analyses.sh           | This code uses the DEXSeq library to generate normalized counts for a DEU file and make a table for VGCC subunit DEU and plotting of individual genes.                                                                                                                                                                                                                                                                                                                      |
| majiq.sh                     | This code uses Majiq to build and analyze GFF3 files from Wilfried Haerty, and then uses Voila to view and dump the results to a tsv file.                                                                                                                                                                                                                                                                                                                                  |
| DEXSEQ_initial_processing.sh | This code is a shell script that takes in an argument from the command line and processes it. It loads modules, creates a variable of the file without the _1. fastq and directory path, trims adapters from raw sequencing data FASTQ files using Trim Galore, generates FastQC reports from trimmed adapter files, aligns trimmed adapter files to genome using HISAT2 and generate BAM files, deletes trimmed adapter files after HISAT2, and generates count files from |

</details>
