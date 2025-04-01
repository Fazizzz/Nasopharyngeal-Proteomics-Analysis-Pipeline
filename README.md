# Nasopharyngeal Proteomics Analysis Pipeline

This repository contains individual shell scripts and a scalable automated Snakemake-based bioinformatics pipeline for analyzing nasopharyngeal samples using STAR for detecting viral pathogens and generating human alignments for performing differential expression analysis on protein expression data. The goal is to identify % of reads aligned to SARS-Cov-2 and subsiquently any significant proteomic changes between two experimental conditions.

Included are a sample data set that has been analyzed using the provided pipline.

_______________________________________________________
## Overview

The Snakemake workflow automates :

1. Trimming of reads using fastp (Version: 0.23.2)
2. QC of reads using fastqc (Version: 0.11.9)
3. (Optional) Down sampling of reads using Seqtk (Version: 1.3)
4. Alignments of reads using STAR (Version: 2.7.9a)
5. Sorting and indexing of alignments using samtools (Version: 1.12)
6. Generating selective coverage using bedtools (Version: 2.30.0)
7. Generating summary report using multiqc (Version: 1.12)
_______________________________________________________
## Table of Contents
- [Installation](#installation)
- [Configuration](#configuration)
- [Environment Setup](#Environment-Setup)
- [Usage](#usage)
- [Sample Input File Formats](#sample-input-file-formats)
- [Outputs](#outputs)
- [Notes](#notes)
- [License](#license)
- [Acknowledgments](#Acknowledgments)

_______________________________________________________
## Installation
Clone this repository to your local machine:

```bash
git clone https://github.com/Fazizzz/Nasopharyngeal-Proteomics-Analysis-Pipeline.git

cd Nasopharyngeal-Proteomics-Analysis-Pipeline


```
_______________________________________________________
## Configuration

Edit `config.yml` to set your input/output directories and reference paths:

```yaml

input_dir: "/path/to/fastqs"
Index: "/path/to/STAR/index"
GTF: "/path/to/annotation.gtf"
output_dir: "/path/to/output"
log_dir: "/path/to/logs"
DS: 3000000  # Downsampling value, optional
Memory_Max: 20000000000  # Max RAM for sorting BAMs

```
________________________________________________________
## Environment Setup

All required tools are managed using Conda. Create the environment as follows:

```bash
conda env create -f New-snakemake-environment.yml
conda activate snake

```
_______________________________________________________
## Usage

The script can be run by editing the config files with the correct paths. the snakefile and config file should be in the same directory.

Once the config is ready and environment is activated:

```bash

snakemake -s Final_snakemake --use-conda --cores 4

```
You can also run a dry-run to preview jobs:

```bash

snakemake -s Final_snakemake -n

```

________________________________________________________

## Sample Input File Formats

To test the pipeline, provide paired-end FASTQ files in the input_dir specified in your config. Example filenames:

```
sample1_R1.fastq.gz
sample1_R2.fastq.gz
sample2_R1.fastq.gz
sample2_R2.fastq.gz
...

```

_________________________________________________________
## Outputs

The pipeline generates:

1. Quality reports (fastqc, multiqc)
2. Aligned and sorted BAM files
3. Count matrices for Human Genes
5. Logs for each step

___________________________________________________________
## Notes

This pipeline is currently in further development to integrate R scripts for deseq2 analysis of the human genes and will have a subsequent shell script added to identify viral reads. For now it is solely focused on human alignments and producing read counts from an annotataed GTF file in an automated and scalable fashion. Further improvements are in the works and can be requested directly.

 **Downsampling via seqtk is optional and controlled via the DS parameter.

 **STAR indexing must be done in advance.

 **The pipeline assumes paired-end data and GTF annotation format compatible with featureCounts.

___________________________________________________________

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html#license-text) file for details.


___________________________________________________________

## Acknowledgments

* **M.Faizan Khalid** - *Author and current maintainer of the code*

This script was developed by Muhammad Faizan Khalid without input or feedback from the developers of the tools used. For information on the aligner please visit the [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) documentation page. The script is intended to be a utility for use in RNA-Seq analysis. Maintenance and usage of this script is currently not supported or validated by any other developers. For any information on the tools used, please visit their Github and pages to use their guide for citing and referencing their tools.

For citing this tool, please use Khalid M.Faizan or Khalid MF. You can follow my research using my [Google Scholar profile](https://scholar.google.com/citations?hl=en&user=qFZQ5wYAAAAJ&sortby=title&view_op=list_works&gmla=AL3_zigRWGX9g8Jc22idbBUMFuy7cVN_pEIyL6_DXSA-qWkJbcaONzhRNSmAwmQXKEm-3-WYGouZZC2pCE6zD9tZLxizbM7jQzzZMOgtkgsuL825u4lvSs9kwsccajhJbBg2Mrc37at_HCQ).

This project is made possible thanks to the open-source bioinformatics community for their resources and support.
