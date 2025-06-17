# Singlecellvariantanalysis

**scRNA-seq Splicing and Variant Analysis Toolkit**

## Introduction

This repository contains a suite of Python command-line utility scripts for analyzing complex genetic events, such as splice junctions and deletions, from 10x Genomics single-cell RNA sequencing (scRNA-seq) BAM files. While existing tools, such as Vartrix, exist, they are not well-maintained and are not well-equipped to handle large structural variants, such as deletions spanning multiple exons. Vartrix also required a predefined VCF file, which is not always readily available, especially if the exact coordinates of genetic variants are unclear. 

These tools are designed to move from unbiased discovery to quantitative analysis, allowing you to:

1. Discover all novel and annotated splice junctions within a gene region.
2. Quantify the number of cells that contain specific junctions (e.g., from a deletion event).
3. Determine Zygosity by classifying cells as wild-type, heterozygous, or homozygous for a given variant.
4. Check for the presence of a known SNP at a specific site.
5. Filter large BAM files to create smaller, focused BAMs containing only cells of interest for easier visualization and downstream analysis.
   
Table of Contents
Prerequisites
Tools Overview
1. analyze_junctions.py
Usage: discover mode
Usage: quantify mode
2. filter_bam_by_barcodes.py
3. check_snp.py
Example Workflow
License

## Prerequisites
Before using these scripts, you need to have Python 3 and the following libraries installed.

1. **Python 3**: Make sure Python 3 is installed on your system.
2. **pysam library**: This is the core dependency for reading BAM files.
```bash

pip install pysam
```
3. **samtools**: While not required to run the Python scripts, samtools is essential for indexing the BAM files (both input and output), which is required by pysam and viewers like IGV.
```bash

# On Linux (using apt)
sudo apt-get update && sudo apt-get install samtools

# On macOS (using Homebrew)
brew install samtools
```

## Tools Overview
1. analyze_junctions.py - **The Main Analysis Tool**

This is a powerful, multi-purpose tool for junction discovery and quantification. It operates in two modes: discover and quantify.

**Usage: discover mode**

This mode performs an unbiased scan of a genomic region to find all splice junctions that exceed a minimum size. It is perfect for the initial exploration of a gene to find deletions or novel splicing events.

Command:

```bash

python analyze_junctions.py discover --bam <BAM_FILE> --region <REGION> [options...]
```

Arguments:

| Argument	| Description | 	Required |
| ---------|--------------|------------|
--bam	|Path to the indexed input BAM file.|	Yes
--region	|Genomic region, e.g., 'chr17:29032000-29042500'.	|Yes
--barcodes|	(Optional) Path to a text file of cell barcodes to limit the analysis to a subset of cells.|	No
--min-gap-size|	(Optional) The minimum size of a gap to be considered a junction. Default: 200.|	No
--output-file|	(Optional) Name of the output file for the junction report. Default: discovered_junctions.tsv.|	No

Example:

```bash

python analyze_junctions.py discover \
  --bam /data/my_experiment.bam \
  --region "17:29032000-29042500" \
  --min-gap-size 1000 \
  --output-file Srsf3_novel_junctions.tsv
```

This will produce a file Srsf3_novel_junctions.tsv listing all junctions larger than 1000 bp and the number of reads supporting each one.

**Usage: quantify mode**

This mode takes specific junction coordinates that you provide (likely discovered using the discover mode) and performs a detailed analysis to count cells and determine their zygosity.

Command:

```bash

python analyze_junctions.py quantify --bam <BAM_FILE> --region <REGION> --junctions [JUNCTIONS...] [options...]
```

Arguments:

|Argument	|Description	|Required|
|---------|-------------|-------|
--bam	|Path to the indexed input BAM file.	|Yes
--region|	Genomic region, e.g., 'chr17:29032000-29042500'.	|Yes
--junctions|	One or more junctions to quantify. Format: 'NAME:START-END'. Can be specified multiple times.	|Yes
--barcodes|	(Optional) Path to a text file of cell barcodes to limit the analysis to a subset of cells.	|No

Example:

This example quantifies three different alleles (Wild-Type, a deletion of Exon 4, and a deletion of Exons 2+3).

```bash

python analyze_junctions.py quantify \
  --bam /data/my_experiment.bam \
  --region "17:29032000-29042500" \
  --junctions "Wild_Type:29036448-29038489" \
              "Deletion_E4:29038624-29040774" \
              "Deletion_E2_3:29032698-29039781"
```

Output:

The script will print a detailed report directly to the console, like this:
```bash
--- Zygosity Report ---
Total cells expressing any target allele: 2115
----------------------------------------------------------
Count		Allele State(s) Detected in Cell
-----		--------------------------------
1050		Wild_Type
1012		Deletion_E4
33		    Wild_Type + Deletion_E4
20		    Deletion_E2_3
----------------------------------------------------------
```
2. filter_bam_by_barcodes.py **- BAM Filtering Utility**
This is a dedicated utility to create a new, smaller BAM file containing only the reads from a specified list of cells. This is very useful for focused visualization in IGV.

Usage:

```bash

python filter_bam_by_barcodes.py --input-bam <INPUT_BAM> --barcodes <BARCODE_FILE> --output-bam <OUTPUT_BAM>
```

Arguments:
```bash
Argument	Flag	Description	Required
--input-bam	-i	Path to the large, original, indexed BAM file.	Yes
--barcodes	-b	Path to a text file containing one cell barcode per line (e.g., AAACCCAAGCATCAGG-1).	Yes
--output-bam	-o	Path for the new, smaller, filtered BAM file.	Yes
```

Example:

```bash

# First, you might generate a list of barcodes from the analysis tool
# Or use your own list, e.g., 'heterozygous_cells.txt'

python filter_bam_by_barcodes.py \
  -i /data/my_experiment.bam \
  -b heterozygous_cells.txt \
  -o heterozygous_only.bam

# CRITICAL NEXT STEP: Index the new BAM file
samtools index heterozygous_only.bam
```

3. check_snp.py **- SNP Pileup Checker**

This is a simple tool to perform a quick analysis at a single genomic coordinate to check the frequencies of different bases. It is useful for verifying a known SNP but is not a tool for discovering new variants.

Usage:

```bash

python check_snp.py --bam <BAM_FILE> --position <POSITION>
```
Arguments:

|Argument|	Description|	Required|
|--------|-------------|----------|
--bam	|Path to the indexed input BAM file.|	Yes
--position	|A single genomic position to check, e.g., '17:29035150'.|	Yes


Example:

```bash

python check_snp.py \
  --bam /data/my_experiment.bam \
  --position "17:29035150"
```

Example Output:

```bash
--- Checking SNP at 17:29035150 ---

Found 84 reads covering the site.
Base counts from reads:
Base	Count	Percentage
----	-----	----------
A	45	53.57%
G	39	46.43%
```

## Example Workflow
A typical analysis might follow these steps:

1. **Discover**: Use analyze_junctions.py discover on your gene of interest to get a list of all novel deletions and their coordinates.
2. **Quantify**: Use analyze_junctions.py quantify with the coordinates you just discovered (and the wild-type junction) to get a full breakdown of cell populations (wild-type, homozygous, heterozygous).
3. **Filter & Visualize**: Use filter_bam_by_barcodes.py to create a small BAM file containing only the heterozygous cells. Index this new BAM file and load it into IGV for focused visual inspection.
4. **Check a SNP**: If you suspect a SNP is associated with one of your cell populations, use check_snp.py on your filtered BAM to see if the allele frequencies differ.

## License
This project is licensed under the MIT License.
