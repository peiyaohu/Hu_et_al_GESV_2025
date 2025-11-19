## GESV: A Simple Pipeline for Amplicon Sequencing-Based Genotyping of CRISPR Edits

### **Project Overview: GESV Pipeline**

This repository contains the complete source code, demo data, and analysis scripts for the GESV (gene edited sequence variants) pipeline, as presented in our publication. GESV is a robust bioinformatics workflow designed for accurate **identification and genotyping of sequence variants (SVs)** from high-throughput amplicon sequencing data, with a primary application in analyzing CRISPR gene-edited populations.

The project is organized into the following main directories, each representing a key component of our study:

------

### **Directory Structure and Description**

**1. [`1_strategy_comparison`](https://github.com/peiyaohu/Hu_et_al_GESV_2025/tree/main/1_strategy_comparison)**

- **Purpose**: To systematically **evaluate and compare different amplicon sequencing strategies** for their effectiveness in generating and detecting sequence variants.
- **Contents**: Includes demo datasets and scripts for SV generations.

**2. [`2_GESV_on_LRR-RLK`](https://github.com/peiyaohu/Hu_et_al_GESV_2025/tree/main/2_GESV_on_LRR-RLK)**

- **Purpose**: Demonstrates the application of the GESV pipeline to a large-scale dataset. It was applied to amplicon sequencing data from **212 gene-edited rice lines** derived from our previously published [receptor-like kinase (LRR-RLK) CRISPR library](https://www.sciencedirect.com/science/article/pii/S1674205221004056).
- **Key Scripts**:
  - [`GESV_LRR.sh`](https://github.com/peiyaohu/Hu_et_al_GESV_2025/blob/main/2_GESV_on_LRR-RLK/GESV_LRR.sh): The main Bash script that executes the core GESV pipeline to generate sequence variants, resulting to File1.fa mentioned in the Figure1C.
  - [`align.R`](https://github.com/peiyaohu/Hu_et_al_GESV_2025/blob/main/2_GESV_on_LRR-RLK/align.R): A R script responsible for the genotyping step, aligning SVs to reference sequences, resulting to File2.csv mentioned in the Figure1C.

**3. [`3_application_on_CRISPResso2Data`](https://github.com/peiyaohu/Hu_et_al_GESV_2025/tree/main/3_application_on_CRISPResso2Data)**

- **Purpose**: To validate the versatility of the GESV pipeline by applying it to data from other CRISPR systems. This analysis uses public datasets processed by the widely-used [CRISPResso2](https://crispresso2.pinellolab.org/) tool.
- **Contents**: Scripts tailored for analyzing different editing outcomes, such as NHEJ, HDR, prime editing and base editing.

**4. [`Figure1D-I`](https://github.com/peiyaohu/Hu_et_al_GESV_2025/tree/main/Figure1D-I)**

- **Purpose**: Contains all the scripts required to regenerate the main figures (Figure 1D-I) and supplementary figures (e.g., Figure S2) presented in the manuscript.
- **Contents**: R  scripts for data visualization and statistical analysis.

### List of pipeline/programs used in the study

| Software and Algorithms | Version   | Purpose                                                      |
| ----------------------- | --------- | ------------------------------------------------------------ |
| seqkit                  | v2.8.0    | Manipulation of FASTA/FASTQ files                            |
| DADA2, R package        | v1.22.0   | Comparative evaluation of SVs calling                        |
| USEARCH                 | v11.0.667 | Comparative evaluation of SVs calling                        |
| VSEARCH                 | v2.30.0   | SVs calling applied in GESV pipeline  (Necessary for GESV pipeline) |
| cutadapt                | v1.9.1    | Primer removal  (Necessary for GESV pipeline)                |
| Biostrings, R package   | v2.62.0   | SVs alignment with reference sequence in GESV pipeline  (Necessary for GESV pipeline) |

## GESV Online Tool: NGS Data Analysis Procedure

### Overview

GESV (Genotyping Edited Sequence Variants) is an integrated online tool for analyzing next-generation sequencing (NGS) data from gene editing experiments. This document provides a step-by-step procedure for using the GESV online platform.

### Input Requirements

### 1. Sample Metadata File

A tab-delimited file containing the following required columns:

| Column                             | Description                                           | Note      |
| ---------------------------------- | ----------------------------------------------------- | --------- |
| **SampleID**                       | Unique identifier for each sequencing file            | Mandatory |
| **MSU locus**                      | Target site identifier (e.g., gene locus)             | Mandatory |
| **WT (5'-3', without primer seq)** | Reference wild-type sequence without primer sequences | Mandatory |
| **Guide_PAM_5prime**               | Guide RNA and PAM sequence (5' end)                   | Mandatory |
| **R1_primer_seq**                  | Forward PCR primer sequence (5'→3')                   | Mandatory |
| **R2_primer_seq**                  | Reverse PCR primer sequence (5'→3')                   | Mandatory |
| **Details**                        | Optional sample information and notes                 | Optimal   |

**Example Metadata:**

| SampleID              | MSU locus      | WT (5'-3')                                                   | Guide_PAM_5prime        | R1_primer_seq             | R2_primer_seq          | Details                |
| --------------------- | -------------- | ------------------------------------------------------------ | ----------------------- | ------------------------- | ---------------------- | ---------------------- |
| SY16011-USR-32598-002 | LOC_Os08g10320 | GTGTTGTCTGAACCATCCTTCGTTCTTTTGTGTAAAAAAAATAGGCATGTGGGGATTCAATGCACTCTCAGGACCAATTCCGA | GTGGGATTCAATGCACTCTCAGG | GTTCTTGCAGTAATAATCAGTTGAT | AAGTTGAGGTTCGTGAGATTCC | RM009-10;R103;LRR-VIII |

### 2. Raw Sequencing Files

- **Format**: FASTQ files (quality encoding: Sanger/Illumina 1.8+)
- **Read Types**:
  - **Paired-end reads**: GESV will automatically merge forward and reverse reads
  - **Single-end reads**: Merging step will be skipped
- **Compression**: GZIP compression supported (.fastq.gz)

### Output Files

After processing, users will receive two main output files:

### 1. `File1.fa` - Sequence Variants per Sample

- **Format**: FASTA format
- **Content**: All detected sequence variants (SVs) with abundance information
- **Naming convention**: `SV{ID};size={abundance_count}`
- **Usage**: Suitable for local multiple sequence alignment and further analysis

**Example:**

```fasta
>SV1;size=1500
GTGTTGTCTGAACCATCCTTCGTTCTTTTGTGTAAAAAAAATAGGCATGTGGGGATTCAATGCACTCTCAGGACCAATTCCGA
>SV2;size=850
GTGTTGTCTGAACCATCCTTCGTTCTTTTGTGTAAAAAAAATAGGCATGTGGGGATTCAATGCACTCTCAGGACCAATTCCGAT
```

### 2. `File2.csv` - Pairwise Alignment Results

- **Format**: CSV (Comma-separated values)
- **Content**: BLASTN-style alignment results comparing SVs against wild-type sequence
- **Columns include**:
  - Sample ID
  - SV identifier
  - Edit type classification(Substitutions (S), Insertions (I), Deletions (D))
  - Mutation details (1D, 2I, 3S)

### Step-by-Step Workflow

1. **Prepare Input Files**
   - Format sample metadata according to the template
2. **Upload to GESV Platform**
   - Access GESV through CRISPR-P website
   - Upload metadata file and sequencing files
   - Select appropriate analysis parameters
3. **Automated Processing**
   - Quality control and adapter trimming
   - Read merging (for paired-end data)
   - Sequence variant calling and genotyping
   - Alignment and edit type classification
4. **Download Results**
   - Retrieve `File1.fa` for variant sequences
   - Download `File2.csv` for detailed alignment results
   - Use outputs for downstream analysis and visualization
