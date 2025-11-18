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


