#!/bin/bash

# =============================================================================
# GESV Amplicon Processing Pipeline (Batch mode)
# Description: Workflow for processing amplicon sequencing data from gene-edited samples
# =============================================================================

# -----------------------------------------------
# Parameter Configuration (Set to your own path)
# -----------------------------------------------

# User provided
PRIMER_TABLE="./source/Metadata212.csv" # metadata containing primer and sample information
SEQ_DIR="/public/home/pyhu/project/ysq/20250528/data/" # raw sequencing data dir
VSEARCH=~/biosoft/vsearch-2.30.0/bin/vsearch # VSEARCH software path 
module load cutadapt/1.9.1 # Load cutadapt module

# Output file/dir
ref_seq="all.WT.fa" # output reference sequence file for SV-reference alignment using R/Biostrings later on.
OUTPUT_DIR="./cutadapt_trimmed_files/" # output dir for Cutadapt results (Step1)
mkdir -p "$OUTPUT_DIR" 
File1_DIR="./GESV_File1/" # output dir for GESV_File1 (Step6)
mkdir -p "$File1_DIR"
File_tmp_DIR="TMP_FILE/" # output dir for Step2-5
mkdir -p "$File_tmp_DIR"

THREADS=5

# ----------------------------
# Data Preprocessing
# ----------------------------

# Extract wild-type reference sequences from primer table
cut -f1,3 -d, "$PRIMER_TABLE" | awk -F, 'NR>1 {print ">"$1"\n"$2}' > "$ref_seq"

# ----------------------------
# Main Processing Loop
# ----------------------------

# Process each sample directly from primer table
tail -n +2 "$PRIMER_TABLE" | while IFS=',' read -r SampleID _ _ R1_primer R2_primer _; do
    echo "Processing sample: $SampleID"
    
    # Raw FASTQ file paths
    R1_FASTQ="${SEQ_DIR}${SampleID}_R1.fastq"
    R2_FASTQ="${SEQ_DIR}${SampleID}_R2.fastq"
    
    # Check if input files exist
    if [[ ! -f "$R1_FASTQ" || ! -f "$R2_FASTQ" ]]; then
        echo "Error: Missing FASTQ file $R1_FASTQ or $R2_FASTQ, skipping sample $SampleID"
        continue
    fi

    ## Step 1: Cutadapt Adapter Trimming
    echo "Cutadapt trimming for $SampleID..."
    R1_OUT="${OUTPUT_DIR}${SampleID}_R1_trimmed.fastq"
    R2_OUT="${OUTPUT_DIR}${SampleID}_R2_trimmed.fastq"
    
    cutadapt -g "$R1_primer" -G "$R2_primer" -o "$R1_OUT" -p "$R2_OUT" "$R1_FASTQ" "$R2_FASTQ"

    ## Step 2: Merge Paired-end Reads
    echo "Merging paired-end reads for $SampleID"
    
    $VSEARCH --threads $THREADS \
             --fastq_mergepairs "$R1_OUT" \
             --reverse "$R2_OUT" \
             --fastq_minovlen 10 \
             --fastq_maxdiffs 15 \
             --fastqout "${File_tmp_DIR}${SampleID}.merged.fastq" \
             --fastq_eeout \
             --fastq_allowmergestagger

    ## Step 3: Quality Filtering
    echo "Quality filtering for $SampleID"
    
    $VSEARCH --threads $THREADS \
             --fastq_filter "${File_tmp_DIR}${SampleID}.merged.fastq" \
             --fastq_maxee 1 \
             --fastq_maxns 0 \
             --fastaout "${File_tmp_DIR}${SampleID}.filtered.fasta" \
             --fasta_width 0

    ## Step 4: Sample-level Dereplication
    echo "Dereplicate at sample level for $SampleID"
    $VSEARCH --fastx_uniques "${File_tmp_DIR}${SampleID}.filtered.fasta" \
             --sizeout \
             --fastaout "${File_tmp_DIR}${SampleID}.unique.fasta" 

    ## Step 5: Filter Low-abundance SVs (abundance threshold = 8)
    $VSEARCH --sortbysize "${File_tmp_DIR}${SampleID}.unique.fasta" \
             --minsize 8 \
             --output "${File_tmp_DIR}${SampleID}.min8.svs.fa"

    ## Step 6: Filter SVs less than 20-bp (potential primer dimer)
    $VSEARCH --fastx_filter "${File_tmp_DIR}${SampleID}.min8.svs.fa" \
             --fastq_minlen 20 \
             --fastaout "${File1_DIR}${SampleID}.min8.svs.length20.fa"

    echo "Completed processing for: $SampleID"
    echo "------------------------------------"
done

mv cutadapt_trimmed_files TMP_FILE
 

