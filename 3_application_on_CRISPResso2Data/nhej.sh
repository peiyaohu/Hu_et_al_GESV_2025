#!/bin/bash
VSEARCH=~/biosoft/vsearch-2.30.0/bin/vsearch

    # Step0: Merge PE reads
    $VSEARCH --fastq_mergepairs ../rawseq/nhej.r1.fastq\
             --reverse ./rawseq/nhej.r2.fastq \
             --fastq_minovlen 10 \
             --fastq_maxdiffs 15 \
             --fastqout "nhej.merged.fastq" \
             --fastq_eeout \
             --fastq_allowmergestagger 
    
    # Step1: Quality control
    echo "Quality filtering for $SampleID"
    $VSEARCH --fastq_filter "nhej.merged.fastq" \
             --fastq_maxee 1 \
             --fastq_maxns 0 \
             --fastaout "nhej.filtered.fasta" \
             --fasta_width 0
    
    # Step 2: Dereplication on sample-level
    echo "Dereplicate at sample level and relabel with sample_n for $SampleID"
    $VSEARCH --fastx_uniques "nhej.filtered.fasta" \
            -sizeout -fastaout "nhej.unique.fa"
    
    # Step 3: SV calling
    echo "Clustering OTUs with cluster_fast for $SampleID"
    $VSEARCH --cluster_fast "nhej.unique.fa" \
             --id 1 \
             --minseqlength 20 \
             --centroids "nhej.svs.fa" \
             --uc "nhej.svs.uc" \
             --sizein \
             --sizeout \
             --relabel SV_ \
             --fasta_width 0 

    # Step 4: Filter low-abundant SV
    $VSEARCH --sortbysize "nhej.svs.fa" \
             --minsize 8 \
             --output "nhej.min8.svs.fa"     
    
    echo "Completed processing for: CR-NHEJ"
