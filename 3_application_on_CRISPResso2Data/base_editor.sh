#!/bin/bash
VSEARCH=~/biosoft/vsearch-2.30.0/bin/vsearch

    # Step 1: Quality control
    echo "Quality filtering for $SampleID"
    $VSEARCH --fastq_filter "./rawseq/base_editor.fastq" \
             --fastq_maxee 1 \
             --fastq_maxns 0 \
             --fastaout "base_editor.filtered.fasta" \
             --fasta_width 0
    
    # Step 2: Dereplication on sample level
    echo "Dereplicate at sample level and relabel with sample_n for $SampleID"
    $VSEARCH --fastx_uniques "base_editor.filtered.fasta" \
            -sizeout -fastaout "base_editor.unique.fa"
    
    # Step 3: SV calling
    echo "Clustering SVs with cluster_fast for $SampleID"
    $VSEARCH --cluster_fast "base_editor.unique.fa" \
             --id 1 \
             --minseqlength 20 \
             --centroids "base_editor.svs.fa" \
             --uc "base_editor.svs.uc" \
             --sizein \
             --sizeout \
             --relabel SV_ \
             --fasta_width 0 
    
    # Step 4: Filter low-abundant SV
    $VSEARCH --sortbysize "base_editor.svs.fa" \
             --minsize 8 \
             --output "base_editor.min8.svs.fa"     
        
    echo "Completed processing for: CR_BE0"
    echo "------------------------------------"
    

