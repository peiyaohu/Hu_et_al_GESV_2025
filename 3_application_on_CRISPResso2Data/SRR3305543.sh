#!/bin/bash
VSEARCH=~/biosoft/vsearch-2.30.0/bin/vsearch
    
    # Step1: Quality control
    echo "Quality filtering for $SampleID"
    $VSEARCH --fastq_filter "./rawseq/SRR3305543.fastq" \
             --fastq_maxee 1 \
             --fastq_maxns 0 \
             --fastaout "SRR3305543.filtered.fasta" \
             --fasta_width 0
    
    # Step 2: Dereplication on sample-level
    echo "Dereplicate at sample level and relabel with sample_n for $SampleID"
    $VSEARCH --fastx_uniques "SRR3305543.filtered.fasta" \
            -sizeout -fastaout "SRR3305543.unique.fa"
    
    # Step 3: SV calling
    echo "Clustering SVs with cluster_fast for $SampleID"
    $VSEARCH --cluster_fast "SRR3305543.unique.fa" \
             --id 1 \
             --minseqlength 20 \
             --centroids "SRR3305543.svs.fa" \
             --uc "SRR3305543.svs.uc" \
             --sizein \
             --sizeout \
             --relabel SV_ \
             --fasta_width 0 
    
    # Step 4: Filter low-abundant SV
    $VSEARCH --sortbysize "SRR3305543.svs.fa" \
             --minsize 8 \
             --output "SRR3305543.min8.svs.fa"     

        
    echo "Completed processing for: SRR3305543 (CR-BE1)"
    

