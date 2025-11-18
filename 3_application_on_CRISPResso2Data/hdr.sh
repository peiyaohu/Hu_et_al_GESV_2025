#!/bin/bash
VSEARCH=~/biosoft/vsearch-2.30.0/bin/vsearch

    # Step1: Quality control
    echo "Quality filtering for $SampleID"
    $VSEARCH --fastq_filter ../rawseq/hdr.fastq \
             --fastq_maxee 1 \
             --fastq_maxns 0 \
             --fastaout "hdr.filtered.fasta" \
             --fasta_width 0
    
    # Step 2: Dereplication on sample-level
    echo "Dereplicate at sample level and relabel with sample_n for $SampleID"
    $VSEARCH --fastx_uniques "hdr.filtered.fasta" \
            -sizeout -fastaout "hdr.unique.fa"
    
    # Step 3: SV calling
    echo "Clustering OTUs with cluster_fast for $SampleID"
    $VSEARCH --cluster_fast "hdr.unique.fa" \
             --id 1 \
             --minseqlength 20 \
             --centroids "hdr.svs.fa" \
             --uc "hdr.svs.uc" \
             --sizein \
             --sizeout \
             --relabel SV_ \
             --fasta_width 0 

    # Step 4: Filter low-abundant SV
    $VSEARCH --sortbysize "hdr.svs.fa" \
             --minsize 8 \
             --output "hdr.min8.svs.fa"     
    
    echo "Completed processing for: CR-HDR"
        
