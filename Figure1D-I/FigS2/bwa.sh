for sample in $(cut -d',' -f1 ../../2_GESV_on_LRR-RLK/source/Metadata212.csv | tail -n +2)
do
    echo "Processing sample: $sample"
    
    # Extract ref seq
    seqkit grep -n -p "$sample" all.WT.fa > ref.fasta
    
    # Check if the ref seq is extracted successfully or not
    if [ ! -s ref.fasta ]; then
        echo "Warning: No reference sequence found for $sample, skipping..."
        continue
    fi
    
    # Indexing
    bwa index ref.fasta
    
    # Alignment
    # TMP_FILE/"${sample}.unique.fasta": output of ../../2_GESV_on_LRR-RLK/GESV_LRR.sh 
    bwa mem ref.fasta ../TMP_FILE/"${sample}.unique.fasta" > "${sample}.aligned.sam"
    
    echo "Completed processing for $sample"
    echo "----------------------------------------"
done

