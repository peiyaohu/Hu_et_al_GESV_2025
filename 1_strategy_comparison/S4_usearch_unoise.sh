#BSUB -J vsearch_unique_for_CRISPR/Cas9-edited_samples
#BSUB -n 3
#BSUB -R span[hosts=1]
#BSUB -q q2680v2
#BSUB -o %J.out
#BSUB -e %J.err
#!/bin/bash
  
# Parameters Configuration
PRIMER_TABLE="./demo_data/metadata.csv"
SEQ_DIR="./demo_data/"
OUTPUT_DIR="./cutadapt_trimmed_files/"
mkdir -p "$OUTPUT_DIR"
LOG_DIR="./log"
mkdir -p "$LOG_DIR"
  
THREADS=5
USEARCH=$(which usearch11.0.667_i86linux32)
LOG_FILE="usearch_processing.log"
  
module load cutadapt/1.9.1
module load seqkit
  
echo "====================================" | tee -a "$LOG_FILE"
echo "Vsearch processing started - $(date)" | tee -a "$LOG_FILE"
echo "====================================" | tee -a "$LOG_FILE"
  
# Create temporary file for loop processing
PRIMER_INFO_FILE="primer_info.tmp"
tail -n +2 "$PRIMER_TABLE" | awk -F',' '{print $1","$4","$5}' > "$PRIMER_INFO_FILE"
  
# Process each sample
while IFS= read -r line; do
    IFS=',' read -r SampleID R1_primer R2_primer <<< "$line"
    {
    echo "Processing sample: $SampleID"
      
    DONE_FLAG="$LOG_DIR/${SampleID}.done"
    # Skip if already processed (checkpoint/restart support)
    if [[ -f "$DONE_FLAG" ]]; then
        echo "[$SampleID] Already finished, skip."
        continue
    fi
  
    # Original fastq file paths
    R1_FASTQ="${SEQ_DIR}${SampleID}_R1.fastq"
    R2_FASTQ="${SEQ_DIR}${SampleID}_R2.fastq"
      
    if [[ ! -f "$R1_FASTQ" || ! -f "$R2_FASTQ" ]]; then
        echo "Error: Missing FASTQ file $R1_FASTQ or $R2_FASTQ, skipping sample $SampleID"
        continue
    fi
  
    ## Step 1: Cutadapt adapter trimming
    echo "Cutadapt trimming for $SampleID..."
    R1_OUT="${OUTPUT_DIR}${SampleID}_R1_trimmed.fastq"
    R2_OUT="${OUTPUT_DIR}${SampleID}_R2_trimmed.fastq"
      
    cutadapt -g "$R1_primer" -G "$R2_primer" -o "$R1_OUT" -p "$R2_OUT" "$R1_FASTQ" "$R2_FASTQ"
      
    # Step 2: Merge paired-end reads
    echo "Merging paired-end reads for $SampleID"
    $USEARCH -fastq_mergepairs "$R1_OUT" \
             -reverse "$R2_OUT" \
             -fastqout "$SampleID.merged.fastq"

      
      
    # Step 3: Quality filtering
    echo "Quality filtering for $SampleID"
    $USEARCH -fastq_filter "$SampleID.merged.fastq" \
            -fastq_maxee 1.0 \
            -fastaout "$SampleID.merged.filtered.fasta"
      
    echo
      
    # Step 4: Dereplicate at sample level
    echo "Dereplicate at sample level and relabel with sample_n for $Sample"
    $USEARCH -fastx_uniques "$SampleID.merged.filtered.fasta" \
             -sizeout -fastaout "$SampleID.unique.fa"
      
    # Step 5: Denoising with UNOISE
    echo "Denoising with UNOISE for $Sample"
    $USEARCH -unoise3 "$SampleID.unique.fa" \
             -zotus "$SampleID.zsvs.min8.fa" \
             -minsize 8 \
             -relabel SV -tabbedout $SampleID.unoise3.txt
    
    ################################################################################################################
    # Zsvs' abundance according to USEARCH/UNOISE Man. (https://drive5.com/usearch/manual/pipe_svtab.html)
    # "Input should be reads before quality filtering and before discarding low-abundance unique sequences, e.g. singletons. "
    
    $USEARCH -otutab $SampleID.merged.fastq -zotus $SampleID.zsvs.min8.fa  -otutabout  $SampleID.zsvtab.txt -mapout $SampleID.zmap.txt

    # Step6.1 Calculate Zsv abundance (rowSum)  
     awk 'BEGIN{FS=OFS="\t"} 
         NR>1 {sum=0; for(i=2;i<=NF;i++) sum+=$i; print $1, sum}' "$SampleID.zsvtab.txt" > "$SampleID.zsv_abundance.tmp"

    # Step6.1 Appened the abundances to sequences' identifiers
    awk 'NR==FNR {abundance[$1]=$2; next} 
         /^>/ {zsv_id=substr($0,2); 
               if(zsv_id in abundance) {print ">" zsv_id ";size=" abundance[zsv_id] ";"} 
               else {print ">" zsv_id ";size=0;"}; 
               next} 
          {print}' "$SampleID.zsv_abundance.tmp" "$SampleID.zsvs.min8.fa" > "$SampleID.zsvs.min8.with_abundance.fa"
    ############################################ This is INACCURATE! for CRISPR/Cas9 amplicon #######################



    ############# Customization procedure ######

    ## Extract SVs from $SampleID.unoise3.txt files

    #awk -F'\t' '$2=="denoise" {seq=$1; sub(/;size=[0-9]+;/, "", seq); print seq}' "$SampleID.unoise3.txt" > "$SampleID.denoise_seqs.tmp"
    
    awk -F'\t' '$2=="denoise" {print $1}' "$SampleID.unoise3.txt" > "$SampleID.denoise_seqs.tmp"

    # Extract sequence using seqkit grep command
    seqkit grep -f "$SampleID.denoise_seqs.tmp" "$SampleID.unique.fa" > "$SampleID.denoise_seqs.fa"

    # remove tmp file
    rm -f "$SampleID.denoise_seqs.tmp"
    #########################################
    
    # No command available in USEARCH11 for fasta filter (https://drive5.com/usearch/manual/cmds_all.html) 
    seqkit seq -m 20 "$SampleID.denoise_seqs.fa" > "$SampleID.denoise_seqs_length20.fa" 


#    seqkit seq -m 20 "$SampleID.zsvs.min8.fa" > "$SampleID.zsvs.min8_length20.fa"
#    seqkit seq -m 20 $SampleID.zsvs.min8.with_abundance.fa > $SampleID.zsvs.min8_length20.with_abundance.fa 
    # Create completion flag file
    touch "$DONE_FLAG"
          
    echo "Completed processing for: $SampleID"
    echo "------------------------------------"
    } >> "$LOG_FILE" 2>&1
      
done < "$PRIMER_INFO_FILE"
  
# Clean up temporary files
rm -f "$PRIMER_INFO_FILE"
  
echo "====================================" | tee -a "$LOG_FILE"
echo "Vsearch processing completed - $(date)" | tee -a "$LOG_FILE"
echo "====================================" | tee -a "$LOG_FILE"
