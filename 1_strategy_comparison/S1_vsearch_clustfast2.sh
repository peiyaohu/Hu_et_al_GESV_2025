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
VSEARCH=~/biosoft/vsearch-2.30.0/bin/vsearch
LOG_FILE="vsearch_processing.log"
  
module load cutadapt/1.9.1

  
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
    $VSEARCH --threads $THREADS \
             --fastq_mergepairs "$R1_OUT" \
             --reverse "$R2_OUT" \
             --fastq_minovlen 10 \
             --fastq_maxdiffs 15 \
             --fastqout "${SampleID}.merged.fastq" \
             --fastq_eeout \
             --fastq_allowmergestagger
      
      
    # Step 3: Quality filtering
    echo "Quality filtering for $SampleID"
    $VSEARCH --threads $THREADS \
             --fastq_filter "${SampleID}.merged.fastq" \
             --fastq_maxee 1 \
             --fastq_maxns 0 \
             --fastaout "${SampleID}.filtered.fasta" \
             --fasta_width 0
      
    echo
      
    # Step 4: Dereplicate at sample level
    echo "Dereplicate at sample level and relabel with sample_n for $Sample"
    $VSEARCH --fastx_uniques "${SampleID}.filtered.fasta" \
            -sizeout -fastaout "${SampleID}.unique.fasta" 
 
    # Step 5: Clustering SVs with cluster_fast
    echo "Clustering SVs with cluster_fast for $Sample"
    $VSEARCH --threads $THREADS \
             --cluster_fast "${SampleID}.unique.fasta" \
             --id 1 \
             --iddef 1 \
             --centroids "${SampleID}.svs.fa" \
             --uc "${SampleID}.svs.uc" \
             --sizein \
             --sizeout \
             --relabel SV_ \
             --fasta_width 0
     
    # Step 5: Filter low-abundance SVs (abundance not lower than 8)
    $VSEARCH --sortbysize "${SampleID}.svs.fa" \
             --minsize 8 \
             --output "${SampleID}.min8.svs.fa"
  
    ## Step 6: Filter SVs less than 20-bp (potential primer dimer)
    $VSEARCH --fastx_filter "${SampleID}.min8.svs.fa" \
             --fastq_minlen 20 \
             --fastaout "${SampleID}.min8.svs.length20.fa"

      
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

