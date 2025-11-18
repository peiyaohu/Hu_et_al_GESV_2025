# Nov13
# R script for SV-WT alignment (for LRR-RLK 212)

setwd("./")

# Check and install Biostrings (Bioconductor)
if (!require("Biostrings", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings")
}

# Check and install seqinr (CRAN)
if (!require("seqinr", quietly = TRUE))
  install.packages("seqinr")


# Load R package
library(Biostrings)
library(seqinr)

allele_dir <- "./GESV_File1" # path created in GESV_LRR.sh
output_dir <- "./GESV_File2"
log_file <- file.path(output_dir, "processing_log.txt")

# creat output dir
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
all_results_file <- file.path(output_dir, "all_edit_results.csv")

# Initialize log file
writeLines("Processing Log", log_file)
# Initialize list to store results from all samples
all_results_list <- list()
# Load sample metadata containing WT sequences, primer information, and experimental details
sampledata<- read.csv("./source/Metadata212.csv", header = T)


# Find all allele sequence files with the specific pattern
allele_files <- list.files(allele_dir, pattern = "min8.svs.length20.fa$", full.names = TRUE)

# Check if any allele files were found
if (length(allele_files) == 0) {
  stop("Error: No .shift_reads.fa files found in ", allele_dir)
}

index212<- allele_files[gsub(".min8.svs.length20.fa","", basename(allele_files)) %in% sampledata$SampleID]


# Process each allele file sequentially
# allele_files<- head(index212)
for (allele_file in index212) {
  # Extract sample name by removing file extension
  sample_name <- sub(".min8.svs.length20.fa$", "", basename(allele_file))
  
  # Find matching sample information in metadata
  sample_row <- sampledata[sampledata$SampleID == sample_name, ]
  if (nrow(sample_row) == 0) {
    warning_msg <- paste(Sys.time(), "Warning: No matching sample found in sampledata for", sample_name, "\n")
    cat(warning_msg, file = log_file, append = TRUE)
    warning(warning_msg)
    next
  }
  
  # Extract sample-specific information from metadata
  wt_sequence <- sample_row$WT[1]
  sample_id <- sample_row$SampleID[1]
  pam <- sample_row$Guide_PAM_5prime[1]
  R1_seq<- sample_row$R1_primer_seq[1]
  R2_seq<- sample_row$R2_primer_seq[1]
  details<- sample_row$Details[1]

  
  # Read allele sequences from FASTA file
  allele_seqs <- tryCatch({
    readDNAStringSet(allele_file)
  }, error = function(e) {
    warning_msg <- paste(Sys.time(), "Error: Failed to read", sample_name, ":", e$message, "\n")
    cat(warning_msg, file = log_file, append = TRUE)
    warning(warning_msg)
    return(NULL)
  })
  
  # Check if allele sequences were successfully read
  if (is.null(allele_seqs) || length(allele_seqs) == 0) {
    warning_msg <- paste(Sys.time(), "Warning: Allele file is empty or invalid:", sample_name, "\n")
    cat(warning_msg, file = log_file, append = TRUE)
    warning(warning_msg)
    next
  }
  
  allele_names <- names(allele_seqs)
  allele_sequences <- as.character(allele_seqs)
  
  # Extract abundance information from sequence headers (format: ...;size=XXX;...)
  sizes <- as.numeric(gsub(".*;size=([0-9]+).*", "\\1", allele_names))
  # Initialize results data frame to store editing analysis for this sample
  results <- data.frame(
    SampleID = character(),
    Sequence_ID = character(),
    SV_ID = character(),
    Edit_Type = character(),
    Edit_Summary = character(),
    Abundance = numeric(),
    Guide_PAM_5prime = character(),
    R1_primer_seq = character(),
    R2_primer_seq= character(),
    Details = character(),
    Alignment_Pattern = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each allele sequence in the current file
  for (i in 1:length(allele_sequences)) {
    allele_seq <- allele_sequences[i]
    allele_id <- allele_names[i]
    sv_id <- paste("SV_",i,sep="")
    
    size <- sizes[i]
    
    # Perform global-local alignment between allele and wild-type sequence
    # Global-local: entire allele aligned to best matching part of WT
    # Gap opening penalty: -10, Gap extension penalty: -1
    alignment <- tryCatch({
      pairwiseAlignment(allele_seq, wt_sequence, type = "global-local", gapOpening = -10, gapExtension = -1)
    }, error = function(e) {
      warning_msg <- paste(Sys.time(), "Error: Alignment failed for", allele_id, "in", sample_name, ":", e$message, "\n")
      cat(warning_msg, file = log_file, append = TRUE)
      warning(warning_msg)
      return(NULL)
    })
    
    # Skip to next sequence if alignment failed
    if (is.null(alignment)) {
      next
    }
    
    # Extract aligned sequences from alignment object
    pattern <- as.character(alignedPattern(alignment))
    subject <- as.character(alignedSubject(alignment))
    
    # Initialize variables for edit analysis
    edit_type <- "None"
    edit_position <- NA
    edit_detail <- character()
    edit_summary <- ""
    alignment_pattern <- ""
    
    # Convert aligned sequences to character vectors for comparison
    pattern_chars <- strsplit(pattern, "")[[1]]
    subject_chars <- strsplit(subject, "")[[1]]
    len <- min(nchar(pattern), nchar(subject))
    
    # Generate BLAST-style match line (| for matches, space for mismatches)
    match_line <- character(len)
    for (j in 1:len) {
      if (pattern_chars[j] == subject_chars[j] && pattern_chars[j] != "-") {
        match_line[j] <- "|"
      } else {
        match_line[j] <- " "
      }
    }
    
    # Initialize counters for different edit types
    deletion_count <- 0
    insertion_count <- 0
    substitution_count <- 0
    edit_positions <- integer()
    edit_details <- character()
    
    # Scan through alignment to identify all edit events
    j <- 1
    while (j <= len) {
      if (pattern_chars[j] != subject_chars[j]) {
        if (pattern_chars[j] == "-") {
          # 删除（Deletion）
          deletion_count <- deletion_count + 1
          edit_positions <- c(edit_positions, j)
          edit_details <- c(edit_details, paste("Deleted", subject_chars[j], "at position", j))
        } else if (subject_chars[j] == "-") {
          # 插入（Insertion）
          insertion_count <- insertion_count + 1
          edit_positions <- c(edit_positions, j)
          edit_details <- c(edit_details, paste("Inserted", pattern_chars[j], "at position", j))
        } else {
          # 替换（Substitution）
          substitution_count <- substitution_count + 1
          edit_positions <- c(edit_positions, j)
          edit_details <- c(edit_details, paste(subject_chars[j], "->", pattern_chars[j], "at position", j))
        }
      }
      j <- j + 1
    }
    
    # Determine the primary edit type and generate summary
    if (deletion_count > 0 || insertion_count > 0 || substitution_count > 0) {
      edit_type <- "Mixed" # Default for complex edits
      if (deletion_count > 0 && insertion_count == 0 && substitution_count == 0) {
        edit_type <- "Deletion"
      } else if (insertion_count > 0 && deletion_count == 0 && substitution_count == 0) {
        edit_type <- "Insertion"
      } else if (substitution_count > 0 && deletion_count == 0 && insertion_count == 0) {
        edit_type <- "Substitution"
      }

      # Generate compact edit summary (e.g., "3D2I1S" for 3 deletions, 2 insertions, 1 substitution)
      edit_summary_parts <- character()
      if (deletion_count > 0) {
        edit_summary_parts <- c(edit_summary_parts, paste0(deletion_count, "D"))
      }
      if (insertion_count > 0) {
        edit_summary_parts <- c(edit_summary_parts, paste0(insertion_count, "I"))
      }
      if (substitution_count > 0) {
        edit_summary_parts <- c(edit_summary_parts, paste0(substitution_count, "S"))
      }
      edit_summary <- paste(edit_summary_parts, collapse = "")
    } else {
      # No edits detected - sequence matches WT
      edit_type <- "None"
      edit_summary <- "0N"
      edit_details <- c("Identical to WT")
    }
    
    # Get first edit position (if any edits exist)
    edit_position <- if (length(edit_positions) > 0) edit_positions[1] else NA
    
    # Combine all edit details into single string
    edit_detail <- paste(edit_details, collapse = "; ")
    
    # Create BLAST-style alignment visualization (three lines: allele, match symbols, WT)
    alignment_pattern <- paste(
      pattern, "\n",
      paste(match_line, collapse = ""), "\n",
      subject, sep = ""
    )
    
    # Add results for current allele to the data frame
    results <- rbind(results, data.frame(
      SampleID = sample_id,
      Sequence_ID = allele_id,
      SV_ID = sv_id,
      Edit_Type = edit_type,
      Edit_Summary = edit_summary,
      Abundance = size,
      Guide_PAM_5prime <- pam,
      R1_primer_seq<- R1_seq,
      R2_primer_seq<- R2_seq,
      Details = details,
      Alignment_Pattern = alignment_pattern
    ))
  }
  
  # Store results for current sample in the global list
  all_results_list[[sample_name]] <- results

  # Calculate summary statistics for edit types in this sample
  # edit_summary <- aggregate(Abundance ~ Edit_Type, data = results, sum)
  
  # Save individual sample results to CSV files
  output_file <- file.path(output_dir, paste0(sample_name, "_edit_results.csv"))
  # summary_file <- file.path(output_dir, paste0(sample_name, "_edit_summary.csv"))
  write.csv(results, file = output_file, row.names = FALSE)
  # write.csv(edit_summary, file = summary_file, row.names = FALSE)
  
  # Log processing completion
  log_msg <- paste(Sys.time(), "Processed", sample_name, "- Results saved to", output_file, "\n")
  cat(log_msg, file = log_file, append = TRUE)
  
  # Print progress to console
  cat(log_msg)
  # cat("Edit type summary saved to", summary_file, "\n")
}


# Combine results from all samples into a single data frame
if (length(all_results_list) > 0) {
  all_results <- do.call(rbind, all_results_list)
  write.csv(all_results, file = all_results_file, row.names = FALSE)
  cat("All results combined and saved to", all_results_file, "\n")
  cat(paste(Sys.time(), "All results combined and saved to", all_results_file, "\n"), file = log_file, append = TRUE)
} else {
  warning_msg <- paste(Sys.time(), "Warning: No results to combine, all_results.csv not generated\n")
  cat(warning_msg, file = log_file, append = TRUE)
  warning(warning_msg)
}


cat("Batch analysis completed.\n")
cat(paste(Sys.time(), "Batch analysis completed\n"), file = log_file, append = TRUE)
