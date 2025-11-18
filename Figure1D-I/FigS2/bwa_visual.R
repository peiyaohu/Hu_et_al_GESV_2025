### step1. bwa_ysq ###
# Date: 17 Nov, 2025
# Corresponding to Figure S2

library(tidyverse)
library(gridExtra)
library(readr)
library(scales)
library(GenomicAlignments)
library(Biostrings)

setwd("~/Documents/Research/LabMember/YSQ/Manuscript/20251103-v4/github/20251111_github/2_GESV_on_LRR-RLK/bwa/")


########## PAM's position relative to WT #######
geneinfo<- read.csv("../../2_GESV_on_LRR-RLK/source/Metadata212.csv", header = T) %>%
  select(SampleID, MSU.locus,WT, Guide_PAM_5prime,  R1_primer_seq, R2_primer_seq)


align_df<- data.frame()
for(i in geneinfo$SampleID){
  gi<- geneinfo %>% filter(SampleID == i)
  alignment<- Biostrings::pairwiseAlignment(gi$WT[1], gi$Guide_PAM_5prime[1], type="local")
  alignment_data <- pattern(alignment)
  # Extract start pos
  start_pos <- start(alignment_data@range)
  end_pos <- end(alignment_data@range)
  width<- width(alignment_data@range)
  df<- data.frame(SampleID= i,
                  PAM_start_pos = start_pos,
                  PAM_end_pos = end_pos, PAM_range=width)
  align_df<- rbind(align_df, df)
}

geneinfo<- geneinfo %>% left_join(align_df)
###########

n_min<- 8

# Run bwa.sh first to generate .sam files
# 3'-most position of the first 5'-terminal match block of the reads
plotlist<- list()
samfile_summary<- list()
samfile_sum_filt<- list()
samples<- paste0(geneinfo$SampleID,".aligned.sam")
for(filename in samples){
  # filename<- "SY16011-USR-32598-001.aligned.sam"
  
  # read sam file
  sam_data <- read_tsv(filename, 
                       comment = "@",
                       col_names = c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", 
                                     "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", 
                                     "QUAL", "TAGS","TYPE","VALUE","V15"),
                       na = "*")
  if (all(is.na(sam_data$CIGAR))) {
    next  
  }
  
  sampleid<- strsplit(filename,"\\.")[[1]][1]
  sampleinfo<- geneinfo %>% filter(SampleID ==sampleid)
  
  sam_data$Abundance <- as.numeric(sapply(strsplit(sam_data$QNAME,"="),"[",3))
  head(sam_data)
  
  sam_data %>% dplyr::count(FLAG) %>%
    mutate(FLAG=as.character(FLAG)) %>%
    ggplot(aes(x=FLAG, y=n,fill=FLAG))+geom_bar(stat = "identity", position = "stack")
  

  
  
  # Ensure pass_POS is factor
  sam_df <- sam_data %>%
    mutate(
      pass_FLAG = ifelse(FLAG == 0, 1, 0)
      # ,
      # pass_POS = as.factor(ifelse(POS == 1, 1, 0))  # transform to factor
    )
  

  
  ######## CIGAR ########
  
  # Parse CIGAR
  cigar_ops <- cigarOpTable(sam_df$CIGAR)
  cigar_width <- cigarWidthAlongReferenceSpace(sam_df$CIGAR)
  
  sam_df$total_length <- rowSums(cigar_ops)
  sam_df$aligned_length <- cigar_width
  sam_df$soft_clip_length <- cigar_ops[,"S"] + cigar_ops[,"H"]
  
  # sam_df_filt<- sam_df  %>% filter(pass_FLAG ==1 & pass_POS ==1) 
  sam_df_filt<- sam_df  %>% filter(pass_FLAG ==1 ) 
  
  sam_df_filt<- sam_df_filt %>% 
    mutate(Components = str_extract_all(CIGAR, "\\d+|[A-Z]"))
  
  sam_df_filt<- sam_df_filt %>%
    mutate(
      #  the number preceding the ‘M’ character in the CIGAR string
      cigar_first_matchend_pos = as.numeric(str_extract(CIGAR, "\\d+(?=M)")),
      cigar_first_matchend_pos = ifelse(is.na(cigar_first_matchend_pos), 1, cigar_first_matchend_pos)
    ) 
  
  sam_df_filt<-  sam_df_filt %>%
    mutate(
      type_com = str_extract_all(TYPE, "\\d+|[A-Z]"),
      # the number following the ‘Z’ character in the TYPE string.
      z_number = sapply(type_com, function(x) {
        z_index <- which(x == "Z")[1]
        if (!is.na(z_index) && z_index < length(x)) {
          # check if the string after 'Z' is number
          if (grepl("^\\d+$", x[z_index + 1])) {
            return(as.numeric(x[z_index + 1]))
          }
        }
        return(NA)
      })
    )
  
  pam_start_pos<-  sampleinfo$PAM_start_pos
  pam_end_pos<- sampleinfo$PAM_end_pos 
  
  sam_df_filt<-  sam_df_filt %>% 
    mutate(first_align_end = pmin(cigar_first_matchend_pos, z_number)) %>%
    mutate(Sampleid = sampleid,
                           pam_start_pos = pam_start_pos,
                           pam_end_pos = pam_end_pos) 
 
  

  filter1<- sam_df_filt %>%  filter(Abundance>= n_min )
  
  sam_df_filt %>% nrow() #356
  filter1 %>% nrow #4
  
  sum(sam_df_filt$Abundance) #[1] 3901
  sum(filter1$Abundance) #[1] 3310
  
  
  # sam_df_filt %>% mutate(first_align_end = pmax(cigar_first_matchend_pos, z_number))  %>% View()
  
  # Plot with dual x-axes and vertical line at target_pos
  Figure2.2<- sam_df_filt %>%
    mutate(Edit = ifelse(Abundance<8,"Artifacts",ifelse(total_length == first_align_end, "Reference seq","Seq with variation"))) %>%
    ggplot(aes(x = first_align_end, y = Abundance, color=Edit)) +
    geom_point(size = 2,alpha=1) +
    geom_vline(xintercept = pam_start_pos, linetype = "dashed", color = "black") +  # Vertical line at target_pos
    geom_vline(xintercept = pam_end_pos, linetype = "dashed", color = "black") +  # Vertical line at target_pos
    geom_hline(yintercept = 8, linetype = "dashed", color = "red") +  # Vertical line at target_pos
    
    scale_y_log10() +  # Log scale for ABUNDANCE (integer reads)
    # scale_x_continuous(
    #   name = "Distance to Target (bp)",  # Primary x-axis
    #   sec.axis = sec_axis(~ . + target_pos, name = "Position from 5' End (bp)")  # Secondary x-axis
    # ) +
    labs(y = "Abundance (Reads)", title=paste0("YSQ-CRISPR/Cas9\n[", sampleid,"]"), x="The 3'-most position of the first 5'-terminal match block") +
    # scale_color_manual(values = c("Editing" = "blue", "Error" = "red", "Uncertain" = "gray")) +
    theme_bw() +
    theme(
      axis.title.x.top = element_text(color = "darkgray"),  # Secondary axis label color
      axis.text.x.top = element_text(color = "darkgray")    # Secondary axis tick color
    )+ scale_color_manual(
      values = c(
        "Artifacts" = "grey",
        "Reference seq" = "blue",      # 改为蓝色
        "Seq with variation" = "black" # 改为黑色
      )
    )
  
  plotlist[[filename]]<- Figure2.2
  
  
  
  if(ncol(sam_df_filt)==28){samfile_summary[[filename]]<- sam_df_filt}
  

  if(ncol(filter1)==28){samfile_sum_filt[[filename]]<- filter1}
  
  
}

# combine all results
combined_filt <- do.call(rbind, samfile_sum_filt)
combined_data <- do.call(rbind, samfile_summary)


combined_data %>%
  group_by(Sampleid) %>%
  slice_max(Abundance, n = 1, with_ties = FALSE) %>%
  ungroup() %>% nrow #212

maxN_min = 500
combined_filtMaxN<- combined_data %>%
  group_by(Sampleid) %>%
  slice_max(Abundance, n = 1, with_ties = FALSE) %>%
  ungroup() %>% filter(Abundance > maxN_min)

nrow(combined_filtMaxN) #212




df<- dplyr::count(combined_data, Sampleid, name = "N_raw") %>% 
  arrange(-N_raw) %>% 
  left_join(combined_filt %>% dplyr::count( Sampleid, name = "N_filt")) %>% 
  mutate(reduce_pct = N_filt/N_raw) 


# sample with maxN >500, which is more reliable
filt_ind<- which(gsub(".aligned.sam","",names(plotlist)) %in% combined_filtMaxN$Sampleid)

# Extract the figures
indices<- sample(filt_ind,16)
selected_plots <- plotlist[indices]

# Arrange the figures
plot20<- grid.arrange(grobs = selected_plots, ncol = 4, nrow = 4)
ggsave("Figure1.bwa_selectplots.pdf", plot20,width=20, height = 16)


combined_data %>% data.frame() %>% dplyr::count(Sampleid) %>% 
  arrange(-n) %>%
  mutate(Sampleid = factor(Sampleid, levels=Sampleid)) %>% 
  ggplot(aes(x=Sampleid, y=n))+geom_bar(stat="identity")
  
window=10
combined_data %>%
  filter(first_align_end >(pam_start_pos-window)|first_align_end<(pam_end_pos+window)) %>%
  data.frame() %>% dplyr::count(Sampleid) %>% 
  arrange(-n) %>% head

combined_data %>% data.frame() %>% dplyr::count(Sampleid) %>% 
  arrange(-n) %>%head

save.image("data_20251118.RData")

