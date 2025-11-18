## Script for Figure 1I
# Date: 17 Nov, 2025

setwd("~/Documents/Research/LabMember/YSQ/Manuscript/20251103-v4/github/20251111_github/3_application_on_CRISPResso2Data/CRISPResso_web_res/")
# a<- read.table("CRISPResso_top20.txt", sep='\t')
# 
# head(a$SampleID)
# seqs<- gsub("-","",a$Aligned_Sequence)
# library(Biostrings)
# 
# seqs<- DNAStringSet(seqs)
# abd<- a$NReads
# names(seqs)<- paste(a$SampleID, rep(paste0("SV",1:20),8))
# writeXStringSet(seqs, "CRISPResso_web_top20/all_top20.fa")
# samp<- unique(a$SampleID)
# for(i in 1:8){
#   subseq<- seqs[((i-1)*20+1):(i*20)]
#   names(subseq)<- paste0("CR",1:20,";size=",abd[((i-1)*20+1):(i*20)])
#   writeXStringSet(subseq, paste0("../MUSCLE/",samp[i],".fa"))
# }
# 
# 
# samp_order<- c("3_GESV-BE0.fa", "2_GESV-HDR.fa", "1_GESV-NHEJ.fa", "8_GESV-PE.fa",
#                "4_GESV-BE1.fa", "5_GESV-BE2.fa","6_GESV-BE3.fa","7_GESV-BEuntreat.fa")
# samp_order<- c("3_GESV_BE0.fa", "2_GESV_HDR.fa", "1_GESV_NHEJ.fa", "8_GESV_PE.fa",
#                "4_GESV_BE1.fa", "5_GESV_BE2.fa","6_GESV_BE3.fa","7_GESV_BEuntreat.fa")
# list.files("../MUSCLE/", pattern = "sv")
# for(i in 1:8){
#   fa<- readDNAStringSet(paste0("../MUSCLE/", list.files("../MUSCLE/", pattern = "sv")[i]))
#   fa2<- fa[1:25]
#   n_abd<- str_extract(names(fa2), "(?<=size=)\\d+")
#   names(fa2)<- paste0("SV",1:25,";size=",n_abd)
#   writeXStringSet(fa2,paste0("../MUSCLE/", samp_order[i]))
# }

reads<- readxl::read_excel("../MUSCLE/Corresponding_Reads.xlsx")

nrow(reads)
CR<- reads$Reads[seq(1, nrow(reads), by = 2)]
SV<- reads$Reads[seq(2, nrow(reads), by = 2)]
CR_SV_df<- data.frame(CR = CR, SV= SV, SampleID= c(rep("CR-NHEJ", 20),
                                                   rep("CR-HDR", 20),
                                                   rep("CR-BE0", 20),
                                                   rep("CR-BE1", 20),
                                                   rep("CR-BE2", 20),
                                                   rep("CR-BE3", 20),
                                                   rep("CR-BEuntreat", 20),
                                                   rep("CR-PE", 20)))
overall_cor <- cor(CR_SV_df$CR, CR_SV_df$SV, use = "complete.obs")
# 在图的右上角添加相关系数
cor_label <- paste("r =", round(overall_cor, 5))

ggplot(CR_SV_df, aes(x = log10(SV), y = log10(CR), color = SampleID)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, 
              color = "black", linetype = "dashed", alpha = 0.2)+
  geom_text(
    x = Inf, y = Inf, 
    label = cor_label,
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black"
  ) +
  labs(
    title = "Correlation between CR and SV",
    y = "CRISPResso2 Reads (CR)",
    x = "GESV Reads (SV)"
  ) +
  theme_light()+theme(axis.text = element_text(color="black"),
                      axis.title = element_text(color="black"))


#### Version2
ggplot(CR_SV_df, aes(x = log10(SV), y = log10(CR), color = SampleID)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, 
              color = "black", linetype = "dashed", alpha = 0.2) +
  geom_text(
    x = Inf, y = Inf, 
    label = cor_label,
    hjust = 1.1, vjust = 1.5,
    size = 5, color = "black"
  ) +
  scale_x_continuous(
    name = "GESV Reads (SV)",
    labels = function(x) scales::comma(10^x)
  ) +
  scale_y_continuous(
    name = "CRISPResso2 Reads (CR)",
    labels = function(x) scales::comma(10^x)
  ) +
  labs(
    title = "Correlation between CR and SV"
  ) +
  theme_light() +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )
ggsave("../MUSCLE/Figure1I-v1.pdf",width=5, height = 3.6)
  