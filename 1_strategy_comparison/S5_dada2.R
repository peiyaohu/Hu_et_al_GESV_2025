### Date: 17 Nov, 2025
# Script for Strategy5

library("dada2")
setwd("/Users/hupeiyao/Documents/Research/LabMember/YSQ/Manuscript/20250818-v2/for_publish/1_Cas9_program/20250821-dada2")
path<-"data/"

list.files(path);

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
saveRDS(fnFs, "dada2_obj/fnFs.RDS")
saveRDS(fnRs, "dada2_obj/fnRs.RDS")

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

saveRDS(filtFs, "dada2_obj/filtFs.RDS")
saveRDS(filtRs, "dada2_obj/filtRs.RDS")

# trim.left, the adaptor sequence is 18nt.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 18, 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)
saveRDS(out, "dada2_obj/out.RDS")


#learn the errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF, "dada2_obj/errF.RDS")
saveRDS(errR, "dada2_obj/errR.RDS")

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, "dada2_obj/dadaFs.RDS")
saveRDS(dadaRs, "dada2_obj/dadaRs.RDS")

dadaFs[[1]]
dadaRs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
saveRDS(mergers, "dada2_obj/mergers.RDS")

head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

saveRDS(seqtab, "dada2_obj/seqtab.RDS")


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track) <- sample.names
head(track)

saveRDS(track, "dada2_obj/track.RDS")
 
################## 20250525 ######################## 
library(Biostrings)
library(dada2)
setwd("/Users/hupeiyao/Documents/Research/LabMember/YSQ/Manuscript/20250818-v2/for_publish/1_Cas9_program/20250821-dada2/")
list.files("dada2_obj/")

dadaFs<- readRDS("dada2_obj/dadaFs.RDS")
dadaRs<- readRDS("dada2_obj/dadaRs.RDS")
mergers <- readRDS("dada2_obj/mergers.RDS")
seqtab_raw <- readRDS("dada2_obj/seqtab.RDS") # 316 5748

track <- readRDS("dada2_obj/track.RDS")
geneinfo<- readRDS("~/Documents/Research/LabMember/YSQ/mosuo/geneinfo.RDS")

seqtab.nochim <- removeBimeraDenovo(seqtab_raw, method="consensus", multithread=TRUE, verbose=TRUE) # 316 5009
saveRDS(seqtab.nochim, "dada2_obj/seqtab.nochim.RDS")

# A DNAStringSet object containing 5009 sequences, with the first three shown (SV_1, SV_2, SV_3), extracted from seqtab.nochim
seq_fa<- DNAStringSet(colnames(seqtab.nochim)) # 316 5009
names(seq_fa)<- paste("SV",1:length(seq_fa),sep="_")
table(colnames(seqtab.nochim) == seq_fa)
colnames(seqtab.nochim)<- names(seq_fa)


library(tidyverse)
datadf<- seqtab.nochim %>% 
  data.frame() %>% rownames_to_column("SampleID") %>%
  pivot_longer(cols = -SampleID, values_to = "count", names_to = "SV") %>%
  filter(count >0)



# 按SampleID分组处理
datadf %>%
  group_by(SampleID) %>%
  group_map(~ {
    # 获取当前样本的SV序列
    sample_seqs <- seq_fa[.x$SV]
    
    # 设置新的序列名称格式
    names(sample_seqs) <- paste0(.x$SV, ";size=", .x$count)
    
    # 写入文件
    writeXStringSet(sample_seqs, filepath = paste0(.y$SampleID, "_SVs.fasta"))
    
  })
