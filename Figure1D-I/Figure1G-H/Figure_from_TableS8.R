###
# Date: Nov 13
# Purpose: TableS8 to generate MainFigure of LRR-RLK gene editing plants (FigureS3 and Figure1G-H )
# setwd("/Users/hupeiyao/Documents/Research/LabMember/YSQ/Manuscript/20251103-v4/github/20251111_github/2_GESV_on_LRR-RLK/")

# SupplementaryTable_v4.xlsx could be downloaded from the paper
allele_data<- readxl::read_excel("../../../SupplementaryTable_v4.xlsx", sheet = 6, skip = 2)
allele_data %>% filter(TYPE2=="Allele") %>% dplyr::count(Edit_Type)

library(gridExtra)
library(tidyverse)

##### Figure S3 #####
# We analyed the lines with abundance > 50:
figure1_data<- allele_data %>%
  dplyr::count(SampleID, TYPE)

allele_data_db<- allele_data %>% 
  group_by(SampleID, TYPE) %>%
  summarise(TYPE_abundance = sum(Abundance)) %>% 
  group_by(SampleID) %>%
  mutate(percentage = TYPE_abundance / sum(TYPE_abundance) * 100) 

row_order<- allele_data_db %>% filter(TYPE=="Allele") %>% arrange(-percentage) %>% pull(SampleID)
Figure2<- allele_data_db %>%
  mutate(
    TYPE = factor(TYPE, levels=c("Allele","Noise")),
         SampleID = factor(SampleID, levels=row_order)) %>%
  ggplot(aes(x = SampleID, y = percentage, fill = TYPE)) +geom_bar(stat="identity")+
  theme_light()+labs(y="Percentage of SV's counts", x="SampleID")+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=c("#5dade2","#f39c12"))+
  theme(axis.text = element_text(color="black"), axis.title = element_text(color="black"))


figure1_data<- figure1_data %>% mutate(
  TYPE = factor(TYPE, levels=c("Noise","Allele")),
  SampleID = factor(SampleID, levels=row_order))
Figure1<- figure1_data %>%
  ggplot(aes(x=SampleID, y=n,fill=TYPE))+geom_bar(stat="identity")+
  geom_hline(yintercept = 2)+geom_hline(yintercept = 4)+
  theme_light()+labs(y="Number of SVs", x="SampleID")+
  theme(axis.text.x = element_blank())+
  scale_fill_manual(values=c("#f39c12","#5dade2"))+
  theme(axis.text = element_text(color="black"), axis.title = element_text(color="black"))

p1<- grid.arrange(Figure1, Figure2, ncol=2)

ggsave("FigureS1.TYPE.pdf",p1, width = 10, height = 4)


##### Figure 1G-H #####

allele_data %>% filter(TYPE2=="Allele") %>% dplyr::count(Edit_Type)

allele_data<- allele_data %>% filter(TYPE2=="Allele")

allele_data %>% filter(grepl("TrueAllele", allele_data$Type)) %>% dplyr::count(Edit_Summary) %>%
  extract(Edit_Summary, into = c("Number", "Letter"), 
          regex = "(\\d+)([A-Z]*)", remove = FALSE) %>% View()

N_Allele_distribution<- allele_data %>% 
  dplyr::count(SampleID,name = "Allele") %>% dplyr::count(Allele, name = "N") 

N_Allele_distribution %>%
  mutate(
    Allele = as.factor(Allele),
    percentage = N / sum(N) * 100,
    label = paste0(Allele, "\n", round(percentage, 1), "%")
  ) %>%
  ggplot(aes(x = 2, y = N, fill = Allele)) +
  geom_col(color = "white") +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5),
            size = 3) +
  coord_polar("y", start = 0) +
  xlim(0.5, 2.5) +
  theme_void() +
  labs(title = "Allele Distribution") +
  theme(legend.position = "none")+scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a"))
write.csv(N_Allele_distribution,"FigureS2.N_Allele_distribution.csv")
ggsave("Figure1G.pdf", width=3,height = 3)

a<- allele_data %>% dplyr::count(Edit_Summary) %>%
  extract(Edit_Summary, into = c("Number", "Letter"), 
          regex = "(\\d+)([A-Z]*)", remove = FALSE) %>% 
  filter(Letter %in% c("I","D")) %>%dplyr::select(Number, Letter, n) %>% 
  group_by(Number, Letter) %>% summarise(N2= sum(n)) %>% 
  arrange(Letter) %>% mutate(Number= as.numeric(Number))  
a %>% mutate(N2=ifelse(Letter =="D",-1*N2,N2)) %>% 
  ggplot(aes(x=Number, y=N2, fill=Letter))+geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#d3a900","#2B55A9"))+theme_bw()+
  ylim(c(-150,150))+
  labs(x="Length", y="Number of Alleles")
write.csv(a, "Figure1H.Indel.csv")
ggsave("Figure1H.InDel_distribution.pdf", width=5,height = 3)
