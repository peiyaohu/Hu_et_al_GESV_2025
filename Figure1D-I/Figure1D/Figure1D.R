# date: 17 Nov
### Figure 1D ####
setwd("/Users/hupeiyao/Documents/Research/LabMember/YSQ/Manuscript/20251103-v4/github/20251111_github/1_strategy_comparison")
library(ComplexHeatmap)
data1d<- readxl::read_excel("Figure1D.xlsx")
data1d %>%
  select(-c("SampleID", "SV_ID")) %>%
  as.matrix() %>%
  Heatmap(cluster_rows = F, cluster_columns = F)



ht_fig1d<- data1d %>%
  select(-c("SampleID", "SV_ID")) %>%
  as.matrix() 


ht_fig1d %>%  pheatmap::pheatmap(cluster_rows = F, cluster_cols = F, display_numbers = T)

library(pheatmap)
library(RColorBrewer)

# Convert data to percentage (assuming ht_fig1d contains raw values)
# Multiply by 100 to represent percentages
ht_fig1d_percent <- ht_fig1d * 100

# Define custom breaks for the color scale
# Assuming your data is now in percentage (0% to 400%)
breaks <- c(0, 0.01, 1,99,100)  # Breaks for 0-1% and 1-400%
colors <- colorRampPalette(brewer.pal(9, "YlGnBu"))(length(breaks))  # Sequential YlGnBu palette


# Plot the heatmap
pheatmap::pheatmap(
  ht_fig1d_percent,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f%%",  # Use custom format for percentages
  breaks = breaks,               # Custom breaks for color scaling
  color = colors               # Custom YlGnBu color palette
  
)
ht_fig1d[,c(2,1,3:7)] %>%
  pheatmap::pheatmap(
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = TRUE,
    number_format = "%.1f%%",
    number_color = "black",  
    color = c("#ffffcc", "#41b6c4", "#0c2c84"),  
    breaks = c(0, 0.01, 1, 4) 
  )

