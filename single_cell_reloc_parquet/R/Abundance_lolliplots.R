# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-01-28
#
# Script Name:    
#
# Script Description:
#
#
# SETUP ------------------------------------
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)
library(dplyr)
library(hrbrthemes)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(scales)
library(stringr)
library(cowplot)
library(forcats)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")
npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)

df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet")
#Todo: read in the column with the localization
df$Protein = df$Protein %>% as.factor()


# Get the maximum value of Frames_post_treatment for each protein and aggregate the median abundance
summary_abund = df %>% 
  group_by(Protein) %>%
  summarise(Frames_post_treatment = max(Frames_post_treatment), Abundance_med = median(Abundance)) %>%
  mutate(Protein = fct_reorder(Protein, Abundance_med))


#Generate a lollipop plot for Protiens and thier abundances
abundance_corr_plot <- ggplot(summary_abund, aes(x = Abundance_med, y = Protein)) +
  geom_point(size = 2.5, alpha = 0.8, color = npg_clrs[2]) +
  geom_segment(aes(x = 0, xend = Abundance_med, y = Protein, yend = Protein), colour = npg_clrs[4], linewidth = 0.5) +
  scale_color_npg()


ggsave(sprintf("%s_Abundance_corr_plot.png", date), abundance_corr_plot, width = 10, height = 10, dpi = 300)
ggsave(sprintf("%s_Abundance_corr_plot.pdf", date), abundance_corr_plot, width = 10, height = 20, dpi = 300)

  


