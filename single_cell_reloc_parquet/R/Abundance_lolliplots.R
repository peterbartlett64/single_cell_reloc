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
library(ggpubr)
library(gghighlight)
library(readxl)
library(tidyr)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")
npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)

# df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet")
df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged.parquet", as_data_frame = T)
#Todo: read in the column with the localization

# df$Protein = as.factor(df$Protein)

# Get the maximum value of Frames_post_treatment for each protein and aggregate the median abundance
summary_abund = df %>% 
  filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>%
  collect() %>% 
  group_by(Protein) %>%
  mutate(Abundance_med = median(Abundance), medLogAbundance = median(LogAbundance)) %>% 
  # summarise(Abundance_med = median(Abundance)) %>%
  ungroup() %>% 
  select(Protein, Abundance_med, medLogAbundance, Single_origin, Double_origin, Single_destination, Double_destination) %>%
  filter(Abundance_med < 180) %>%
  mutate(zScoreGlobal_Abundance = (Abundance_med - mean(Abundance_med, na.rm  = T))/sd(Abundance_med, na.rm = T)) %>%
  mutate(Protein = fct_reorder(Protein, zScoreGlobal_Abundance, .desc = FALSE)) %>%
  drop_na() %>%
  distinct()

summary_abund = df %>% 
  filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>%
  collect() %>% 
  group_by(Protein) %>%
  mutate(Abundance_med = median(Abundance), medLogAbundance = median(LogAbundance)) %>% 
  # summarise(Abundance_med = median(Abundance)) %>%
  ungroup() %>% 
  select(Protein, Abundance_med, medLogAbundance, Single_origin, Double_origin, Single_destination, Double_destination) %>%
  filter(Abundance_med < 180) %>%
  mutate(zScoreGlobal_Abundance = (medLogAbundance - mean(medLogAbundance, na.rm  = T))/sd(medLogAbundance, na.rm = T)) %>%
  mutate(Protein = fct_reorder(Protein, zScoreGlobal_Abundance, .desc = FALSE)) %>%
  drop_na() %>%
  distinct()


#Generate a global plot of the abundances
abundance_plot <- ggplot(summary_abund, aes(y = zScoreGlobal_Abundance, x = Protein)) +
  geom_point(color = npg_clrs[2], size = 1.5, alpha = 0.8) +
  # geom_segment(aes(y = 0, yend = Abundance_med, x = Protein, xend = Protein), colour = npg_clrs[4], linewidth = 0.5) +
  scale_color_npg()+
  theme_minimal() +
  ylab("Median Median Abundance (last frame)")+
  xlab("Proteins")+
  theme_minimal()+
  theme(axis.text.x = element_blank())#, legend.position = "none")
abundance_plot

abundance_plot_fg <- ggplot(filter(summary_abund, Double_destination %in% c('Cytoplasm', 'Cytoplasm and Cytoplasm foci', 'Nuclear Foci', 'Nucleus')), aes(y = zScoreGlobal_Abundance, x = Protein)) +
  geom_point(color = npg_clrs[2], size = 1.5, alpha = 0.8) +
  # geom_segment(aes(y = 0, yend = Abundance_med, x = Protein, xend = Protein), colour = npg_clrs[4], linewidth = 0.5) +
  scale_color_npg()+
  theme_minimal() +
  ylab("Median Median Abundance (last frame)")+
  xlab("Proteins")+
  theme(axis.text.x = element_blank())+#, legend.position = "none")
  facet_wrap(~Double_destination)
abundance_plot_fg

# Read in the abundance data from Brandon's Paper
Real_Abundance <- read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/Abundance_TableS8_MOD.xlsx")
Real_foldChange <- read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/foldAbundance_TableS9_MOD.xlsx")




#Summarize the abundances by time (+/- 15 minutes of target)
twoHours <- df %>%
  filter(Frames_post_treatment >= 14 & Frames_post_treatment <=18) %>% 
  group_by(Protein) %>%
  summarise(Abundance_med = median(LogAbundance)) %>% 
  left_join(Real_Abundance, join_by(Protein == Standard_Name)) %>%
  select(Protein, Abundance_med, TKA_2h_03)
  
fourHours <- df %>%
  filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>% 
  group_by(Protein) %>% 
  summarise(Abundance_med = median(LogAbundance)) %>% 
  left_join(Real_Abundance, join_by(Protein == Standard_Name)) %>%
  select(Protein,Abundance_med, MAZ_4h_03)

abundance_2hr_corr <- ggscatterstats(filter(twoHours, Abundance_med <180 & TKA_2h_03 < 30000),
               x = "TKA_2h_03",
               y = "Abundance_med", 
               title = "Abundance at 2 hours",
               xlab = "Molecule Count (2h 0.03% MMS)",
               ylab = "Median Abundance (2 hours post treatment)",
               point.args = list(colour= npg_clrs[1]),
               label.var = Protein,
               # label.expression = Protein == 'RAD51'| Protein == 'FLR1' | Protein == 'LSM7',
               xsidehistogram.args = list(fill = npg_clrs[2]),
               ysidehistogram.args = list(fill = npg_clrs[3]))
  # ggthemes::theme_clean()
ggsave(sprintf("%s_2hrLogAbundance_plot.png", date), abundance_2hr_corr, width = 10, height = 10, dpi = 300)
ggsave(sprintf("%s_2hrLogAbundance_plot.pdf", date), abundance_2hr_corr, width = 10, height = 20, dpi = 300)


abundance_4hr_corr <- ggscatterstats(filter(fourHours, Abundance_med <180, MAZ_4h_03 < 30000),
               x = "MAZ_4h_03",
               y = "Abundance_med",
               title = "Abundance at 4 hours",
               xlab = "Molecule Count (4h 0.03% MMS)",
               ylab = "Median Abundance (2 hours post treatment)",
               point.args = list(colour= npg_clrs[1]),
               label.var = Protein,
               # label.expression = Protein == 'RAD51'| Protein == 'FLR1' | Protein == 'LSM7',
               xsidehistogram.args = list(fill = npg_clrs[2]),
               ysidehistogram.args = list(fill = npg_clrs[3]))

ggsave(sprintf("%s_4hrLogAbundance_plot.png", date), abundance_4hr_corr, width = 10, height = 10, dpi = 300)
ggsave(sprintf("%s_4hrLogAbundance_plot.pdf", date), abundance_4hr_corr, width = 10, height = 20, dpi = 300)


twoHours <- df %>%
  filter(Frames_post_treatment >= 14 & Frames_post_treatment <=18) %>% 
  group_by(Protein) %>%
  summarise(Abundance_med = median(Abundance)) %>% 
  left_join(Real_Abundance, join_by(Protein == Standard_Name)) %>%
  select(Protein, Abundance_med, TKA_2h_03)

fourHours <- df %>%
  filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>% 
  group_by(Protein) %>% 
  summarise(Abundance_med = median(Abundance)) %>% 
  left_join(Real_Abundance, join_by(Protein == Standard_Name)) %>%
  select(Protein,Abundance_med, MAZ_4h_03)

abundance_2hr_corr <- ggscatterstats(filter(twoHours, Abundance_med <180 & TKA_2h_03 < 30000),
                                     x = "TKA_2h_03",
                                     y = "Abundance_med", 
                                     title = "Abundance at 2 hours",
                                     xlab = "Molecule Count (2h 0.03% MMS)",
                                     ylab = "Median Abundance (2 hours post treatment)",
                                     point.args = list(colour= npg_clrs[1]),
                                     # label.var = Protein,
                                     # label.expression = Protein == 'RAD51'| Protein == 'FLR1' | Protein == 'LSM7',
                                     xsidehistogram.args = list(fill = npg_clrs[2]),
                                     ysidehistogram.args = list(fill = npg_clrs[3]))
# ggthemes::theme_clean()
ggsave(sprintf("%s_2hr_Abundance_plot.png", date), abundance_2hr_corr, width = 10, height = 10, dpi = 300)
ggsave(sprintf("%s_2hr_Abundance_plot.pdf", date), abundance_2hr_corr, width = 10, height = 20, dpi = 300)


abundance_4hr_corr <- ggscatterstats(filter(fourHours, Abundance_med <180, MAZ_4h_03 < 30000),
                                     x = "MAZ_4h_03",
                                     y = "Abundance_med",
                                     title = "Abundance at 4 hours",
                                     xlab = "Molecule Count (4h 0.03% MMS)",
                                     ylab = "Median Abundance (2 hours post treatment)",
                                     point.args = list(colour= npg_clrs[1]),
                                     # label.var = Protein,
                                     # label.expression = Protein == 'RAD51'| Protein == 'FLR1' | Protein == 'LSM7',
                                     xsidehistogram.args = list(fill = npg_clrs[2]),
                                     ysidehistogram.args = list(fill = npg_clrs[3]))

ggsave(sprintf("%s_4hr_Abundance_plot.png", date), abundance_4hr_corr, width = 10, height = 10, dpi = 300)
ggsave(sprintf("%s_4hr_Abundance_plot.pdf", date), abundance_4hr_corr, width = 10, height = 20, dpi = 300)


abundance_Loc_corr <- df %>% 
  select(Protein, MedianProtCorr, Single_origin, Double_origin, Single_destination, Double_destination) %>% 
  distinct()


local_distinct <-  df %>% 
  select("Protein", "Single_origin", "Double_origin", "Single_destination", "Double_destination") %>%
  distinct()


compartmental_Spearmans <- ggbetweenstats(filter(abundance_Loc_corr, Double_destination %in% c('Cytoplasm', 'Cytoplasm and Cytoplasm foci', 'Nuclear Foci', 'Nucleus')),
               x = "Double_destination",
               y = "MedianProtCorr",
               title = "A comparison of Loc-Abundance Spearmans between target localizations",
               xlab = "Localization",
               ylab = "Median Protein Correlation",
               type = 'np')

ggsave(sprintf("%s_compartmental_Spearmans.png", date), compartmental_Spearmans, width = 10, height = 10, dpi = 300)
ggsave(sprintf("%s_compartmental_Spearmans.pdf", date), compartmental_Spearmans, width = 10, height = 20, dpi = 300)




#Plot the points in space and color by localization
ggplot(filter(fourHours, Abundance_med <180, MAZ_4h_03 < 30000)) +
  geom_point(aes(x = MAZ_4h_03, y = Abundance_med))
  
  



