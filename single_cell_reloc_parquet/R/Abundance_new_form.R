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
library(ggExtra)
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
library(robustbase)
library(emmeans)
library(tidymodels)
library(arrow)
library(robcor)
library(performance)
library(lme4)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")
npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)

# df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet")
# df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)

# df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>%
#   filter(Protein == 'SAE2') %>%
#   collect() %>%
#   group_by(Protein, Frames_post_treatment) %>%
#   nest()
#
# df2 <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = T) %>%
#   filter(Protein == 'SAE2') %>%
#   filter(Frames_post_treatment == 20)


# lin_corr <- function(df){
#   # return(cor.test(df$z_score_Abund, df$z_score_Abund, method = 'kendall'))
#   lmer(z_score_Loc ~z_score_logAbund + (Unique_pos), data = df)
# }



# tal_m <- lmer(z_score_Loc ~z_score_logAbund + (Unique_pos | Frames_post_treatment), data = df2)
# m <-lm(data = df2, z_score_Loc ~z_score_logAbund)
# check_model(m)
#
# x <- check_model(tal_m, panel = F)
# ggsave("Check_model.pdf",x, width = 10, height = 10, dpi = 300)
#
# lmer(z_score_Loc ~z_score_logAbund + (Unique_pos | Frames_post_treatment), data = df$data[[47]])
#
# df_nested <- df %>%
#   mutate(model = map(data, lin_corr))
#
# View(df_nested$model[[1]])
# performance(df_nested$model[[1]])


#Todo: read in the column with the localization

# df$Protein = as.factor(df$Protein)
# post_t <- df %>%
#   filter(MMS_localization_class, HU_localization_class, MMS_HU_merged_class, Cell_Barcode, Date,
#          Unique_Frame, Protein, Is_treated, Frames_post_treatment, Unique_pos,
#          Loc_score, Relocalized, Abundance, log_Abundance, log_Loc_score, z_score_Loc,
#          z_score_Abund, z_score_logLoc, z_score_logAbund, Progen_bud, Yet, Does, When, Yes_yet,
#          No_yet, pres_end, LogAbundance, CurrNot, CurrYes, currProportion, RelocVelocity,
#          RelocAcceleration, CurrNotYet, CurrYet, currYetProportion, YetVelocity, YetAcceleration,
#          countDoesNot, countDoes, DoesProportion, DoesVelocity, Does_FinDiff) %>%
#   collect()

#generate individual plots for each protein

# cols <- c('MMS_localization_class', 'HU_localization_class', 'MMS_HU_merged_class', 'Cell_Barcode', 'Date',
#           'Unique_Frame', 'Protein', 'Is_treated', 'Frames_post_treatment', 'Unique_pos',
#           'Loc_score', 'Relocalized', 'Progen_bud', 'Yet', 'Does', 'When')


# prot_t <- df %>%
#   filter(Frames_post_treatment >= 0) %>%
#   select(MMS_localization_class, HU_localization_class, MMS_HU_merged_class, Cell_Barcode, Date,
#           Unique_Frame, Protein, Is_treated, Frames_post_treatment, Unique_pos,
#           Loc_score, Relocalized, Abundance, log_Abundance, log_Loc_score, z_score_Loc,
#           z_score_Abund, z_score_logLoc, z_score_logAbund, Progen_bud, Yet, Does, When, Yes_yet,
#           No_yet, pres_end, LogAbundance, CurrNot, CurrYes, currProportion, RelocVelocity,
#           RelocAcceleration, CurrNotYet, CurrYet, currYetProportion, YetVelocity, YetAcceleration,
#           countDoesNot, countDoes, DoesProportion, DoesVelocity, Does_FinDiff) %>%
#   collect() %>%
#   mutate_at(.vars = cols, as_factor) %>%
#   mutate(model = map(splits, ~lm(LocScore ~ logAbundnace, data = .)),
#          augmented = map(model, augment))



  # cor.test(~ + Frames_post_treatment,
  #           data= .,
  #           method = "spearman",
  #           continuity = FALSE,
  #           conf.level = 0.95)
  # model_Abund_vs <- lm(formula = log_Abundance ~ Loc_score, data = post_t)
  # model_Abund_vs %>% plot()
  #
  # model_Abund_vs %>% emmeans(pairwise ~ Loc_score) %>% plot()
  #
  #
  # glm(data = post_t, formula)
  # model_conf <- check_model(model_Abund_vs)
# }
#
# aov(Loc_score ~ log_Abundance + Frames_post_treatment, data = df %>%
#       filter(Protein == 'SAE2') %>%
#       collect())%>%
#   summary()



# Get the maximum value of Frames_post_treatment for each protein and aggregate the median abundance
# summary_abund = df %>%
#   filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>%
#   collect() %>%
#   group_by(Protein) %>%
#   mutate(Abundance_med = median(Abundance), medLogAbundance = median(LogAbundance)) %>%
#   # summarise(Abundance_med = median(Abundance)) %>%
#   ungroup() %>%
#   select(Protein, Abundance_med, medLogAbundance, Single_origin, Double_origin, Single_destination, Double_destination) %>%
#   filter(Abundance_med < 180) %>%
#   mutate(zScoreGlobal_Abundance = (Abundance_med - mean(Abundance_med, na.rm  = T))/sd(Abundance_med, na.rm = T)) %>%
#   mutate(Protein = fct_reorder(Protein, zScoreGlobal_Abundance, .desc = FALSE)) %>%
#   drop_na() %>%
#   distinct()
#
# summary_abund = df %>%
#   filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>%
#   collect() %>%
#   group_by(Protein) %>%
#   mutate(Abundance_med = median(Abundance), medLogAbundance = median(LogAbundance)) %>%
#   # summarise(Abundance_med = median(Abundance)) %>%
#   ungroup() %>%
#   select(Protein, Abundance_med, medLogAbundance, Single_origin, Double_origin, Single_destination, Double_destination) %>%
#   filter(Abundance_med < 180) %>%
#   mutate(zScoreGlobal_Abundance = (medLogAbundance - mean(medLogAbundance, na.rm  = T))/sd(medLogAbundance, na.rm = T)) %>%
#   mutate(Protein = fct_reorder(Protein, zScoreGlobal_Abundance, .desc = FALSE)) %>%
#   drop_na() %>%
#   distinct()


#Generate a global plot of the abundances
# abundance_plot <- ggplot(summary_abund, aes(y = zScoreGlobal_Abundance, x = Protein)) +
#   geom_point(color = npg_clrs[2], size = 1.5, alpha = 0.8) +
#   # geom_segment(aes(y = 0, yend = Abundance_med, x = Protein, xend = Protein), colour = npg_clrs[4], linewidth = 0.5) +
#   scale_color_npg()+
#   theme_minimal() +
#   ylab("Median Median Abundance (last frame)")+
#   xlab("Proteins")+
#   theme_minimal()+
#   theme(axis.text.x = element_blank())#, legend.position = "none")
# abundance_plot
#
# abundance_plot_fg <- ggplot(filter(summary_abund, Double_destination %in% c('Cytoplasm', 'Cytoplasm and Cytoplasm foci', 'Nuclear Foci', 'Nucleus')), aes(y = zScoreGlobal_Abundance, x = Protein)) +
#   geom_point(color = npg_clrs[2], size = 1.5, alpha = 0.8) +
#   # geom_segment(aes(y = 0, yend = Abundance_med, x = Protein, xend = Protein), colour = npg_clrs[4], linewidth = 0.5) +
#   scale_color_npg()+
#   theme_minimal() +
#   ylab("Median Median Abundance (last frame)")+
#   xlab("Proteins")+
#   theme(axis.text.x = element_blank())+#, legend.position = "none")
#   facet_wrap(~Double_destination)
# abundance_plot_fg

# Read in the abundance data from Brandon's Paper
Real_Abundance <- read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/TableS8_Final_MOD.xlsx")
Real_foldChange <- read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/foldAbundance_TableS9_MOD.xlsx")

df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>% 
  select("Standard Name", Protein, MMS_localization_class, HU_localization_class, MMS_HU_merged_class, Cell_Barcode, Date,
                   Unique_Frame, Protein, Is_treated, Frames_post_treatment, Unique_pos,
                   Loc_score, Relocalized, Abundance, LogAbundance, log_Abundance, log_Loc_score, z_score_Loc,
                   z_score_Abund, z_score_logLoc, z_score_logAbund, currProportion, currYetProportion) %>%
  filter(!(Protein %in% c("PPH22d0210", "LCD1", "RAD5", "RFA1d0210", "LSM3d0210r1", "RRP17d0210", "SLX8", "XRS2d0210", "Contam."))) %>%
  collect() %>%
  rename(Standard_Name = "Standard Name")
  # group_by(Protein, Frames_post_treatment)
  # summarise(Abundance_med = median(Abundance), .by_group = T)


#Summarize the abundances by time (+/- 15 minutes of target)
zeroHours <- df %>%
  filter(Frames_post_treatment <= 0, Frames_post_treatment > -4) %>%
  group_by(Protein) %>%
  summarise(Abundance_med = median(LogAbundance)) %>% 
  # summarise(Abundance_med = median(log_Abundance)) %>%
  left_join(Real_Abundance, by = c("Protein" = "Standard_Name")) %>%
  select(Protein, Abundance_med, UNT_MEAN, UNT_PercentileR) %>% 
  mutate(UNT_PercentileRn = percent_rank(UNT_MEAN))


twoHours <- df %>%
  filter(Frames_post_treatment >= 14 & Frames_post_treatment <=18) %>%
  group_by(Protein) %>%
  summarise(Abundance_med = median(LogAbundance)) %>% 
  # summarise(Abundance_med = median(log_Abundance)) %>%
  left_join(Real_Abundance, by = c("Protein" = "Standard_Name")) %>%
  select(Protein, Abundance_med, TKA_2h) %>% 
  mutate(TKA_PercentileRn = percent_rank(TKA_2h))

fourHours <- df %>%
  filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>%
  group_by(Protein) %>%
  summarise(Abundance_med = median(LogAbundance)) %>% 
  # summarise(Abundance_med = median(log_Abundance)) %>%
  left_join(Real_Abundance, by = c("Protein" = "Standard_Name")) %>%
  select(Protein,Abundance_med, MAZ_4h) %>% 
  mutate(MAZ_PercentileRn = percent_rank(MAZ_4h))



abundance_0hr_corr <- ggscatterstats(filter(zeroHours, UNT_MEAN < 100000),
                                     x = "UNT_MEAN",
                                     y = "Abundance_med",
                                     title = "Abundance at 0 hours",
                                     xlab = "Molecule Count (2h 0.03% MMS)",
                                     ylab = "Median Abundance (0 hours post treatment)",
                                     point.args = list(colour= npg_clrs[1]),
                                     type = "robust",
                                     label.var = Protein,
                                     label.expression = Protein == 'RAD51'| Protein == 'RMI1' | Protein == 'SAE2',  #UNT_MEAN > 100000,
                                     xsidehistogram.args = list(fill = npg_clrs[2]),
                                     ysidehistogram.args = list(fill = npg_clrs[3]))
abundance_0hr_corr
# ggthemes::theme_clean()
ggsave(sprintf("%s_0hrLogAbundance_plot.png", date), abundance_0hr_corr, width = 20, height = 10)
ggsave(sprintf("%s_0hrLogAbundance_plot.pdf", date), abundance_0hr_corr, width = 20, height = 10)
ggsave(sprintf("%s_0hrLogAbundance_plot.eps", date), abundance_0hr_corr, width = 20, height = 10)

abundance_2hr_corr <- ggscatterstats(filter(twoHours, TKA_2h < 100000),
               x = "TKA_2h",
               y = "Abundance_med",
               title = "Abundance at 2 hours",
               type = 'robust',
               xlab = "Molecule Count (2h 0.03% MMS)",
               ylab = "Median Abundance (2 hours post treatment)",
               point.args = list(colour= npg_clrs[1]),
               label.var = Protein,
               label.expression = Protein == 'RAD51'| Protein == 'RMI1' | Protein == 'SAE2',# TKA_2h > 100000,
               xsidehistogram.args = list(fill = npg_clrs[2]),
               ysidehistogram.args = list(fill = npg_clrs[3]))
abundance_2hr_corr
  # ggthemes::theme_clean()
ggsave(sprintf("%s_2hrLogAbundance_plot.png", date), abundance_2hr_corr, width = 20, height = 10)
ggsave(sprintf("%s_2hrLogAbundance_plot.pdf", date), abundance_2hr_corr, width = 20, height = 10)
ggsave(sprintf("%s_2hrLogAbundance_plot.eps", date), abundance_2hr_corr, width = 20, height = 10)


abundance_4hr_corr <- ggscatterstats(filter(fourHours, MAZ_4h < 100000),
               x = "MAZ_4h",
               y = "Abundance_med",
               title = "Abundance at 4 hours",
               xlab = "Molecule Count (4h 0.03% MMS)",
               type = 'robust',
               ylab = "Median Abundance (2 hours post treatment)",
               point.args = list(colour= npg_clrs[1]),
               label.var = Protein,
               label.expression = Protein == 'RAD51'| Protein == 'RMI1' | Protein == 'SAE2',# MAZ_4h > 100000,
               xsidehistogram.args = list(fill = npg_clrs[2]),
               ysidehistogram.args = list(fill = npg_clrs[3]))
abundance_4hr_corr

ggsave(sprintf("%s_4hrLogAbundance_plot.png", date), abundance_4hr_corr, width = 20, height = 10)
ggsave(sprintf("%s_4hrLogAbundance_plot.pdf", date), abundance_4hr_corr, width = 20, height = 10)
ggsave(sprintf("%s_4hrLogAbundance_plot.eps", date), abundance_4hr_corr, width = 20, height = 10)


df_corr_pen <- df %>% 
  group_by(Frames_post_treatment) %>% 
  mutate(Abundance_med = median(LogAbundance)) %>% 
  ungroup() %>% 
  mutate(Frames_post_treatment = as.integer(Frames_post_treatment)) %>% 
  mutate(currProportion = currProportion * 100)

#* Attempt to see if there is a correlation between abundance and penetrance
abundance_curr_penetrance<- grouped_ggscatterstats(filter(df_corr_pen, Frames_post_treatment %in% c(0, 8, 16, 24, 32)),
                                     x = "currProportion",
                                     y = "Abundance_med",
                                     grouping.var = "Frames_post_treatment",
                                     xlab = "Current Proportion",
                                     type = 'robust',
                                     ylab = "Median Abundance (2 hours post treatment)"
                                     # point.args = list(colour= npg_clrs[1]),
                                     # label.expression = Protein == 'RAD51'| Protein == 'RMI1' | Protein == 'SAE2',# MAZ_4h > 100000,
                                     # label.var = Protein)
)
                                     # xsidehistogram.args = list(fill = npg_clrs[2]),
   
abundance_curr_penetrance

df_g <- df %>% 
  group_by(Protein, Frames_post_treatment) %>% 
  mutate(log_Abundance = log(Abundance)) %>% 
  mutate(Abundance_med = median(log_Abundance)) %>% 
  ungroup() %>% 
  mutate(currProportion = currProportion *100) %>% 
  mutate(currYetProportion = currYetProportion * 100)

                                  # ysidehistogram.args = list(fill = npg_clrs[3]))
abundance_curr_penetrance <- ggplot(data = filter(df_g, Frames_post_treatment %in% c(0,8,16,24,32)),
                               mapping = aes(x = currProportion, y = Abundance_med))+
  geom_point() +
  # geom_density2d()+
  sm_statCorr(corr_method = 'pearson', show_text = TRUE, color = 'black')+
  facet_row(~Frames_post_treatment)

abundance_yet_penetrance <- ggplot(data = filter(df_g, Frames_post_treatment %in% c(0,8,16,24,32)),
                               mapping = aes(x = currYetProportion, y = Abundance_med))+
  geom_point() +
  # geom_density2d()+
  sm_statCorr(corr_method = 'pearson', show_text = TRUE, color = 'black')+
  facet_row(~Frames_post_treatment)

ggsave(sprintf("%s_Abundance_curr_plot.png", date), abundance_curr_penetrance, width = 30, height = 10)
ggsave(sprintf("%s_Abundance_curr_plot.pdf", date), abundance_curr_penetrance, width = 30, height = 10)
ggsave(sprintf("%s_Abundance_curr_plot.eps", date), abundance_curr_penetrance, width = 30, height = 10)

ggsave(sprintf("%s_Abundance_yet_plot.png", date), abundance_yet_penetrance, width = 30, height = 10)
ggsave(sprintf("%s_Abundance_yet_plot.pdf", date), abundance_yet_penetrance, width = 30, height = 10)
ggsave(sprintf("%s_Abundance_yet_plot.eps", date), abundance_yet_penetrance, width = 30, height = 10)

abundance_curr_penetrance
abundance_yet_penetrance






#
#

# twoHours <- df %>%
#   filter(Frames_post_treatment >= 14 & Frames_post_treatment <=18) %>%
#   group_by(Protein) %>%
#   summarise(Abundance_med = median(Abundance)) %>%
#   left_join(Real_Abundance, join_by(Protein == Standard_Name)) %>%
#   select(Protein, Abundance_med, TKA_2h_03)
#
# fourHours <- df %>%
#   filter(Frames_post_treatment >= 30 & Frames_post_treatment <=34) %>%
#   group_by(Protein) %>%
#   summarise(Abundance_med = median(Abundance)) %>%
#   left_join(Real_Abundance, join_by(Protein == Standard_Name)) %>%
#   select(Protein,Abundance_med, MAZ_4h_03)
#
# abundance_2hr_corr <- ggscatterstats(filter(twoHours, Abundance_med <180 & TKA_2h_03 < 30000),
#                                      x = "TKA_2h_03",
#                                      y = "Abundance_med",
#                                      title = "Abundance at 2 hours",
#                                      xlab = "Molecule Count (2h 0.03% MMS)",
#                                      ylab = "Median Abundance (2 hours post treatment)",
#                                      point.args = list(colour= npg_clrs[1]),
#                                      # label.var = Protein,
#                                      # label.expression = Protein == 'RAD51'| Protein == 'FLR1' | Protein == 'LSM7',
#                                      xsidehistogram.args = list(fill = npg_clrs[2]),
#                                      ysidehistogram.args = list(fill = npg_clrs[3]))
# # ggthemes::theme_clean()
# ggsave(sprintf("%s_2hr_Abundance_plot.png", date), abundance_2hr_corr, width = 10, height = 10, dpi = 300)
# ggsave(sprintf("%s_2hr_Abundance_plot.pdf", date), abundance_2hr_corr, width = 10, height = 20, dpi = 300)
#
#
# abundance_4hr_corr <- ggscatterstats(filter(fourHours, Abundance_med <180, MAZ_4h_03 < 30000),
#                                      x = "MAZ_4h_03",
#                                      y = "Abundance_med",
#                                      title = "Abundance at 4 hours",
#                                      xlab = "Molecule Count (4h 0.03% MMS)",
#                                      ylab = "Median Abundance (2 hours post treatment)",
#                                      point.args = list(colour= npg_clrs[1]),
#                                      # label.var = Protein,
#                                      # label.expression = Protein == 'RAD51'| Protein == 'FLR1' | Protein == 'LSM7',
#                                      xsidehistogram.args = list(fill = npg_clrs[2]),
#                                      ysidehistogram.args = list(fill = npg_clrs[3]))
#
# ggsave(sprintf("%s_4hr_Abundance_plot.png", date), abundance_4hr_corr, width = 10, height = 10, dpi = 300)
# ggsave(sprintf("%s_4hr_Abundance_plot.pdf", date), abundance_4hr_corr, width = 10, height = 20, dpi = 300)
#
#
# abundance_Loc_corr <- df %>%
#   select(Protein, MedianProtCorr, Single_origin, Double_origin, Single_destination, Double_destination) %>%
#   distinct()
#
#
# local_distinct <-  df %>%
#   select("Protein", "Single_origin", "Double_origin", "Single_destination", "Double_destination") %>%
#   distinct()
#
#
# compartmental_Spearmans <- ggbetweenstats(filter(abundance_Loc_corr, Double_destination %in% c('Cytoplasm', 'Cytoplasm and Cytoplasm foci', 'Nuclear Foci', 'Nucleus')),
#                x = "Double_destination",
#                y = "MedianProtCorr",
#                title = "A comparison of Loc-Abundance Spearmans between target localizations",
#                xlab = "Localization",
#                ylab = "Median Protein Correlation",
#                type = 'np')
#
# ggsave(sprintf("%s_compartmental_Spearmans.png", date), compartmental_Spearmans, width = 10, height = 10, dpi = 300)
# ggsave(sprintf("%s_compartmental_Spearmans.pdf", date), compartmental_Spearmans, width = 10, height = 20, dpi = 300)
#
#
#
#
# #Plot the points in space and color by localization
# ggplot(filter(fourHours, Abundance_med <180, MAZ_4h_03 < 30000)) +
#   geom_point(aes(x = MAZ_4h_03, y = Abundance_med))





