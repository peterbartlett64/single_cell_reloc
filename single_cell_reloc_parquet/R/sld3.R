# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-05-30
#
# Script Name: Manual per-frame and timecourse pentrance plot
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
library(ggforce)
library(statsExpressions)
library(purrr)
library(ggthemr)
library(extrafont)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)
show_col(npg_clrs)



sld3_df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>% 
  filter(Protein == 'SLD3') %>% 
  select("Frames_post_treatment","CurrNot", "CurrYes", "currProportion", "CurrNotYet", "CurrYet", "currYetProportion", "YetVelocity", "YetAcceleration") %>% 
  collect()

sld_graph_df <- sld3_df %>%
  group_by(Frames_post_treatment) %>% 
  summarise(currProportion = mean(currProportion),
    currYetProportion = mean(currYetProportion)) %>%
  ungroup() %>%
  replace(is.na(.),0) %>% 
  # mutate(currProportion = currProportion *100,
         # currYetProportion = currYetProportion * 100) %>% 
  filter(Frames_post_treatment >= -10, Frames_post_treatment <= 32)
  

sld3_pen_compare <- ggplot(data = sld_graph_df, mapping = aes(x = Frames_post_treatment)) +
  geom_point(aes(y = currProportion), inherit.aes = T, colour = npg_clrs[2], size = 2.5) +
  geom_point(data= sld3_df %>% filter(Frames_post_treatment >= 0, Frames_post_treatment <= 32),
             mapping = aes(x = Frames_post_treatment, y = currYetProportion), color = npg_clrs[4], size = 2.5)+
  scale_x_continuous(breaks = round(seq(min(-10), max(32), by = 2),1),
                     labels = scales::label_number())+
  scale_y_continuous(trans = log2_trans(),
                     limits = c(0.1, 1),
                     labels = scales::label_percent())+
  theme_cowplot()

sld3_pen_compare

ggsave("Sld3_pen_compare.png", sld3_pen_compare, width = 12, height = 9, device = cairo_pdf)
ggsave("Sld3_pen_compare.eps", sld3_pen_compare, width = 12, height = 9, device = cairo_pdf)
ggsave("Sld3_pen_compare.pdf", sld3_pen_compare, width = 12, height = 9, device = cairo_pdf)
  


