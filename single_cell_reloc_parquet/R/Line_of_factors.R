# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-02-06
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
library(ggforce)
library(statsExpressions)
library(purrr)
library(ggthemr)
library(extrafont)
library(directlabels)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

#Read in the data, in the arrow format. Done for greater speed
data <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)

npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)
# rq_hetviolin <- function(protein){
#   
#   
#   
# }
smaller<- data %>% 
  filter(Frames_post_treatment >= 0 & Protein == 'ECO1') %>%
  # filter(Frames_post_treatment %in% c('0', '10', '20', '30', '42')) %>%
  collect()
# select(ImageID, Loc_score, Does, No_yet, Yes_yet, pres_end, CurrNot,
#        CurrYes, currProportion, RelocVelocity,RelocAcceleration, CurrNotYet,
#        CurrYet, currYetProportion, YetVelocity, YetAcceleration, countDoesNot,
#        countDoes, DoesProportion, DoesVelocity, Does_FinDiff)


graph <- smaller %>% 
  select(Frames_post_treatment, CurrNot, CurrYes, currProportion, RelocVelocity,RelocAcceleration, CurrNotYet,
       CurrYet, currYetProportion, YetVelocity, YetAcceleration, countDoesNot,
       countDoes, DoesProportion, DoesVelocity) %>% 
  distinct() %>% 
  tidyr::pivot_longer(!Frames_post_treatment, names_to = "Variable", values_to = "Value") %>%
  mutate(label = if_else(Value == max(Value), as.character(Variable), NA_character_)) %>%
  {ggplot(., aes(x = Frames_post_treatment, y = Value))+
      geom_line(aes(color = Variable))+
      scale_colour_discrete(guide = 'none') +
      theme(legend.position = "none")+
      # scale_x_continuous(expand=c(0, 1)) +
      # scale_x_continuous(breaks = scales::pretty_breaks(20)) +
      # # geom_dl(aes(label = Variable), method = list(dl.combine("first.points", "last.points")), cex = 0.8)+
      # geom_dl(aes(label = Variable), method = list(dl.trans(x = x + 0.2), "last.points", cex = 0.8)) +
      # geom_dl(aes(label = Variable), method = list(dl.trans(x = x - 0.2), "first.points", cex = 0.8))
      geom_labelline(aes(group = Variable, colour = Variable, label = Variable), hjust = .7, gap = T, straight = T)
      
  }

simp <- data %>% 
  select(Protein, Frames_post_treatment, RelocVelocity, RelocAcceleration) %>% 
  filter(Frames_post_treatment >= 0) %>% 
  collect() %>% 
  group_by(Protein) %>% 
  #Get the time and maximum of relocAceleration and relocVelocity
  summarise(maxRelocA = max(RelocAcceleration, na.rm = T),
            maxRelocV = max(RelocVelocity, na.rm = T), 
            maxRelocVtime = Frames_post_treatment[which.max(RelocVelocity)], 
            maxRelocAtime = Frames_post_treatment[which.max(RelocAcceleration)])
  # tidyr::pivot_longer(!Protein, names_to = "Variable", values_to = "Value")

# 
# ggplot(data = simp, aes(x = Frames_post_treatment, y = RelocAcceleration))+
#   geom_smooth(aes(group = Protein))


# ggplot(data = simp,aes(x = Protein, y = Value, colour = Protein))+
#   geom_histogram()+
#   facet_wrap(~Variable)

ggbarstats(simp, x = Protein, y = maxRelocA)
  
