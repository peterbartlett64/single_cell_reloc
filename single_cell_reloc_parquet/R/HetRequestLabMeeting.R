# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-02-02
#
# Script Name:    
#
# Script Description:
#
#
# SETUP ------------------------------------
library(ggplot2)
library(ggExtra)
library(arrow)library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)
library(dplyr)
library(hrbrthemes)
library(ggthemes)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(scales)
library(stringr)
library(cowplot)
library(ggforce)
library(statsExpressions)
library(purrr)
# library(ggthemr)
library(extrafont)
# library(showtext)
# showtext_auto()

# Module Code--------------------------------------------
# This is the requested plot from Grant. The ggstats plot will not 


date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")
# font_import(path="C:/Windows/Fonts", prompt=FALSE)
# choose_font("Arial")
# use_font <- "Arial"


# Load the data. This file is very large so bringing in with the arrow backend and will do smallerting
data <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)

npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)
# rq_hetviolin <- function(protein){
#   
#   
#   
# }
smaller<- data %>% 
  filter(Frames_post_treatment >= 0 & Protein == 'ECO1') %>%
  filter(Frames_post_treatment %in% c('0', '10', '20', '30', '42')) %>%
  collect()
  # smaller(ImageID, Loc_score, Does, No_yet, Yes_yet, pres_end, CurrNot,
  #        CurrYes, currProportion, RelocVelocity,RelocAcceleration, CurrNotYet,
  #        CurrYet, currYetProportion, YetVelocity, YetAcceleration, countDoesNot,
  #        countDoes, DoesProportion, DoesVelocity, Does_FinDiff)
  
smaller$Relocalized <- as.factor(smaller$Relocalized)
smaller$Yet <- as.factor(smaller$Yet)
smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
# smaller <- filter(smaller, Frames_post_treatment == "20")

# ggthemr('pale', sour)
# ggthemr_reset()
facet_loc_graph <- smaller %>%
  group_by(Frames_post_treatment, Yet) %>%
  {ggplot(smaller, aes(x = Yet, y = Loc_score)) +
    # geom_violin(aes(fill = Yet), scale = 'width')+
    # geom_point(aes(fill = Relocalized), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.2)) +
    geom_boxplot(aes(group = Yet, fill = Yet), notch = T, alpha = 0.4)+
    geom_violin(aes(group = Yet), scale = 'width', color = npg_clrs[4], fill = NA) +
    geom_sina(aes(color = Relocalized, group = Yet), scale = 'width', jitter_y = F) +
    scale_fill_npg()+
    scale_color_npg()+
    # scale_y_log10(limits = c(0.9,2)) +
    # theme_ipsum() +
    # theme_minimal()+
    theme(legend.position = "bottom") +
    labs(title = "ECO1", y = "Loc_score")+
    scale_x_discrete("Frames_post_treatment", labels = c(
      "0" = "Has not",
      "1" = "Has relocalized far"
    ))+
    theme_light(base_family = "Arial")+
      # labels = c(paste0("Not Yet (", .$CurrNotYet, ")"), paste0("Yes (", .$CurrYet), ")"))+
    # subtitle = results_data$expression[[1]] +
    facet_grid(rows = NULL, cols = vars(Frames_post_treatment), scale = "free_x")}

library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)
library(dplyr)
library(hrbrthemes)
library(ggthemes)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(scales)
library(stringr)
library(cowplot)
library(ggforce)
library(statsExpressions)
library(purrr)
# library(ggthemr)
library(extrafont)
# library(showtext)
# showtext_auto()

# Module Code--------------------------------------------
# This is the requested plot from Grant. The ggstats plot will not 


date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")
# font_import(path="C:/Windows/Fonts", prompt=FALSE)
# choose_font("Arial")
# use_font <- "Arial"


# Load the data. This file is very large so bringing in with the arrow backend and will do smallerting
data <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)

npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)
# rq_hetviolin <- function(protein){
#   
#   
#   
# }
smaller<- data %>% 
  filter(Frames_post_treatment >= 0 & Protein == 'ECO1') %>%
  filter(Frames_post_treatment %in% c('0', '10', '20', '30', '42')) %>%
  collect()
  # smaller(ImageID, Loc_score, Does, No_yet, Yes_yet, pres_end, CurrNot,
  #        CurrYes, currProportion, RelocVelocity,RelocAcceleration, CurrNotYet,
  #        CurrYet, currYetProportion, YetVelocity, YetAcceleration, countDoesNot,
  #        countDoes, DoesProportion, DoesVelocity, Does_FinDiff)
  
smaller$Relocalized <- as.factor(smaller$Relocalized)
smaller$Yet <- as.factor(smaller$Yet)
smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
# smaller <- filter(smaller, Frames_post_treatment == "20")

# ggthemr('pale', sour)
# ggthemr_reset()
facet_loc_graph <- smaller %>%
  group_by(Frames_post_treatment, Yet) %>%
  {ggplot(smaller, aes(x = Yet, y = Loc_score)) +
    # geom_violin(aes(fill = Yet), scale = 'width')+
    # geom_point(aes(fill = Relocalized), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.2)) +
    geom_boxplot(aes(group = Yet, fill = Yet), notch = T, alpha = 0.4)+
    geom_violin(aes(group = Yet), scale = 'width', color = npg_clrs[4], fill = NA) +
    geom_sina(aes(color = Relocalized, group = Yet), scale = 'width', jitter_y = F) +
    scale_fill_npg()+
    scale_color_npg()+
    # scale_y_log10(limits = c(0.9,2)) +
    # theme_ipsum() +
    # theme_minimal()+
    theme(legend.position = "bottom") +
    labs(title = "ECO1", y = "Loc_score")+
    scale_x_discrete("Frames_post_treatment", labels = c(
      "0" = "Has not",
      "1" = "Has relocalized far"
    ))+
    theme_light(base_family = "Arial")+
      # labels = c(paste0("Not Yet (", .$CurrNotYet, ")"), paste0("Yes (", .$CurrYet), ")"))+
    # subtitle = results_data$expression[[1]] +
    facet_grid(rows = NULL, cols = vars(Frames_post_treatment), scale = "free_x")}

ggsave("ECO1_loc_score.pdf", plot = facet_loc_graph, width = 4, height = 3, units = 'in', device = cairo_pdf)

facet_loc_graph
  
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 
# npst <- ggbetweenstats(
#   data = smaller,
#   x    = Yet,
#   y    = Loc_score,
#   type = "np",
#   var.equal = FALSE,
#   outlier.tagging = TRUE,
#   point.args = list(alpha = 0),
#   violin.args = list(width = 0, linewidth = 0)
# ) +
#   geom_sina(aes(color = Relocalized, group = Yet), scale = 'width') +
#   scale_fill_npg()+
#   scale_color_npg()+
#   # scale_y_log10(limits = c(0.9,2)) +
#   theme_ipsum() +
#   theme(legend.position = "top") +
#   labs(title = "FLR1", y = "Loc_score")
# npst

  
# ggwithinstats(
#   data = smaller,
#   x    = Frames_post_treatment,
#   y    = Loc_score,
#   type = "robust"
# )  
#   
#   
# 
# data_mods <- data %>% 
#   group_by(Frames_post_treatment) %>%
#   mutate()
#   mutate(Yet_count = n(filter(., Yet == "Yes"))) %>% 
#   mutate(Yet_display = paste0(, " (n = ", Yet_count, ")")) %>% 
#   mutate(yet)
# data_mod
# 
# 
# rq_plot <- ggpt
#   
# 
# 
# facet_grid(~Frames_post_treatment)
# 
# 
