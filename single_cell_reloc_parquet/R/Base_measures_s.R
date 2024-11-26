
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(corrplot)
library(smplot2)
library(corrplot)
library(scales)
library(ggsci)
library(ggstatsplot)
library(ggsci)
library(ggforce)

# Module Code--------------------------------------------
# This script is for comparing the difference between the base measures as inputs for the reloc decision. 
# It is to show the distance apart of the percentiles and the normalization. 



setwd("D:/Second Mind/Academic/Project Stuff/Figures")
# knitr::opts_knit$set(root.dir = "D:/Second Mind/Academic/Project Stuff/Figures")
date <- Sys.Date()
npg_clrs <-  pal_npg("nrc", alpha = 0.7)(15)


df <- read_parquet("D:/ALL_FINAL/Quant_d0220r1p170300_primary.parquet")

#Select only the factors and the 
base_measures <- df %>% 
  select('Cell_Barcode', 'Frame','factor_median_OBJ_GFP', 'factor_mean_OBJ_GFP', 'factor_total_OBJ_GFP', 'factor_GFP_background_Med', 'factor_GFP_background_Avg', 'factor_GFP_background_Tot', 'averageIntensity_GFP_Frame', 'averageIntesntiy_GFP_Background', 'averageIntensity_GFP_Object', 'x80thPercentile_GFP_RAW', 'x60thPercentile_GFP_RAW', 'x90thPercentile_GFP_RAW', 'x95thPercentile_GFP_RAW', 'x99thPercentile_GFP_RAW','max_GFP_RAW') %>%
  tidyr::pivot_longer(cols = (!'Cell_Barcode')&(!'Frame'), cols_vary ="fastest", names_to = 'Base', values_to = 'value')

base_measures$value <- as.integer(base_measures$value)

base_between <- ggplot(base_measures, aes(x = Base, y = value, colors = Base)) +
  # col_bin(palette = npg_clrs)+
  geom_sina()

ggsave(sprintf("base_measures_%s.pdf", date), base_between, width = 30, height = 18)
ggsave(sprintf("base_measures_%s.png", date), base_between, width = 30, height = 18)

base_between <- ggbetweenstats(
  data = base_measures,
  x    = Base,
  y    = value,
  type = "np",
  var.equal = FALSE,
  outlier.tagging = TRUE
)
ggsave(sprintf("base_measures_%s.pdf", date), base_between, width = 30, height = 18)
ggsave(sprintf("base_measures_%s.png", date), base_between, width = 30, height = 18)

