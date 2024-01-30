# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-01-18
#
# Script Name:    C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc_parquet/R/control_comparison.R
#
# Script Description:
#This is an updatded modificaition from previous file (lollipop_sort_single)
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
library(yaml)
# library(duckplyr)
# library(conflicted)
library(stringr)
library(purrr)
library(RColorBrewer)
library(dichromat)
library(paletteer)
# conflicts_prefer(duckplyr::filter)
library(Polychrome)
#--------------------------------------------

date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

#read_yaml() # This can be corrected later on to take in the global variables script.
subset_by = 'Unique_pos'
subset = list(c(''))
exclude = list(c('d0214r2p160300', 'd0214r2p160200', 'd0214r2p150200', 'd0214r2p150300', 'd0214r2p130300', 'd0214r2p130200', 'd0218r2p1010200'))
control = 'DDC2'

# df_read <- duckplyr_df_from_file(
#   "D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet",
#   "read_parquet"
# )

df_read <- read_parquet("D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet")

df <- df_read %>% 
  dplyr::filter(!subset_by %in% exclude) %>% 
  dplyr::filter(str_like(Protein, "DDC2d0218r2%")) %>%  # This will subset by the protein prefix, because there will be multiple days/suffixes
  dplyr::filter(Frames_post_treatment == 20)

Q <- quantile(df$Loc_score, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(df$Loc_score)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
eliminated<- subset(df, df$Loc_score > (Q[1] - 1.5*iqr) & df$Loc_score < (Q[2]+1.5*iqr))
# max_same_post <-  aggregate(df, by = Protein, FUN = 'max')

Q <- quantile(df$cell_area, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(df$cell_area)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
eliminated<- subset(df, df$cell_area > (Q[1] - 1.5*iqr) & df$cell_area < (Q[2]+1.5*iqr))
dropped<- subset(df, df$cell_area < (Q[1] - 1.5*iqr) | df$cell_area > (Q[2]+1.5*iqr))

colorCount <-  length(unique(df$Protein))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
control_between<- ggbetweenstats(
                      data = dropped,
                      x = Protein,
                      y = cell_area,
                      type = "nonparametric",
                      p.adjust.method = 'bonferroni',
                      pairwise.display = "none",
                      messages = FALSE,
                      package = "Polychrome",
                      palette = "glasbey",
                      # palette = colorRampPalette(brewer.pal(27, "Accent"))(colorCount),
                      var.equal = FALSE,
                      outlier.tagging = TRUE
                      # results.subtitle
)
ggsave(sprintf("Control_comparison_%s.png", date), control_between, width = 30, height = 18)


control_same_day <- ggbetweenstats(
                        data = test,
                        x = Protein,
                        y = cell_area,
                        type = "nonparametric",
                        p.adjust.method = 'bonferroni',
                        pairwise.display = "none",
                        messages = FALSE,
                        package = "Polychrome",
                        palette = "glasbey",
                        # palette = colorRampPalette(brewer.pal(27, "Accent"))(colorCount),
                        var.equal = FALSE,
                        outlier.tagging = TRUE
                        # results.subtitle
)

  

