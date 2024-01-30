library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(duckplyr)
library(conflicted)
library(RColorBrewer)
conflicts_prefer(duckplyr::filter)

df_og <- duckplyr_df_from_file(
  "D:/ALL_FINAL/Raw_quant/2023-09-14/*.parquet",
  "read_parquet",
  list(hive_partitioning = TRUE)
)

df <- df_og %>% 
  filter(ImageID == "d0224r1p050200f0001") %>% 
  select('Cell_Barcode', 'Frame','factor_median_OBJ_GFP', 'factor_mean_OBJ_GFP', 'factor_GFP_background_Med', 'factor_GFP_background_Avg', 'averageIntensity_GFP_Frame', 'averageIntesntiy_GFP_Background', 'averageIntensity_GFP_Object', 'x80thPercentile_GFP_RAW', 'x60thPercentile_GFP_RAW', 'x90thPercentile_GFP_RAW', 'x95thPercentile_GFP_RAW', 'x99thPercentile_GFP_RAW','max_GFP_RAW') %>%
  tidyr::pivot_longer(cols = (!'Cell_Barcode')&(!'Frame'), cols_vary ="fastest", names_to = 'Base', values_to = 'value')

df$value <- as.integer(df$value)
df$Base <- as.factor(df$Base)

low_comp <- df %>%
  filter(Base %in% c('averageIntesntiy_GFP_Background', 'factor_GFP_background_Avg', 'factor_GFP_background_Med'))

high_comp <- df %>% 
  filter(!Base %in% c('averageIntensity_GFP_Frame', 'averageIntesnity_GFP_Background', 'factor_GFP_background_Avg', 'factor_GFP_background_Med'))

ggplot(high_comp, aes(x = Base, y = value))+
  geom_violin()


sp <- ggbetweenstats(
  data = high_comp,
  x = Base,
  y = value,
  ggplot.component = list(
    scale_color_manual(values = mypal),
    theme(axis.text.x = element_text(angle = 90))
  )
)

ggstatsplot:::.is_palette_sufficient("ggsci", "npg", 10L)
ggstatsplot:::.is_palette_sufficient("RColorBrewer", "Dark2", 6L)
