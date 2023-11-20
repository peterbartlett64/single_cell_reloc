library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)



#This was deprecated for similar code in python using the plotly.
#If I decide to go back to using R for this visualization, then it is mostly written.
#The only change that needs to be made is using the the protein pen specific columns

theme_set(
  theme_minimal() +
    theme(legend.position = "right")
)


#df <- read_parquet("D:/ALL_FINAL/Final_combined_comparison.parquet")
#df <-filter(df, Protein == "FLR1")

#This takes in the two dataframes and merges. 
#This would require the suffix to be removed from all the protein names
df_cross <- read_parquet("D:/ALL_FINAL/95th_percentile/All_pos_pt_percentages_melt.parquet")
df_cross %>% 
  mutate(across('Protein', str_replace, '--perc_t', ''))

df_longitude <- read_parquet("D:/ALL_FINAL/95th_percentile/All_pos_pt_percentages_melt.parquet")
df_longitude %>% 
  mutate(across('country', str_replace, 'India', 'Albania'))

df_merged <- merge(df_cross,df_longitude, by="Protein")
inta

df <- read_parquet("D:/ALL_FINAL/95th_percentile/All_pos_pt_percentages_melt.parquet")


df %>% 
  select(Protein, Percentage_reloc, Time_post_treatment) %>%
  filter(Protein == 'FLR1-perc_t', Time_post_treatment<300) %>% 
  ggplot(aes(x = Time_post_treatment, y = Percentage_reloc, group = Protein)) +
    geom_line() +
    scale_color_npg() +
    scale_fill_aaas()
  #geom_line(aes(x = 'Time_post_treatment', y = 'Percentage_reloc'))
