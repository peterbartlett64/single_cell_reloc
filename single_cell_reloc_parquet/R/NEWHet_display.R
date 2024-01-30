# require(pacman) #* Load the require package and then the remaining dependancies
#pacman:: p_load(ggdist,dplyr, GGally, ggplot2, ggthemes, ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr, ggpmisc)
library(ggplot2)
library(tidyverse)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(corrplot)
library(smplot2)
library(corrplot)
library(ggsci)
library(car)
library(dlookr)
library(wesanderson)

setwd("D:/Second Mind/Academic/Project Stuff/Figures") #Save to folder in vault that will be copied to Dropbox
df = read_parquet("D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet")
# df = read_parquet("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc_parquet/Visualizations/FLR1_with_stuff.parquet")
date <- Sys.Date() # Capture a date element so that plots can be saved automatically


#* Add a new column that denotes whether the cell ultimately does show reloc.
#* This needs to be done before the number is switched to a str factor
smaller <- df %>%
  filter(Protein == 'RAD51') %>% 
  select('Protein', 'Frames_post_treatment', 'Cell_Barcode', 'Loc_score', "Yet", "Relocalized", 'Frame', 'Yes_yet', 'No_yet') %>% #* Do the subletting here rather than later
  group_by(Cell_Barcode) %>% 
  mutate(min_value = max(Yet)) %>% 
  ungroup()
max_frame = max(smaller$Frames_post_treatment)
smallerProtein = unique(smaller$Protein)

#Need to change the values before changing to factor class
smaller$Yet <- replace(smaller$Yet, smaller$Yet == 0, 'No')
smaller$Yet <- replace(smaller$Yet, smaller$Yet == 1, 'Yes')

#The read in file needs to have some columns switched to categorical
smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
smaller$Relocalized <- as.factor(smaller$Relocalized)
smaller$Frame<- as.factor(smaller$Frame)
# smaller$Yes_yet <- as.factor(smaller$Yes_yet)
# smaller$No_yet <- as.factor(smaller$No_yet)
smaller$Yet <- as.factor(smaller$Yet)

#Get a subset of the data to see if normally distributed. Can we use a parametric test or not


#* Get the final value to compare between groups
single <- smaller %>% 
  filter(Frames_post_treatment == max_frame) #changed this to 39 Manually. Could find the highest for Protein


#Check to see if data is normally distributed. In this case, it is not, so must do a non-parametric test
single %>% 
  group_by(Yet) %>%
  normality(Loc_score)

#The levene test to check if the variances are the same between groups
leveneTest(Loc_score ~ Yet, single)


# Generate a plot to see the non-parametric difference at the final frame between the groups (Does and Doesn't {Yet = 1, Yet = 0}).
# Under the assumption that the distributions are not normally distributed. Besed on the normality call above
npst <- ggbetweenstats(
  data = single,
  x    = Yet,
  y    = Loc_score,
  type = "np",
  var.equal = FALSE,
  outlier.tagging = TRUE
)

ggsave(sprintf("npst_%s_%s.pdf", date, smallerProtein), npst, width = 30, height = 18)
ggsave(sprintf("npst_%s_%s.png", date, smallerProtein), npst, width = 30, height = 18)

# Generate a plot to the the Wilk's t-test comparison between the groups.
abst <- ggbetweenstats(
  data = single,
  x    = Yet,
  y    = Loc_score,
  type = "robust"
)
ggsave(sprintf("abst_%s_%s.pdf", date,smallerProtein), abst, width = 30, height = 18)
ggsave(sprintf("abst_%s_%s.png", date, smallerProtein), abst, width = 30, height = 18)


sp <- grouped_ggbetweenstats(
  data = smaller,
  x = Yet,
  y = Loc_score,
  grouping.var = Frames_post_treatment,
)
ggsave(sprintf("grouped_between_%s_%s.pdf", date, smallerProtein), sp, width = 30, height = 18)
ggsave(sprintf("grouped_between_%s_%s.png", date, smallerProtein), sp, width = 30, height = 18)



# pstg <- pstg %>% 
#   pivot_wider(
#     id_cols = Cell_Barcode,
#     names_from = Yet,
#     values_from = Loc_score) %>% 
#   mutate(difference = Yes - No)

shift <- ggwithinstats(
  data = pstg,
  x = Frames_post_treatment,
  y = Loc_score, 
  type = 'np'
)
ggsave(sprintf("shiftLoc_st_end_%s_%s.pdf", date, smallerProtein), shift, width = 30, height = 18)
ggsave(sprintf("shiftLoc_st_end_%s_%s.png", date, smallerProtein), shift, width = 30, height = 18)






