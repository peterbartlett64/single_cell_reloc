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
# df = read_parquet("D:/ALL_FINAL/Combined_by_perc/")
df = read_parquet("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc_parquet/Visualizations/FLR1_with_stuff.parquet")
date <- Sys.Date() # Capture a date element so that plots can be saved automatically

#Create a new df that has just the cells at time +0 and +42 (last frame of treatment).
# This will be used to see if the same cells are moving
pstg <- df %>%
  filter(Frames_post_treatment %in% c("42", "0"))

#This reveals that there are non-unique values of Yet
#There are cells which are pre-loc and there are cells which are not present at
#both timepoints
pstg %>%
  group_by(Cell_Barcode, Yet) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  filter(n > 1L)

start_yes <- pstg %>% 
  filter(Frames_post_treatment == '0' & Yet == "Yes") %>% 
  select(Cell_Barcode)

start_no <- pstg %>% 
  filter(Frames_post_treatment == '0' & Yet == "No") %>% 
  select(Cell_Barcode)

#Ending
end_yes <- pstg %>% 
  filter(Frames_post_treatment == '42' & Yet == "Yes") %>% 
  select(Cell_Barcode)

end_no <- pstg %>% 
  filter(Frames_post_treatment == '42' & Yet == "No") %>% 
  select(Cell_Barcode)

start_yes_g <- ggwithinstats(
  data = filter(pstg, Cell_Barcode %in% start_yes$Cell_Barcode),
  x = Frames_post_treatment,
  y = Loc_score,
  type = 'np'
)
start_yes_g
ggsave(sprintf("startYES_g_%s.pdf", date), start_yes_g, width = 30, height = 18)
ggsave(sprintf("startYES_g_%s.png", date), start_yes_g, width = 30, height = 18)


start_no_g <- ggwithinstats(
  data = filter(pstg, Cell_Barcode %in% start_no$Cell_Barcode),
  x = Frames_post_treatment,
  y = Loc_score,
  type = 'np'
)
start_no_g
ggsave(sprintf("startNO_g_%s.pdf", date), start_no_g, width = 30, height = 18)
ggsave(sprintf("startNO_g_%s.png", date), start_no_g, width = 30, height = 18)


#* Compare whether there is a difference in the groups that start above threshold and those that start below
both_no_g <- ggwithinstats(
  data = filter(pstg, Cell_Barcode %in% start_no$Cell_Barcode)%>%
    filter(Cell_Barcode %in% end_no$Cell_Barcode),
  x = Frames_post_treatment,
  y = Loc_score,
  type = 'np')
both_no_g

both_yes_g <- ggwithinstats(
  data = filter(pstg, Cell_Barcode %in% start_yes$Cell_Barcode & Cell_Barcode %in% end_yes$Cell_Barcode),
  x = Frames_post_treatment,
  y = Loc_score,
  type = 'np'
)
both_yes_g

start_yes_end_no <- ggwithinstats(
  data = filter(pstg, Cell_Barcode %in% start_yes$Cell_Barcode & Cell_Barcode %in% end_no$Cell_Barcode),
  x = Frames_post_treatment,
  y = Loc_score,
  type = 'np'
)
start_yes_end_no


counts <- df %>%
  group_by(Frames_post_treatment) %>% 
  count()


ggplot(df, aes(x = Frames_post_treatment, y = Loc_score))+
  geom_boxplot(outlier.colour = 'red',
               outlier.shape = 3,
               outlier.size= 2,
               notch = TRUE,
               varwidth = TRUE)+
  scale_x_discrete(labels = paste(counts$Frames_post_treatment, "\n n = ", counts))




ggsave(sprintf("startNO_g_%s.pdf", date), start_no_g, width = 30, height = 18)
ggsave(sprintf("startNO_g_%s.png", date), start_no_g, width = 30, height = 18)









