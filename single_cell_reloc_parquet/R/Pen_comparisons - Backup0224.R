# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-02-01
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
library(ggstatsplot)
library(tidylog)
library(ggrepel)
library(forcats)
library(CGPfunctions)
# library(networkD3)
library(ggalluvial)
library(ggsankey)
library(alluvial)







# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")


#Read in the dataframe and correct the concatonated protein name

name = "All_pos_pt_t_less_percentages_melt"
fix = "-perc_yet"

read_correct <- function(name, fix){
  #Read in the data. There were corrections made so that the all this selecting stuff does not have to done in the future
  penetrance_values_95 <- read_parquet(sprintf("D:/ALL_FINAL/95th_percentile/%s.parquet", name), as_data_frame = T) %>%
    mutate(Protein = stringr::str_replace(Protein, fix, "")) %>%
    mutate(Protein = as.factor(Protein))
  # penetrance_values_95$From  <- as.factor("95")

  penetrance_values_99 <- read_parquet(sprintf("D:/ALL_FINAL/99th_percentile/%s.parquet", name), as_data_frame = T)  %>%
    mutate(Protein = stringr::str_replace(Protein, fix, "")) %>%
    mutate(Protein = as.factor(Protein))
  # penetrance_values_99$From  <- as.factor("99")

  #Read in the table of values to use.
  selector <- read_parquet("D:/ALL_FINAL/percentage_to_use_lib.parquet",  as_data_frame = T) %>%
    mutate(Protein = as.factor(Protein))

  selector_95  <- selector %>% filter(Selected_series == "95")
  selector_99  <- selector %>% filter(Selected_series == "99")

  #Join the tables together
  penetrance_values_95 <- left_join(selector_95, penetrance_values_95, by = "Protein")
  penetrance_values_99 <- left_join(selector_99, penetrance_values_99, by = "Protein")

  penetrance_values <- penetrance_values_95 %>%
    bind_rows(penetrance_values_99)
  return(penetrance_values)
}

penetrance_values <- read_correct("All_pos_pt_t_less_percentages_melt", "-perc_yet") %>%
  group_by(Frames_post_treatment) %>% 
  mutate(hexile = ntile(Percentage_reloc_less, 6)) %>% 
  drop_na()
  

subset_controler <- function(df){
  min_val <- df %>% 
    group_by(Protein) %>%
    summarise(max = max(Frames_post_treatment)) %>%
    filter(max >= 30) %>% #Had to add this temporarily because there were two proteins with very low post_count
    ungroup() %>% 
    summarise(min = min(max))
  subbed <- filter(penetrance_values, Frames_post_treatment %in% c(0, min_val)) # min_val%/%2, min_val))
  return(subbed)
}

new <- subset_controler(penetrance_values) %>%
  arrange(Time_post_treatment) %>%
  mutate(Frames_post_treatment = as.factor(Frames_post_treatment)) %>%
  mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>% 
  mutate(Time_post_treatment = as.factor(Time_post_treatment)) %>% 
  mutate(hexile = as.factor(hexile))
  # mutate(Time_post_treatment = as.factor(Time_post_treatment))
# new$Frames_post_treatment <- as.character(new$Frames_post_treatment)


# newggslopegraph(dataframe = new,
#   Times = Frames_post_treatment,
#   Measurement = Percentage_reloc_less,
#   Grouping = Protein)
# 


# is_alluvia_form(as.data.frame(new), axes = 3:6, silent = TRUE)
# ggplot(as.data.frame(new),
#        aes(y = Percentage_reloc_less, axis1 = Frames_post_treatment, axis2 = hexile)) +
#   geom_alluvium(aes(fill = Admit), width = 1/12) +
#   geom_stratum(width = 1/12, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#   scale_fill_brewer(type = "qual", palette = "Set1") +
#   ggtitle("UC Berkeley admissions and rejections, by sex and department")
sankeyed_data <- subset_controler(penetrance_values)%>%
  mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>% 
  group_by(Protein) %>% 
  arrange(Time_post_treatment, .by_group = T) %>%
  # group_by(Protein) %>% 
  mutate(Frames_post_treatment = as.factor(Frames_post_treatment)) %>%
  mutate(nxt_cat = hexile[row_number()+1]) %>%
  mutate(diff = hexile - nxt_cat) %>%
  mutate(diff = replace_na(diff, 0)) %>%
  mutate(diff_dup = diff) %>% 
  mutate(hexile_cpy = hexile) %>%
  mutate(hexile = as.factor(hexile))%>%
  mutate(hexile = recode_factor(hexile, '1' = 'First', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
  mutate(nxt_cat = as.factor(nxt_cat)) %>%
  mutate(hexile = fct_reorder(hexile, hexile_cpy, .desc = T)) %>%
  mutate(nxt_cat = fct_reorder(nxt_cat, Percentage_reloc_less, .desc = T)) %>%
  mutate(nxt_x = Frames_post_treatment[row_number()+1]) %>% 
  mutate(diff = as.factor(diff)) %>%
  mutate(diff = recode_factor(diff, `0` = "No Change", `1` = "Up 1", `2` = "Up 2", `3` = "Up 3", `4` = "Up 4", `5` = "Up 5", `-1` = "Down 1", `-2` = "Down 2", `-3` = "Down 3", `-4` = "Down 4", `-5` = "Down 5")) %>% 
  mutate(diff = fct_reorder(diff, diff_dup, .desc = T)) %>% 
  ungroup()
  # drop_na()#This is to get rid of the last row for each protein. As written this can be used for serval timepoins instead of just the first and last
sankeyed_data

#try something new
ggplot(data = sankeyed_data, aes(x = Frames_post_treatment,
                                 stratum = hexile,
                                 alluvium = Protein,
                                 fill = hexile,
                                 label = hexile))+
  scale_fill_npg()+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = 'darkgray', width = 1/10)+
  geom_stratum(width = .1)+
  theme(legend.position = "bottom")+
  ggtitle('Testing crap')


ggplot(data = sankeyed_data, aes(x = Frames_post_treatment,
                                 stratum = hexile,
                                 alluvium = Protein,
                                 fill = hexile,
                                 label = hexile))+
  scale_fill_npg()+
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = 'darkgray', width = 1/10)+
  geom_stratum(width = .1)+
  theme(legend.position = "bottom")+
  ggtitle('Testing crap')



# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)

ggplot(data = sankeyed_data %>% make_long(hexile, diff),
       mapping = aes(x = Frames_post_treatment,
                     next_x = nxt_x,
                     node = hexile,
                     next_node = nxt_cat,
                     group = Protein,
                     fill = factor(hexile)))+
  geom_sankey()+
  # geom_sankey_label()+
  theme_sankey()


sankey_sub <- sankeyed_data %>% 
  filter((diff_dup > 0 )| (diff_dup < 0)) #Remove the proteins which did not move
  
  

long_sankey <- sankey_sub %>%   #sankeyed_data %>%
  ungroup() %>% 
  filter(Frames_post_treatment == '0')%>%
  make_long(., hexile, diff)
long_sankey$node <- as.factor(long_sankey$node)
long_sankey$node <- factor(long_sankey$node, 
                           levels=c('Sixth','Fifth','Fourth','Third',
                                    'Second','First', "Down 5", "Down 4",
                                    "Down 3", "Down 2", "Down 1", "Up 1",
                                    "Up 2", "Up 3", "Up 4", "Up 5")) # Got rid of the No Change for subset
                           # levels=c('Sixth','Fifth','Fourth','Third','Second','First', "Down 5", "Down 4", "Down 3", "Down 2", "Down 1", "No Change", "Up 1", "Up 2", "Up 3", "Up 4", "Up 5"))

ggplot(data = long_sankey,
       mapping = aes(x = x,
                     next_x = next_x,
                     node = node,
                     next_node = next_node,
                     fill = node,
                     label = node))+
  geom_sankey()+
  # geom_sankey_label()+
  theme_sankey()+
  geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.5)
  # scale_fill_npg()

df <- mtcars %>%
  make_long(cyl, vs, am, gear, carb)



penetrance_values$Frames_post_treatment <- as.factor(penetrance_values$Frames_post_treatment)


#Function to compare the starting penetrance of each protein. Get an idea of the extent of priming.

pen_comparison <- function(df, time_sub){
  ggplot(data = filter(df, Time_post_treatment == time_sub) %>%
                         distinct() %>%
                         mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
         mapping = aes(x = Protein, y = Percentage_reloc_less))+
    geom_point(mapping = aes(color = Selected_series), alpha = 0.5) + #inherit the aesthetics from the parent ggplot
    theme(legend.position = "none",
          axis.text.x=element_blank())+
    scale_fill_npg()
    # geom_text_repel(aes(label = Protein))#remove the legend
}



pen_comparison(penetrance_values, 0)
pen_comparison(penetrance_values, 1)


#Testing the 

df_testing <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)





#Create the plot of proteins which have a GO term of interest
# ggplot(data = filter(df, Protein %in% Protein %in% c("RFA1", "HTA2", "RDH54", "POL30",
#                                               "NPL4", "CHK1", "MRC1", "MSH3", "MMS21",
#                                               "RAD51", "RAD24", "RTT107", "RRD1", "SLD2",
#                                               "DOA1", "ECO1", "RPN4", "SUB2", "RAD57",
#                                               "DBF4", "PPH3", "RAD55", "CDC1", "RAD9",
#                                               "XRS2", "LRS4", "SLD3", "INO80", "RAD54",
#                                               "SAE2", "ZIP2", "ARP4", "DPB11", "SRS2",
#                                               "CDC6", "RFC2", "SLX4", "TOP3", "NEJ1",
#                                               "RAD33", "RAD52", "NAM7", "YKU80", "PSO2",
#                                               "SGS1", "MRE11", "YKU70", "MGS1", "RAD50",
#                                               "RFC3", "RFC4", "RTS1", "EXO1", "ULS1",
#                                               "REV1", "RAD53", "DDC1", "RIM1")),
#        mapping = aes(x = Frames_post_treatment, y = Loc_score, group = 'Protein'))+
#   geom_line() #inherit the aesthetics from the header