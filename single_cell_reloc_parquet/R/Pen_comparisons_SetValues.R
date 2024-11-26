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
# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(alluvial)
library(gt)
library(data.table)
library(purrr)
library(binr)
library(ggforce)
library(Hmisc)
library(ggbeeswarm)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

#Read in the dataframe and correct the concatonated protein name

name = "All_pos_pt_t_less_percentages_melt"
fix = "-perc_yet"
# npg_clrs <-  pal_npg("nrc", alpha = 0.7)
# show_col(npg_clrs)

#, The below reads in the file to be used based on the dictionary of used values and then get the file.
#. Because there was updated penetrance values, those are being used instead.
# read_correct <- function(name, fix){
#   #Read in the data. There were corrections made so that the all this selecting stuff does not have to done in the future
#   penetrance_values_95 <- read_parquet(sprintf("D:/ALL_FINAL/95th_percentile/%s.parquet", name), as_data_frame = T) %>%
#     mutate(Protein = stringr::str_replace(Protein, fix, "")) %>%
#     mutate(Protein = as_factor(Protein))
#   # penetrance_values_95$From  <- as_factor("95")
# 
#   penetrance_values_99 <- read_parquet(sprintf("D:/ALL_FINAL/99th_percentile/%s.parquet", name), as_data_frame = T)  %>%
#     mutate(Protein = stringr::str_replace(Protein, fix, "")) %>%
#     mutate(Protein = as_factor(Protein))
#   # penetrance_values_99$From  <- as_factor("99")
# 
#   #Read in the table of values to use.
#   selector <- read_parquet("D:/ALL_FINAL/percentage_to_use_lib.parquet",  as_data_frame = T) %>%
#     mutate(Protein = as_factor(Protein))
# 
#   selector_95  <- selector %>% filter(Selected_series == "95")
#   selector_99  <- selector %>% filter(Selected_series == "99")
# 
#   #Join the tables together
#   penetrance_values_95 <- left_join(selector_95, penetrance_values_95, by = "Protein")
#   penetrance_values_99 <- left_join(selector_99, penetrance_values_99, by = "Protein")
# 
#   penetrance_values <- penetrance_values_95 %>%
#     bind_rows(penetrance_values_99)
#   return(penetrance_values)
# }


space_rank = function(value, space, bins){
  # return((ceiling(value/(space/bins))))
  return((bins + 1) - (ceiling(value/(space/bins))))
}

#* This does not seem to be in use. Comment in cleaning
# pentile_rank <- function(df, var){
#   df <- df %>%
#     group_by(Frames_post_treatment) %>%
#     mutate(pentile_rank = ntile(Percentage_reloc_less, 6)) %>%
#     mutate(bin_full_space = space_rank(Percentage_reloc_less, 100, 6)) %>%
#     #* THis had to be turned off below
#     # mutate(bin_full_space = ceiling(Percentage_reloc_less / 10)) %>% #This is to place each protein the possible space
#     drop_na()
#   return(df)
# }


#Rank the penetrance values. This is right now only being done for the less
pen_ranks <- function(df, space, bins){
  df <- df %>%
    group_by(Frames_post_treatment) %>%
    mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>%
    
    #Generate the bin_ranks (Ranked based on library of proteins with equal membership size)
    mutate(cumm_rank = cume_dist(Percentage_reloc_less)) %>%
    
    #Rank based on the total possible size (bins of equal span rather than membership). Because penetrance,0-100 percent 
    mutate(bin_full_space = space_rank(Percentage_reloc_less, 100, bins)) %>%
    
    
    #Rank based on the occupied space of the library (bins of equal span).
    # This is the version that ranks the proteins based on the used space rather than 100. 
    #* Had to move to the cut_interval function rather than the custom space rank when range is between 40-100
    mutate(bin_use_space = cut_interval(Percentage_reloc_less, n = 5, labels = c('5','4','3','2','1'))) %>% 
    # mutate(ten_bin_use_space = cut_interval(Percentage_reloc_less, n = 10, labels = c('10', '9', '8', '7', '6', '5', '4', '3','2', '1'))) %>% 
           
    # mutate(bin_use_space = space_rank(Percentage_reloc_less, space = max(Percentage_reloc_less) - min(Percentage_reloc_less), bins)) %>% 
    # mutate(ten_bin_use_space = space_rank(Percentage_reloc_less, space = max(Percentage_reloc_less) - min(Percentage_reloc_less), bins)) %>% 
    
    #* This is breaks the proteins into bins of equal membership
    mutate(pentile_rank = space_rank(cumm_rank, space = 1, bins = 5)) %>%
    mutate(decile_rank = space_rank(cumm_rank, space = 1, bins = 10)) %>% 
    ungroup()
    # mutate(bin_full_space = as_factor(bin_full_space))
    # mutate(bin_full_space = fct_reorder(bin_full_space, integer()))


    # mutate(bin_full_space = ceiling(Percentage_reloc_less / (space / bins))) #%>% #This is to place each protein the possible space
    # drop_na()
  return(df)
}

penetrance_values_pre <-  read_parquet("D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet", as_data_frame = F) %>% 
  select(Protein, Cell_Barcode, Frames_post_treatment, Loc_score, Relocalized) %>% 
  collect() %>% 
  filter(Frames_post_treatment == '-1') %>% 
  mutate(reloc_cpy = Relocalized) %>% 
  group_by(Protein, Relocalized) %>% 
  summarise(Group = length(reloc_cpy)) %>% 
  mutate(Frames_post_treatment = -1) %>% 
  ungroup() %>% 
  group_by(Protein) %>% 
  mutate(Percentage_reloc_less = Group/sum(Group)) %>% # This isn't really the timecourse penetrance
  mutate(Percentage_reloc_less = Percentage_reloc_less *100) %>% 
  filter(Relocalized == 1) %>% 
  select(Protein, Percentage_reloc_less, Frames_post_treatment) %>% 
  ungroup()


penetrance_values <- read_parquet("D:/ALL_FINAL/Combined_by_perc/new_percs.parquet") %>% #read_correct("All_pos_pt_t_less_percentages_melt", "-perc_yet") %>%
  rename(Percentage_reloc_less = 'Yet_perc')%>%
  select(Protein, Frames_post_treatment, Percentage_reloc_less) %>% 
  bind_rows(penetrance_values_pre) %>% 
  filter(!(Protein %in% c("LCD1", "RAD5", 'DDC2d0224r1', 'Contam.'))) %>%
  filter(!(Protein %like% 'DDC2')) %>%
  filter(!(Protein %like% 'd0210')) %>%
  filter(!(Protein %in% c("SLX4", "MSB3", "SIP5", "FGV2", "ZPR1", "HOS2", "RAD53", "CTR86", "MKT1", "SRP68", "BMH2", "EXO70"))) %>% 
  mutate(Frames_post_treatment = as.integer(Frames_post_treatment)) %>% 
  arrange(Protein, Frames_post_treatment) %>% 
  subset(Frames_post_treatment <= 32) %>% 
  pen_ranks(space = 100, bins = 5) #rank the proteins on the total possible space which is 0-100


penetrance_values
write_parquet(penetrance_values, "mod_penetrancevalues.parquet")

  # mutate(bin_full_space = ((5 + 1) - (ceiling(Percentage_reloc_less/(100/5)))))
  # pen_ranks(., space = 100, bins = 5)
  # group_by(Frames_post_treatment) %>%
  # mutate(pentile_rank = ntile(Percentage_reloc_less, 6)) %>%
  # mutate(bin_full_space = space_rank(Percentage_reloc_less, 100, 6)) %>%
  # mutate(bin_full_space = ceiling(Percentage_reloc_less / 10)) %>% #This is to place each protein the possible space
  # drop_na()



#* Subset controller to get only the two timepoints that we are interested in
subset_controller <- function(df, min_val = NULL, max_val = NULL){
  #* The automated calculation of min(max) was turned off
  # max_val <- ifelse(is.null(min_val), df %>%
  #                     group_by(Protein) %>%
  #                     summarise(max = max(Frames_post_treatment)) %>%
  #                     filter(max >= 30) %>%
  #                     ungroup() %>%
  #                     summarise(min = min(max)), min_val)
  
  subbed <- filter(penetrance_values, Frames_post_treatment %in% c(min_val, max_val))
  return(subbed)
}

#. I don't think that this is being used
# new <- subset_controller(penetrance_values) %>%
#   mutate(Time_post_treatment = Frames_post_treatment *7.5) %>%
#   arrange(Time_post_treatment) %>%
#   mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
#   mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>%
#   mutate(Time_post_treatment = as_factor(Time_post_treatment)) %>%
#   mutate(pentile_rank = as_factor(pentile_rank)) %>%
#   #Put the bin_full_space in ascending order regardless of when they show up in the data
#   mutate(bin_full_space = as_factor(bin_full_space))
  #Put these in ascending order regardless of when they show up in the data

  # mutate(bin_full_space = fct_relevel(bin_full_space, sort(integer())))



# mutate(Time_post_treatment = as_factor(Time_post_treatment))
# new$Frames_post_treatment <- as.character(new$Frames_post_treatment)


# newggslopegraph(dataframe = new,
#   Times = Frames_post_treatment,
#   Measurement = Percentage_reloc_less,
#   Grouping = Protein)

# is_alluvia_form(as.data.frame(new), axes = 3:6, silent = TRUE)
# ggplot(as.data.frame(new),
#        aes(y = Percentage_reloc_less, axis1 = Frames_post_treatment, axis2 = pentile_rank)) +
#   geom_alluvium(aes(fill = Admit), width = 1/12) +
#   geom_stratum(width = 1/12, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#   scale_fill_brewer(type = "qual", palette = "Set1") +
#   ggtitle("UC Berkeley admissions and rejections, by sex and department")

direction_sankey <- function(start, end){
  return(paste0(start, "->", end))
}


#* Separate the three bids into Low, Mid and High bins relative to library
#* This only works if it is using a multiple of three in the reference bin
three_bin <- function(val){
  #Don't forget that this is a rank, so high values are low in rankings
  if(val >= 5){
    return("Low")
  }else if(val <= 4 & val>= 3){
    return("Mid")
  }else if(val <= 2){
    return("High")
  }
  else{
    return(NA)

  }
}


#* Create the function to group proteins based on whether they were above or below in priming/end
#* This was changed by request of Grnat to the make it more simple. The thresholds we decided upon were 20 and 95.
yet_threshold <- function(val, fr_post, target_s = 20, target_e = 95){
  if (fr_post <= 0){
    if (val >= target_s){
      return("Above")
    }
    if (val < target_s){
      return("Below")
    }
  }
  if (fr_post > 0){
    if (val >= target_e){
      return("Above")
    }
    if (val < target_e){
      return("Below")
    }
  }
}
  

#* Create the secondary secondary binning
sankey_it <- function(df){
  sankeyed_data <- df%>%
    mutate(Percentage_reloc_less = round(Percentage_reloc_less)) %>% 
    # mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>%
    group_by(Protein) %>%
    arrange(Frames_post_treatment, .by_group = T) %>%
    
    #* Pure threshold by set values. Likely set to >=15% for priming and < 95 for incomplete penetrance
    #* This was added by request from Grant for simplification of calls for schema grouping
    mutate(time_thresh = map2(Percentage_reloc_less, Frames_post_treatment, .f = yet_threshold)) %>%
    mutate(time_thresh = unlist(time_thresh)) %>% 
    mutate(nxt_thresh = time_thresh[row_number()+1]) %>% 
    mutate(thresh_changes = map2(time_thresh, nxt_thresh, direction_sankey)) %>%
    mutate(thresh_changes = unlist(thresh_changes)) %>% 
    mutate(thresh_changes = as_factor(thresh_changes)) %>%
    mutate(time_thresh = as_factor(time_thresh)) %>% 
    
    
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
    
    #* Seperate bins into High, Medium, Low. This only work if divisible by 3 for equal size. 
    #* Moving to 5 will mean decision must be made about whether 3 is the middle or 2-4 is
    mutate(three_bin = map(pentile_rank, three_bin)) %>%
    mutate(nxt_three = three_bin[row_number()+1]) %>%

    #* Generate the next equal member rank
    mutate(nxt_pent = pentile_rank[row_number()+1]) %>%
    mutate(diff_pent = pentile_rank - nxt_pent) %>%
    mutate(diff_pent = replace_na(diff_pent, 0)) %>%
    mutate(diff_pent_dup = diff_pent) %>%
    mutate(pentile_cpy = pentile_rank) %>%
    mutate(pentile_rank = as_factor(pentile_rank))%>%
    mutate(pentile_rank = recode_factor(pentile_rank, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth')) %>%
    mutate(diff_pent = as_factor(diff_pent)) %>%
    mutate(diff_pent = recode_factor(diff_pent, '0' = "No Change", '1' = "Up 1", '2' = "Up 2", '3' = "Up 3", '4' = "Up 4", '-1' = "Down 1", '-2' = "Down 2", '-3' = "Down 3", '-4' = "Down 4")) %>%
    mutate(nxt_pent = as_factor(nxt_pent)) %>%
    mutate(nxt_pent  = recode_factor(nxt_pent, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
    mutate(Changes_full = map2(pentile_rank, nxt_pent, direction_sankey)) %>%
    mutate(Changes_full = unlist(Changes_full)) %>%
    mutate(Changes_full = as_factor(Changes_full)) %>%
    
    #* Generate the equal space rank
    mutate(nxt_decSp = bin_full_space[row_number()+1]) %>%
    #Get the difference between the prior frame present (there can be spaced for subseting) and the current
    mutate(diff_decSp = bin_full_space - nxt_decSp) %>%
    mutate(diff_decSp = replace_na(diff_decSp, 0)) %>% #This is to deal with the final value where there is nothing to compare against
    mutate(bin_full_space = as_factor(bin_full_space)) %>%
    
    #Generate the next used space rank
    mutate(nxt_useSp = bin_use_space[row_number()+1]) %>% 
    mutate(diff_useSp = as.integer(bin_use_space) - as.integer(nxt_useSp)) %>%
    mutate(diff_useSp = replace_na(diff_useSp, 0)) %>%
    mutate(bin_use_space = recode_factor(bin_use_space, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth')) %>%
    mutate(diff_useSp = as_factor(diff_useSp)) %>%
    mutate(diff_useSp = recode_factor(diff_useSp, '0' = "No Change", '1' = "Up 1", '2' = "Up 2", '3' = "Up 3", '4' = "Up 4", '-1' = "Down 1", '-2' = "Down 2", '-3' = "Down 3", '-4' = "Down 4")) %>%
    
    #* Get the nxt_x values asuming that Frames_post_treatment is on the x
    mutate(nxt_x = Frames_post_treatment[row_number()+1]) %>%
   

    mutate(Changes_three = map2(three_bin, nxt_three, direction_sankey)) %>%
    mutate(Changes_three = unlist(Changes_three)) %>%
    mutate(Changes_three = as_factor(Changes_three)) %>%


    #Ungroup the dataframe
    ungroup()
    # mutate(diff_pent = fct_reorder(diff_pent, diff_pent_dup, .desc = T)) %>%
    # drop_na() #This is to get rid of the last row for each protein. As written this can be used for several time points instead of just the first and last

  return(sankeyed_data)
}


# sankey_var_it <- function(df, var){
#   sankeyed_data <- df%>%
#     mutate(var = as.integer(var)) %>%
#     group_by(Protein) %>%
#     arrange(Frames_post_treatment, .by_group = T) %>%
#     # group_by(Protein) %>%
#     mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
#     mutate(Time_post_treatment = as_factor(Time_post_treatment)) %>%
#     mutate(nxt_pent = pentile_rank[row_number()+1]) %>%
#     mutate(nxt_decSp = bin_full_space[row_number()+1]) %>%
#     mutate(diff_decSp = bin_full_space - nxt_decSp) %>%
#     mutate(changes = map2_chr(direction_sankey, bin_full_space, nxt_decSp)) %>%
#     mutate(diff_decSp = replace_na(diff_decSp, 0)) %>%
# 
#     mutate(diff_pent = pentile_rank - nxt_pent) %>%
#     mutate(diff_pent = replace_na(diff_pent, 0)) %>%
#     # mutate(diff_pent_dup = diff_pent) %>%
# 
#     # mutate(pentile_cpy = pentile_rank) %>%
#     mutate(pentile_rank = as_factor(pentile_rank))%>%
#     mutate(pentile_rank = recode_factor(pentile_rank, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
#     # mutate(bin_full_space = as_factor(bin_full_space)) %>%
#     #Put these in ascending order regardless of when they show up in the data
#     # mutate(bin_full_space = fct_relevel(bin_full_space, sort(integer()))) %>%
#     mutate(nxt_pent = as_factor(nxt_pent)) %>%
#     # mutate(pentile_rank = fct_reorder(pentile_rank, integer(), .desc = T)) %>%
#     # mutate(nxt_pent = fct_reorder(nxt_pent, var, .desc = T)) %>%
#     mutate(nxt_x = Frames_post_treatment[row_number()+1]) %>%
#     mutate(diff_pent = as_factor(diff_pent)) %>%
#     mutate(diff_pent = recode_factor(diff_pent, '0' = "No Change", '1' = "Up 1", '2' = "Up 2", '3' = "Up 3", '4' = "Up 4", '5' = "Up 5", '-1' = "Down 1", '-2' = "Down 2", '-3' = "Down 3", '-4' = "Down 4", '-5' = "Down 5")) %>%
#     # mutate(diff_pent = fct_reorder(diff_pent, integer(), .desc = T)) %>%
#     # mutate(nxt_pent  = recode_factor(pentile_rank, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
#     # mutate(changes = map2(pentile_rank, nxt_pent, direction_sankey)) %>%
#     # mutate(changes = unlist(changes)) %>%
#     # mutate(changes = as_factor(changes)) %>%
#     ungroup()
#   # drop_na() #This is to get rid of the last row for each protein. As written this can be used for several time points instead of just the first and last
# 
#   return(sankeyed_data)
# }

sankeyed_data <- sankey_it(subset_controller(penetrance_values,0, 32) %>%
                             mutate(Time_post_treatment = Frames_post_treatment *7.5)  #. This would normally be present. If it is not, then make it again
                             # mutate(time_thresh = map2(Percentage_reloc_less, Frames_post_treatment, .f = yet_threshold))
                             )


#### GENERATE A PAIRED T-TEST BETWEEN THE STARTING AND ENDING FOR PROTEINS
sankeyed_data <- sankeyed_data %>% 
  mutate(nxt_thresh = as.character(nxt_thresh)) %>% 
  mutate(thresh_changes = as.character(thresh_changes)) %>% 
  group_by(Protein) %>% 
  mutate(thresh_changes = ifelse(is.na(nxt_thresh), Lag(thresh_changes,1),thresh_changes)) %>%
  ungroup() %>% 
  mutate(thresh_changes = as_factor(thresh_changes)) %>% 
  filter(thresh_changes %in% c('Above->Below', 'Below->Below', 'Below->Above', 'Above->Above')) # confirm that that all the proteins have the combination


without_stats_schemadist <- ggplot(data = sankeyed_data, aes(
  x = Frames_post_treatment,
  y = Percentage_reloc_less))+
  geom_sina(aes(color = thresh_changes, group = Frames_post_treatment), size = 2, alpha = 0.6, scale = 'count')+
  geom_hline(yintercept = 95)+
  geom_hline(yintercept = 20)+
  scale_color_npg()+
  scale_fill_npg()+
  theme_minimal()
ggsave(paste0(date, "_protiens_schema_distribution.pdf"), plot = without_stats_schemadist, width = 12, height = 9, device = cairo_pdf)
ggsave(paste0(date, "_protiens_schema_distribution.eps"), plot = without_stats_schemadist, width = 12, height = 9, device = cairo_pdf)


stats_schemadist <- ggstatsplot::ggwithinstats(data = sankeyed_data,
                           x = Frames_post_treatment,
                           y = Percentage_reloc_less,
                           type = 'np',
                           violin.args = list(width = 0),
                           point.path = FALSE,
                           # violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
                           point.args = list(size = 0))+
  geom_hline(yintercept = 95)+
  geom_hline(yintercept = 20)+
  geom_sina(data = sankeyed_data, mapping = aes(color = thresh_changes, group = Frames_post_treatment), size = 2, alpha = 0.6, scale = 'count')+
  scale_color_npg()+
  scale_fill_npg()
ggsave(paste0(date, "_protiens_schema_distribution_wStats.pdf"), plot = stats_schemadist, width = 12, height = 9, device = cairo_pdf)
ggsave(paste0(date, "_protiens_schema_distribution_wStats.eps"), plot = stats_schemadist, width = 12, height = 9, device = cairo_pdf)

  # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
  # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+

head(sankeyed_data) #just a quick look





bees_swarm <- ggplot(data = sankeyed_data,
       mapping = aes(x = Frames_post_treatment, y = Percentage_reloc_less, color = thresh_changes))+
  geom_beeswarm(method = 'hex', cex =1.2)+
  # geom_pointrange(mapping = aes(x = Frames_post_treatment, y = Percentage_reloc_less, group = Frames_post_treatment),
  #                 stat = "summary",
  #                 fun.min = function(z) { quantile(z,0.25)},
  #                 fun.max = function(z) { quantile(z,0.75)},
  #                 fun = median)+
  geom_boxplot(mapping = aes(group = Frames_post_treatment), width = 0.25, alpha = 0.2, outlier.shape = NA)+
  scale_color_nejm()+
  geom_hline(yintercept = 95)+
  geom_hline(yintercept = 20)+
  scale_y_continuous(breaks = c(0,20,40,60,80,100))+
  theme_minimal()+
  theme(legend.position = "bottom")
bees_swarm
  
ggsave(paste0(date, "_protiens_schema_distribution_BeeSwarm.pdf"), plot = bees_swarm, width = 12, height = 9, device = cairo_pdf)
ggsave(paste0(date, "_protiens_schema_distribution_BeeSwarm.eps"), plot = bees_swarm, width = 12, height = 9, device = cairo_pdf)


#, Funciton to generage requested Alluvium Flow with variable name
# protein_alluvium_pentile <- function(df, mod_name = ""){
#   date <- Sys.Date()
#   df <- na.omit(df,)
#   alluvium_plot <- ggplot(data = df %>%
#            distinct() %>%
#            mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
#          mapping = aes(x = Frames_post_treatment,
#                        stratum = pentile_rank,
#                        alluvium = Protein,
#                        fill = pentile_rank,
#                        label = pentile_rank))+
#     scale_fill_npg()+
#     geom_flow(stat = "alluvium", lode.guidance = "frontback", color = 'darkgray', width = 1/10)+
#     geom_stratum(width = .05)+
#     theme_sankey()+
#     theme(legend.position = "bottom")+
#     ggtitle('Global movment of protein through pentile ranks')
# 
#   ggsave(sprintf("%s_%s_pent_rank_Alluvium_proteins.png", date, mod_name), plot = alluvium_plot, width = 12, height = 6)
#   ggsave(sprintf("%s_%s_pent_rank_Alluvium_proteins.eps", date, mod_name), plot = alluvium_plot, width = 12, height = 6)
#   return(alluvium_plot)
# }

protein_alluvium_thresh <- function(df){
  alluvium_plot <- ggplot(data = df,
                          # distinct() %>%
                          # mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
                          mapping = aes(x = Frames_post_treatment,
                                        stratum = time_thresh,
                                        alluvium = Protein,
                                        label = time_thresh,
                                        fill = time_thresh))+
    geom_flow(aes(fill = thresh_changes), aes.bind = "flows")+
    geom_stratum(width = .05)+
    geom_text(stat = "stratum")+
    theme(legend.position = "bottom")
  #   ggtitle(sprintf('Global movment of protein through %s', var_name))
  # 
  ggsave(sprintf("%s_thresh_Alluvium_proteins.png", date), plot = alluvium_plot, width = 12, height = 6)
  ggsave(sprintf("%s_thresh_Alluvium_proteins.eps", date), plot = alluvium_plot, width = 12, height = 6)
  ggsave(sprintf("%s_thresh_Alluvium_proteins.pdf", date), plot = alluvium_plot, width = 12, height = 6)
}



protien_alluvium_subspace <- function(df){
  
    
    
    
  alluvium_plot <- ggplot(data = df %>% 
                            filter(thresh_changes == "Below->Below" | thresh_changes == "Below->Above") %>% 
                          mutate(sub_zero_use_space = fct_reorder(sub_zero_use_space, Percentage_reloc_less)),
                          mapping = aes(x = Frames_post_treatment,
                                        stratum = sub_zero_use_space,
                                        y = Percentage_reloc_less,
                                        alluvium = Protein,
                                        label = sub_zero_use_space,
                                        fill = sub_zero_use_space))+
    geom_flow(aes(fill = thresh_changes), aes.bind = "flows")+
    geom_stratum(width = .05, )+
    geom_text(stat = "stratum")+
    theme(legend.position = "bottom")
  alluvium_plot
  
testing_max <- read_parquet("D:/ALL_FINAL/Combined_by_perc/penetrance_updated_trimmed.parquet", as_data_frame = T) %>% 
  mutate(diff_pen = updated_yet_perc - Percentage_reloc) %>% 
  arrange(diff_pen) %>% 

View(testing_max)  col
  # 
  # ggplot(data = df %>% 
  #          filter(thresh_changes == "Below->Below" | 	
  #                   thresh_changes == "Below->Above")
  #        mapping = aes(x = Frames_post_treatment,
  #                      stratum = sub_zero_use_space,
  #                      alluvium = Protein,
  #                      fill = sub_zero_use_space,
  #                      label = sub_zero_use_space))+
  #   scale_fill_npg()+
  #   geom_flow(stat = "alluvium", lode_forward())+
  #   geom_stratum(width = .05)+
  #   theme_sankey()+
  #   theme(legend.position = "bottom")

}





protein_alluvium_thresh_full <- function(df, target_si, target_ei){
  
  yet_threshold <- function(val, fr_post, target_s = target_si, target_e = target_ei){
    if (fr_post <= 0){
      if (val >= target_s){
        return("Above")
      }
      if (val < target_s){
        return("Below")
      }
    }
    if (fr_post > 0){
      if (val >= target_e){
        return("Above")
      }
      if (val < target_e){
        return("Below")
      }
    }
  }
  
  
  sankey_it <- function(df){
    sankeyed_data <- df%>%
      mutate(Percentage_reloc_less = as.integer(Percentage_reloc_less)) %>%
      group_by(Protein) %>%
      arrange(Frames_post_treatment, .by_group = T) %>%
      
      
      #* Pure threshold by set values. Likely set to >=15% for priming and < 95 for incomplete penetrance
      #* This was added by request from Grant for simplification of calls for schema grouping
      mutate(time_thresh = map2(Percentage_reloc_less, Frames_post_treatment, .f = yet_threshold)) %>%
      mutate(time_thresh = unlist(time_thresh)) %>% 
      mutate(nxt_thresh = time_thresh[row_number()+1]) %>% 
      mutate(thresh_changes = map2(time_thresh, nxt_thresh, direction_sankey)) %>%
      mutate(thresh_changes = unlist(thresh_changes)) %>% 
      mutate(thresh_changes = as_factor(thresh_changes)) %>%
      mutate(time_thresh = as_factor(time_thresh)) %>% 
      
      
      mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
      
      #* Seperate bins into High, Medium, Low. This only work if divisible by 3 for equal size. 
      #* Moving to 5 will mean decision must be made about whether 3 is the middle or 2-4 is
      mutate(three_bin = map(pentile_rank, three_bin)) %>%
      mutate(nxt_three = three_bin[row_number()+1]) %>%
      
      #* Generate the next equal member rank
      mutate(nxt_pent = pentile_rank[row_number()+1]) %>%
      mutate(diff_pent = pentile_rank - nxt_pent) %>%
      mutate(diff_pent = replace_na(diff_pent, 0)) %>%
      mutate(diff_pent_dup = diff_pent) %>%
      mutate(pentile_cpy = pentile_rank) %>%
      mutate(pentile_rank = as_factor(pentile_rank))%>%
      mutate(pentile_rank = recode_factor(pentile_rank, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
      mutate(diff_pent = as_factor(diff_pent)) %>%
      mutate(diff_pent = recode_factor(diff_pent, '0' = "No Change", '1' = "Up 1", '2' = "Up 2", '3' = "Up 3", '4' = "Up 4", '5' = "Up 5", '-1' = "Down 1", '-2' = "Down 2", '-3' = "Down 3", '-4' = "Down 4", '-5' = "Down 5")) %>%
      mutate(nxt_pent = as_factor(nxt_pent)) %>%
      mutate(nxt_pent  = recode_factor(nxt_pent, '1' = 'Top', '2' = 'Second', '3' = 'Third', '4' = 'Fourth', '5' = 'Fifth', '6' = 'Sixth')) %>%
      mutate(Changes_full = map2(pentile_rank, nxt_pent, direction_sankey)) %>%
      mutate(Changes_full = unlist(Changes_full)) %>%
      mutate(Changes_full = as_factor(Changes_full)) %>%
      
      #* Generate the equal space rank
      mutate(nxt_decSp = bin_full_space[row_number()+1]) %>%
      #Get the difference between the prior frame present (there can be spaced for subseting) and the current
      mutate(diff_decSp = bin_full_space - nxt_decSp) %>%
      mutate(diff_decSp = replace_na(diff_decSp, 0)) %>% #This is to deal with the final value where there is nothing to compare against
      mutate(bin_full_space = as_factor(bin_full_space)) %>%
      
      
      #Generate the next used space rank
      mutate(nxt_useSp = bin_use_space[row_number()+1]) %>% 
      mutate(diff_useSp = as.integer(bin_use_space) - as.integer(nxt_useSp)) %>%
      mutate(diff_useSp = replace_na(diff_useSp, 0)) %>%
      
      
      #* Get the nxt_x values asuming that Frames_post_treatment is on the x
      mutate(nxt_x = Frames_post_treatment[row_number()+1]) %>%
      
      
      mutate(Changes_three = map2(three_bin, nxt_three, direction_sankey)) %>%
      mutate(Changes_three = unlist(Changes_three)) %>%
      mutate(Changes_three = as_factor(Changes_three)) %>%
      
      
      #Ungroup the dataframe
      ungroup()
    # mutate(diff_pent = fct_reorder(diff_pent, diff_pent_dup, .desc = T)) %>%
    # drop_na() #This is to get rid of the last row for each protein. As written this can be used for several time points instead of just the first and last
    
    return(sankeyed_data)
  }

  df = sankey_it(subset_controller(penetrance_values, 31) %>%
                   mutate(Time_post_treatment = Frames_post_treatment *7.5)  #. This would normally be present. If it is not, then make it again
                 # mutate(time_thresh = map2(Percentage_reloc_less, Frames_post_treatment, .f = yet_threshold))
  )
  date <- Sys.Date()
  alluvium_plot <- ggplot(data = df,
                            # distinct() %>%
                            # mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
                          mapping = aes(x = Frames_post_treatment,
                                        stratum = time_thresh,
                                        alluvium = Protein,
                                        label = time_thresh,
                                        fill = time_thresh))+
    geom_flow(aes(fill = thresh_changes), aes.bind = "flows")+
    geom_stratum(width = .05)+
    geom_text(stat = "stratum")+
    theme(legend.position = "bottom")
  #   ggtitle(sprintf('Global movment of protein through %s', var_name))
  # 
  ggsave(sprintf("%s_thresh(%s,%s)_Alluvium_proteins.png", date, target_si, target_ei), plot = alluvium_plot, width = 12, height = 6)
  ggsave(sprintf("%s_thresh(%s,%s)_Alluvium_proteins.eps", date, target_si, target_ei), plot = alluvium_plot, width = 12, height = 6)
  ggsave(sprintf("%s_thresh(%s,%s)_Alluvium_proteins.pdf", date, target_si, target_ei), plot = alluvium_plot, width = 12, height = 6)
  
  
  test <- sankeyed_data %>% 
    mutate(nxt_thresh = as.character(nxt_thresh)) %>% 
    mutate(thresh_changes = as.character(thresh_changes)) %>% 
    group_by(Protein) %>% 
    mutate(thresh_changes = ifelse(is.na(nxt_thresh), Lag(thresh_changes,1),thresh_changes)) %>%
    ungroup() %>% 
    mutate(thresh_changes = as_factor(thresh_changes)) %>% 
    filter(thresh_changes %in% c('Above->Below', 'Below->Below', 'Below->Above', 'Above->Above'))
  
  
  distribution <- ggplot(data = test, aes(
    x = Frames_post_treatment,
    y = Percentage_reloc_less))+
    geom_sina(aes(color = thresh_changes, group = Frames_post_treatment), size = 2, alpha = 0.6, scale = 'width')+
    scale_color_npg()+
    scale_fill_npg()+
    theme_minimal()
  
  
  ggsave(sprintf("%s_thresh(%s,%s)_Distribution_sankey_proteins.png", date, target_si, target_ei), plot = distribution, width = 12, height = 6)
  ggsave(sprintf("%s_thresh(%s,%s)_Distribution_sankey_proteins.eps", date, target_si, target_ei), plot = distribution, width = 12, height = 6)
  ggsave(sprintf("%s_thresh(%s,%s)_Distribution_sankey_proteins.pdf", date, target_si, target_ei), plot = distribution, width = 12, height = 6)
  
  return(distribution)
}



protein_alluvium_var <- function(df, var_name){
  date <- Sys.Date()
  alluvium_plot <- ggplot(data = df %>%
    
                            distinct() %>%
                            mutate(Protein = fct_reorder(Protein, Percentage_reloc_less)),
                          mapping = aes(x = Frames_post_treatment,
                                        stratum = .data[[var_name]],
                                        
                                        alluvium = Protein,
                                        fill = .data[[var_name]],
                                        label = .data[[var_name]]))+
    scale_fill_npg()+
    geom_flow(stat = "alluvium", lode.guidance = "frontback")+
    geom_stratum(width = .05)+
    theme_sankey()+
    theme(legend.position = "bottom")+
    ggtitle(sprintf('Global movment of protein through %s', var_name))

  ggsave(sprintf("%s_%s_pent_rank_Alluvium_proteins.png", date, var_name), plot = alluvium_plot, width = 12, height = 6)
  ggsave(sprintf("%s_%s_pent_rank_Alluvium_proteins.pdf", date, var_name), plot = alluvium_plot, width = 12, height = 6)
  return(alluvium_plot)
}

# sankeyed_data %>%
#   map(.x = ., .y = c('pentile_rank', 'bin_full_space'),
#        ~protein_alluvium(.x, .y)
#   )

protein_alluvium_thresh(sankeyed_data)


#### Testing
# alluvium_plot <-
  
df <- sankeyed_data %>% 
  mutate(Frames_post_treatment = as.integer(as.character(Frames_post_treatment))) %>%
  mutate(quick_group = case_when((Frames_post_treatment <= 0 & time_thresh == "Below") ~ "B",
                                 (Frames_post_treatment <= 0 & time_thresh == "Above") ~ "A",
                                 (Frames_post_treatment > 0 ) ~ "L")) %>%
  group_by(quick_group) %>%
  mutate(sub_zero_use_space = cut_interval(Percentage_reloc_less, n = 5, labels = c('Fifth','Fourth','Third','Second','Top'))) %>%
  mutate(zero_pentile = ntile(Percentage_reloc_less,5)) %>%
  mutate(zero_pentile = 6 - zero_pentile) %>%
  mutate(zero_pentile = paste0(as.character(zero_pentile),quick_group)) %>%
  select(Protein, Frames_post_treatment, Percentage_reloc_less, thresh_changes, time_thresh, zero_pentile, sub_zero_use_space, bin_use_space, bin_full_space, nxt_useSp) %>%
  ungroup() %>%
  mutate(zero_pentile = as_factor(zero_pentile)) %>% 
  mutate(zero_pentile = fct_reorder(zero_pentile, Percentage_reloc_less)) %>% 
  mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>% 
  mutate(sub_zero_use_space = as_factor(sub_zero_use_space)) %>% 
  mutate(Frames_post_treatment = recode_factor(Frames_post_treatment, '0' = 'Start', '32' = 'End'))#'31' = 'End'))


diff_slice <- function(series,quick_group){
  if (quick_group[1] == 'L'){
    use_space = cut_interval(series, n = 10, labels = c('Tenth', 'Ninth', 'Eighth', 'Seventh', 'Sixth', 'Fifth','Fourth','Third','Second','Top'))
  } else{
    use_space = cut_interval(series, n = 5, labels = c('Fifth','Fourth','Third','Second','Top'))
  }
  return(use_space)
}

diff_slice_top <- function(series, top_group){
  if (top_group[1] == 'Y'){
    use_space = cut_interval(series, n = 4, labels = c('Top-D', 'Top-C', 'Top-B', 'Top-A'))
  } else if (top_group[1] == 'N'){
    use_space = cut_interval(series, n = 4, labels = c('Fifth','Fourth','Third','Second'))
  } else {
    use_space = cut_interval(series, n = 5, labels = c('Fifth','Fourth','Third','Second','Top'))
  }
  return(use_space)
}

df<- df %>% 
  mutate(top_group = case_when((quick_group == 'L' & sub_zero_use_space == "Top") ~ "Y",
                               (quick_group == 'L' & sub_zero_use_space != "Top") ~ "N",
                               (quick_group == 'A') ~ "A",
                               (quick_group == 'B') ~ "B")) %>%
  mutate(quick_cpy = top_group) %>% 
  mutate(quick_cpy = as.character(quick_cpy)) %>% 
  group_by(top_group) %>% 
  mutate(sub_zero_use_space = diff_slice_top(Percentage_reloc_less,quick_cpy)) %>%
  ungroup() %>% 
  mutate(Frames_post_treatment = fct_relevel(Frames_post_treatment,"Start", "End")) %>% 
  # mutate(sub_zero_use_space = fct_reorder(sub_zero_use_space, Percentage_reloc_less)) %>% 
  mutate(sub_zero_use_space = fct_relevel(sub_zero_use_space, 'Fifth', 'Fourth', 'Third', 'Second', 'Top', 'Top-D', 'Top-C', 'Top-B', 'Top-A'))
  
write_csv(df, "sankey_data_move_TopSplit32.csv")
write_parquet(df, "sankey_data_move_32.parquet")


df <- read_csv("sankey_data_move.csv") %>% 
  mutate(zero_pentile = as_factor(zero_pentile)) %>% 
  mutate(zero_pentile = fct_reorder(zero_pentile, Percentage_reloc_less)) %>% 
  mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>% 
  mutate(sub_zero_use_space = as_factor(sub_zero_use_space)) %>% 
  mutate(Frames_post_treatment = recode_factor(Frames_post_treatment, '0' = 'Start', '31' = 'End'))

  
  # filter(Frames_post_treatment == 0)

# df_alu <- to_alluvia_form(df, key = 'Frames_post_treatment', value = "sub_zero_use_space", id = "Protein")  

# is_alluvia_form(df_alu, key = "Frames_post_treatment", value = "sub_zero_use_space", id = "Protein")


# ggplot(df,
#        aes(x = Frames_post_treatment,
#            stratum = sub_zero_use_space,
#            alluvium = Protein)) +
#   geom_alluvium(aes(fill = thresh_changes)) +
#   geom_stratum()+
#   stat_stratum(reverse = FALSE)+
#   stat_stratum(geom = 'text', aes(label = sub_zero_use_space), reverse = FALSE)

#,This looks like it is unused
df <- read_csv_arrow("sankey_data_move.csv") %>%
  mutate(Frames_post_treatment = fct_relevel(Frames_post_treatment,"Start", "End")) %>%
  # mutate(sub_zero_use_space = fct_reorder(sub_zero_use_space, Percentage_reloc_less)) %>%
  mutate(sub_zero_use_space = fct_relevel(sub_zero_use_space, 'Tenth', 'Ninth', 'Eighth', 'Seventh', 'Sixth', 'Fifth', 'Fourth', 'Third', 'Second', 'Top')) %>%
  mutate(zero_pentile = fct_reorder(zero_pentile, Percentage_reloc_less))


#Added After July2024 Meeting
#Both bins to final
df <- read.csv("sankey_data_move_TopSplit.csv") %>%
  mutate(zero_pentile = as_factor(zero_pentile)) %>% 
  mutate(zero_pentile = fct_reorder(zero_pentile, Percentage_reloc_less)) %>% 
  mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>% 
  mutate(sub_zero_use_space = as_factor(sub_zero_use_space)) %>% 
  mutate(Frames_post_treatment = recode_factor(Frames_post_treatment, '0' = 'Start', '31' = 'End'))


f_choose <- function(sel_col, s){
  if (sel_col == 'L'){
    return(s)
  }else if (sel_col == 'A'){
    return(sel_col)
  }else if (sel_col == 'B'){
    return(sel_col)
  }
}

df_simpStart <- df %>% 
  mutate(sub_zero_use_space = as.character(sub_zero_use_space)) %>% 
  mutate(simpFlow = map2_chr(quick_group, sub_zero_use_space, f_choose)) %>% 
  mutate(simpFlow = as_factor(simpFlow)) %>%
  mutate(simpFlow = fct_reorder(simpFlow, Percentage_reloc_less)) %>% 
  mutate(sub_zero_use_space = as_factor(sub_zero_use_space)) %>% 
  mutate(sub_zero_use_space = fct_relevel(sub_zero_use_space, 'Fifth', 'Fourth', 'Third', 'Second', 'Top', 'Top-D', 'Top-C', 'Top-B', 'Top-A'))
  

simpFlows <- ggplot(df_simpStart,
                        aes(x = Frames_post_treatment,
                            stratum = simpFlow,
                            alluvium = Protein,
                            label = simpFlow,
                            fill = simpFlow,
                            y = Percentage_reloc_less))+
  geom_flow(aes(fill = simpFlow), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = simpFlow), reverse = FALSE)+
  theme(legend.position = "bottom")+
  scale_color_aaas()
ggsave(sprintf("%s_SubAbove_use_Alluvium_proteins.png", date), plot = sub_above_use, width = 12, height = 6)
ggsave(sprintf("%s_SubAbove_use_Alluvium_proteins.pdf", date), plot = sub_above_use, width = 12, height = 6)


#####

sub_above_use <- ggplot(df %>% filter(thresh_changes %in% c('Above->Above', 'Above->Below')),
       aes(x = Frames_post_treatment,
           stratum = sub_zero_use_space,
           alluvium = Protein,
           label = sub_zero_use_space,
           fill = sub_zero_use_space))+
  geom_flow(aes(fill = sub_zero_use_space), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = sub_zero_use_space), reverse = FALSE)+
  theme(legend.position = "bottom")+
  scale_color_aaas()
ggsave(sprintf("%s_SubAbove_use_Alluvium_proteins.png", date), plot = sub_above_use, width = 12, height = 6)
ggsave(sprintf("%s_SubAbove_use_Alluvium_proteins.pdf", date), plot = sub_above_use, width = 12, height = 6)

sub_above_pent<- ggplot(df %>% filter(thresh_changes %in% c('Above->Above', 'Above->Below')),
       aes(x = Frames_post_treatment,
           stratum = zero_pentile,
           alluvium = Protein,
           label = zero_pentile,
           fill = zero_pentile))+
  geom_flow(aes(fill = zero_pentile), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = zero_pentile), reverse = FALSE)+
  theme(legend.position = "bottom")
ggsave(sprintf("%s_SubAbove_pent_Alluvium_proteins.png", date), plot = sub_above_pent, width = 12, height = 6)
ggsave(sprintf("%s_SubAbove_pent_Alluvium_proteins.pdf", date), plot = sub_above_pent, width = 12, height = 6)


sub_below_use <- ggplot(df %>% filter(thresh_changes %in% c('Below->Above', 'Below->Below')),
                        aes(x = Frames_post_treatment,
                            stratum = sub_zero_use_space,
                            alluvium = Protein,
                            label = sub_zero_use_space,
                            fill = sub_zero_use_space))+
  geom_flow(aes(fill = sub_zero_use_space), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = sub_zero_use_space), reverse = FALSE)+
  theme(legend.position = "bottom")
ggsave(sprintf("%s_SubBelow_use_Alluvium_proteins.png", date), plot = sub_below_use, width = 12, height = 6)
ggsave(sprintf("%s_SubBelow_use_Alluvium_proteins.pdf", date), plot = sub_below_use, width = 12, height = 6)


sub_below_pent<- ggplot(df %>% filter(thresh_changes %in% c('Below->Above', 'Below->Below')),
                        aes(x = Frames_post_treatment,
                            stratum = zero_pentile,
                            alluvium = Protein,
                            label = zero_pentile,
                            fill = zero_pentile))+
  geom_flow(aes(fill = zero_pentile), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = zero_pentile), reverse = FALSE)+
  theme(legend.position = "bottom")
ggsave(sprintf("%s_SubBelow_pent_Alluvium_proteins.png", date), plot = sub_below_pent, width = 12, height = 6)
ggsave(sprintf("%s_SubBelow_pent_Alluvium_proteins.pdf", date), plot = sub_below_pent, width = 12, height = 6)



############ This plot is not showing the final in the correct order.
df <- df %>% 
  mutate(sub_zero_use_space = paste0(quick_group," - ", as.character(sub_zero_use_space))) %>% 
  mutate(sub_zero_use_space = as_factor(sub_zero_use_space)) %>% 
  mutate(sub_zero_use_space = fct_relevel(sub_zero_use_space, 'L - Fifth', 'L - Fourth', 'L - Third', 'L - Second', 'L - Top-D', 'L - Top-C', 'L - Top-B', 'L - Top-A', 'B - Fifth', 'B - Fourth', 'B - Third', 'B - Second', 'B - Top', 'A - Fifth', 'A - Fourth', 'A - Third', 'A - Second', 'A - Top'))

sub_g_use <- ggplot(df,
                        aes(x = Frames_post_treatment,
                            stratum = sub_zero_use_space,
                            alluvium = Protein,
                            label = sub_zero_use_space,
                            fill = sub_zero_use_space))+
  geom_flow(aes(fill = sub_zero_use_space), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = sub_zero_use_space), reverse = FALSE)+
  theme(legend.position = "bottom")
ggsave(sprintf("%s_SubALL_use_Alluvium_proteins.png", date), plot = sub_g_use, width = 6, height = 12)
ggsave(sprintf("%s_SubALL_use_Alluvium_proteins.pdf", date), plot = sub_g_use, width = 6, height = 12)

sub_g_use_perc <- ggplot(df,
                    aes(x = Frames_post_treatment,
                        stratum = sub_zero_use_space,
                        alluvium = Protein,
                        label = sub_zero_use_space,
                        fill = sub_zero_use_space,
                        y = Percentage_reloc_less))+
  geom_flow(aes(fill = sub_zero_use_space), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = sub_zero_use_space), reverse = FALSE)+
  theme(legend.position = "bottom")
ggsave(sprintf("%s_SubALL_use_perc_Alluvium_proteins.png", date), plot = sub_g_use_perc, width = 6, height = 12)
ggsave(sprintf("%s_SubALL_use_perc_Alluvium_proteins.pdf", date), plot = sub_g_use_perc, width = 6, height = 12)


###############

sub_g_pent<- ggplot(df,
                        aes(x = Frames_post_treatment,
                            stratum = zero_pentile,
                            alluvium = Protein,
                            label = zero_pentile,
                            fill = zero_pentile))+
  geom_flow(aes(fill = zero_pentile), aes.bind = "flows", reverse = FALSE) +
  geom_stratum(width = 0.25, reverse = FALSE)+
  # geom_text(stat = "stratum", reverse = FALSE)+
  stat_stratum(reverse = FALSE)+
  stat_stratum(geom = 'text', aes(label = zero_pentile), reverse = FALSE)+
  theme(legend.position = "bottom")
ggsave(sprintf("%s_SubALL_pent_Alluvium_proteins.png", date), plot = sub_g_pent, width = 12, height = 6)
ggsave(sprintf("%s_SubALL_pent_Alluvium_proteins.pdf", date), plot = sub_g_pent, width = 12, height = 6)


# ggplot(df_alu, aes(x = Frames_post_treatment,
#                stratum = sub_zero_use_space,
#                alluvium = thresh_changes,
#                fill = sub_zero_use_space,
#                label = sub_zero_use_space))+
#   scale_fill_aaas()+
#   geom_flow(stat = 'alluvium', lode.guidance = 'frontback',
#             color = 'darkgray')+
#   geom_stratum()+
#   theme(legend.position = 'bottom')




x <- to_alluvia_form(df,
              key = "Frames_post_treatment", value = "sub_zero_use_space", id = "thresh_changes")



#+ 
  # Vanilla GGplot here onwards
  ggtitle("Majors opted across semesters")+
  scale_y_discrete() +
  ylab("Number of students enrolled")+
  theme_bw()+
  theme(axis.text = element_text(size = 7))

# 
# ggplot(data = df %>% 
#          filter(time_thresh == "Below") %>% 
#          mutate(sub_zero_use_space = fct_reorder(sub_zero_use_space, Percentage_reloc_less)),
#                         mapping = aes(y = Percentage_reloc_less,
#                                       axis1 = sub_zero_use_space,
#                                       axis2 = nxt_useSp))+
#   geom_alluvium(aes(fill = thresh_changes), aes.bind=TRUE, width = 1/12)+
#   geom_stratum(width = 1/4, fill = "white", color = "black")+
#   geom_text(stat = "stratum", label.strata = TRUE) +
#   scale_x_discrete(limits = c("Time threshold scaled space", "Full library use space"),
#                    expand = c(.05, .05)) +
#   # scale_fill_manual(values = c("red", "orange", "blue")) +
#   # labs(y = "Percentage Reloc Less") +
#   theme_minimal()+
#   theme(legend.position = "bottom")
# 
# 
# 
# 
#   
# ggplot(aes(y = Percentage_reloc_less, axis1 = origin, axis2 = carrier, axis3 = dest)) +
#   geom_alluvium(aes(fill = origin), aes.bind=TRUE, width = 1/12) +
#   geom_stratum(width = 1/4, fill = "white", color = "black") +
#   geom_text(stat = "stratum", label.strata = TRUE) +
#   scale_x_discrete(limits = c("Origin", "Carrier", "Destination"),
#                    expand = c(.05, .05)) +
#   scale_fill_manual(values = c("red", "orange", "blue")) +
#   labs(y = "Cases") +
#   theme_minimal() +
#   theme(legend.position = "none") +
#   ggtitle("NYC flights volume for top destinations and airlines")  
# 
# 
# 
# ggplot(as.data.frame(UCBAdmissions),
#        aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
#   geom_alluvium(aes(fill = Admit), width = 1/12) +
#   geom_stratum(width = 1/12, fill = "black", color = "grey") +
#   geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.05, .05)) +
#   scale_fill_brewer(type = "qual", palette = "Set1") +
#   ggtitle("UC Berkeley admissions and rejections, by sex and department")
# 
# 


  ggtitle(sprintf('Global movment of protein through %s', b))



protein_alluvium_var(sankeyed_data %>% filter(time_thresh == "Below", Frames_post_treatment == '-1'), 'bin_use_space')


protein_alluvium_var(sankeyed_data, 'time_thresh')
protein_alluvium_var(sankeyed_data, 'pentile_rank')
protein_alluvium_var(sankeyed_data, 'bin_full_space')


thresh_member_table <- function(df){
  df <- df %>%
    mutate(Frames_post_treatment = as.numeric(as.character(Frames_post_treatment))) %>% 
    filter(Frames_post_treatment <= 0) %>%
    select(Protein, thresh_changes) %>% 
    drop_na() %>%
    as.data.table()
  # filter((diff_pent_dup > 0 )| (diff_pent_dup < 0))
  # mutate(pentile_rank = fct_relevel(pentile_rank, sort(integer())))#Remove the proteins which did not move
  table <- df[, .("Protein Members" = list(unlist(unique(Protein))), N = .N), by = thresh_changes] %>%
    arrange(thresh_changes) %>%
    # arrange(Changes_full) %>% 
    gt()
  
  gtsave(data = table, filename = sprintf("%s_%s_thresh_table.html", date, by))
  gtsave(data = table, filename =  sprintf("%s_%s_thresh_table.pdf", date, by))
  
  # row_group_order(groups = levels(df$pentile_rank))
  return(table)
}

df <- test %>%
  mutate(Frames_post_treatment = as.numeric(as.character(Frames_post_treatment))) %>% 
  filter(Frames_post_treatment <= 0) %>% 
  select(Protein, thresh_changes) %>% 
  drop_na() %>%
  as.data.table()
# filter((diff_pent_dup > 0 )| (diff_pent_dup < 0))
# mutate(pentile_rank = fct_relevel(pentile_rank, sort(integer())))#Remove the proteins which did not move
table <- df[, .("Protein Members" = list(unlist(unique(Protein))), N = .N), by = thresh_changes] %>%
  arrange(thresh_changes) %>%
  # arrange(Changes_full) %>% 
  gt()
table
by= 'priming_confluecy'
gtsave(data = table, filename = sprintf("%s_%s_thresh_table.html", date, by))
gtsave(data = table, filename =  sprintf("%s_%s_thresh_table.pdf", date, by))



get_member_table <- function(df, by = 'pentile_rank'){
  df <- df %>%
    filter(!(Protein %like% 'DDC2')) %>%
    filter(!(Protein %like% 'd0210')) %>%
    # filter((as.integer(diff_pent) < -1 )| (as.integer(diff_pent) > 1)) %>%
    filter(Frames_post_treatment == '0') %>%
    drop_na() %>%
    mutate(Changes_three = fct_relevel(Changes_three, "High->High", "High->Mid", "High->Low", "Mid->High", "Mid->Mid", "Mid->Low", "Low->High", "Low->Mid", "Low->Low")) %>% 
    mutate(Changes_full = fct_reorder(Changes_full, pentile_cpy)) %>%
    as.data.table()
    # filter((diff_pent_dup > 0 )| (diff_pent_dup < 0))
    # mutate(pentile_rank = fct_relevel(pentile_rank, sort(integer())))#Remove the proteins which did not move
  table <- df[, .("Protein Members" = list(unlist(unique(Protein))), N = .N), by = by] %>%
    # arrange(Changes_three) %>%
    # arrange(Changes_full) %>% 
    gt()
    

  gtsave(data = table, filename = sprintf("%s_%s_penHet_table.html", date, by))
  gtsave(data = table, filename =  sprintf("%s_%s_penHet_table.pdf", date, by))

    # row_group_order(groups = levels(df$pentile_rank))
  return(table)
}

  
thresh_member_table(test)
thresh_member_table(sankeyed_data)

get_member_table(sankeyed_data, by = 'diff_pent')
get_member_table(sankeyed_data, by = 'bin_full_space')
get_member_table(sankeyed_data, by = 'Changes_full')
get_member_table(sankeyed_data, by = 'Changes_three')



#testing again
df <- sankeyed_data %>%
  filter(!(Protein %like% 'DDC2')) %>%
  filter(!(Protein %like% 'd0210')) %>%
  # filter((as.integer(diff_pent) < -1 )| (as.integer(diff_pent) > 1)) %>%
  filter(Frames_post_treatment == '0') %>%
  mutate(changes = fct_reorder(changes, pentile_cpy)) %>%
  drop_na() %>%
  as.data.table()


is_lodes_form(df,)


sankey_thresh_shift <- ggplot(data = sankeyed_data %>% ungroup()%>% make_long(time_thresh, thresh_changes),
                       mapping = aes(x = Frames_post_treatment,
                                     next_x = nxt_x,
                                     node = time_thresh,
                                     next_node = nxt_thresh,
                                     group = Protein,
                                     fill = factor(time_thresh)))+
  geom_sankey()+
  geom_sankey_label()+
  theme_sankey()
sankey_thresh_shift





sankey_shift <- ggplot(data = sankeyed_data %>% make_long(pentile_rank, diff_pent),
       mapping = aes(x = Frames_post_treatment,
                     next_x = nxt_x,
                     node = pentile_rank,
                     next_node = nxt_pent,
                     group = Protein,
                     fill = factor(pentile_rank)))+
  geom_sankey()+
  geom_sankey_label()+
  theme_sankey()
sankey_shift

sankey_sub <- sankeyed_data %>%
  filter((diff_pent_dup > 0 )| (diff_pent_dup < 0)) #Remove the proteins which did not move


get_types <- function(df){
  bet_het_n <- df %>%
    filter(diff_pent_dup < 0) %>%
    nrow()
  bet_het_members <- df %>%
    filter(diff_pent_dup < 0) %>%
    unique(.$Protein)
  #   reframe(N = nrow(.), Proteins = unique(.$Protein))
  # misc_het <- df %>%
  #   filter(diff_pent_dup > 0) %>%
  #   reframe(N = nrow(.), Proteins = unique(.$Protein))
  # homogenous <- df %>%
  #   filter(diff_pent_dup == 0) %>%
  #   reframe(N = nrow(.), Proteins = unique(.$Protein))

  return(table_new)
}

x <- get_types(sankeyed_data)
sankeyed_data %>%
  filter(diff_pent_dup < 0) %>%
  distinct(Protein) %>%
  as.list()

table_results <- gt(sankeyed_data) %>%
  tab_header(
    title = "Assigned Schemas",
    subtitle = "Movement of proteins through pentile ranks"
  ) %>%
  tab_source_note(
    source_note = "Data from the penetrance values non-shifters removed"
  ) %>%
  tab_stubhead(label = "Relocation Schema") %>%
  tab_row_group(label = 'Homogenous',
                rows = matches(filter(sankeyed_data, diff_pent_dup < 0)$Protein)) #%>%
  # tab_row_group(label = 'Bet-hedging Relocation',
  #               rows = slice(filter(sankeyed_data, diff_pent_dup < 0))) %>%
  # tab_row_group(label = 'Miscellaneous Relocation',
  #               rows = slice(filter(sankeyed_data, diff_pent_dup > 0)))
  # tab_spanner(
  #   label = "Description of Members",
  #   columns = vars(Protein),
  #   locations = cells_body(rows())
  # )




long_sankey <- sankey_sub %>%   #sankeyed_data %>%
  ungroup() %>%
  filter(Frames_post_treatment == '0')%>%
  make_long(., pentile_rank, diff_pent) %>%
  drop_na(node) %>%
  mutate(node = as_factor(node))

long_sankey$node <- factor(long_sankey$node,
                           levels=c('Sixth','Fifth','Fourth','Third',
                                    'Second','Top', "Down 5", "Down 4",
                                    "Down 3", "Down 2", "Down 1", "Up 1",
                                    "Up 2", "Up 3", "Up 4", "Up 5")) # Got rid of the No Change for subset
                           # levels=c('Sixth','Fifth','Fourth','Third','Second','First', "Down 5", "Down 4", "Down 3", "Down 2", "Down 1", "No Change", "Up 1", "Up 2", "Up 3", "Up 4", "Up 5"))

sankey_shift <- ggplot(data = long_sankey,
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
sankey_shift

ggsave(sprintf('%s_Sankey_shift.png', date), plot = sankey_shift, width = 12, height = 6)




df <- mtcars %>%
  make_long(cyl, vs, am, gear, carb)



penetrance_values$Frames_post_treatment <- as_factor(penetrance_values$Frames_post_treatment)


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