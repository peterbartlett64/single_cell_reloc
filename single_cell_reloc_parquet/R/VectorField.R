# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-06-07
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
# library(ggstatsplot)
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
# library(statsExpressions)
library(purrr)
library(ggthemr)
library(extrafont)
library(grid)
library(gganimate)
library(rsample)
library(tidyr)
library(rlang)
library(forcats)
library(ggbeeswarm)
library(ggside)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

### Load Dataset
#This is just a subset of the data so that I can see if this actually works
#, This is the slowest read in

#####
# df_for_facet_abundance = read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>%
#   select(LogAbundance, Loc_score, z_score_Loc, z_score_logAbund, Protein, Frames_post_treatment, Cell_Barcode, Unique_pos, Relocalized, Yet) %>%
#   # filter(Frames_post_treatment %in% c(0,8,16,24,32)) %>%
#   filter(Protein %in% c('SAE2', 'RAD51', 'RMI1')) %>%
#   filter(Frames_post_treatment >= 0, Frames_post_treatment <= 32) %>%
#    #* This could be changed to have more than just the post
#   collect()
# 
# 
# df <- df_for_facet_abundance %>%
#   filter(Protein == 'RAD51') %>%
#   unique() %>%
#   group_by(Cell_Barcode) %>%
#   arrange(Frames_post_treatment) %>%
#   mutate(v_Loc = Loc_score[row_number()+1]- Loc_score[row_number()]) %>%
# 
#   mutate(v_z_Loc = z_score_Loc[row_number()+1]- z_score_Loc[row_number()]) %>% #Velocity
#   mutate(a_z_Loc = v_z_Loc[row_number()+1]- v_z_Loc[row_number()]) %>% #acceleration
#   
#   mutate(v_Abund = as.numeric(LogAbundance[row_number()+1]- LogAbundance[row_number()])) %>%
#   mutate(v_z_Abund = as.numeric(z_score_logAbund[row_number()+1]- z_score_logAbund[row_number()])) %>%
#   mutate(a_z_Abund = as.numeric(v_z_Abund[row_number()+1]- v_z_Abund[row_number()])) %>%
# 
#   mutate(hyp = as.numeric(sqrt((v_z_Loc)^2 + (v_z_Abund)^2))) %>% #hyp for scaling
# 
#   mutate(v_y_sc_hyp = as.numeric(v_z_Abund/hyp)) %>% #normalize by hyp to get vectors of the same size
#   mutate(v_x_sc_hyp = as.numeric(v_z_Loc/hyp)) %>%
#   ungroup() %>%
#   drop_na()
# 
# norm = FALSE
# if (norm){
#     x_axis_s <- 'z_score_Loc'
#     y_axis_s <- 'z_score_logAbund'
#     vx_axis_s <- 'v_z_Loc'
#     vy_axis_s <- 'v_z_Abund'
#   } else{
#     x_axis_s <- 'Loc_score'
#     y_axis_s <- 'LogAbundance'
#     vx_axis_s <- 'v_Loc'
#     vy_axis_s <- 'v_Abund'
#   }
#   
#   x_axis <- rlang::ensym(x_axis_s)
#   y_axis <- rlang::ensym(y_axis_s)
#   vx_axis <- rlang::ensym(vx_axis_s)
#   vy_axis <- rlang::ensym(vy_axis_s)
#   
#   df <- df_input %>%
#     filter(Protein == 'RAD51') %>%
#     unique() %>%
#     group_by(Cell_Barcode) %>%
#     arrange(Frames_post_treatment) %>%
#     mutate(v_Loc = Loc_score[row_number()+1]- Loc_score[row_number()]) %>% #Velocity
#     
#     mutate(v_z_Loc = z_score_Loc[row_number()+1]- z_score_Loc[row_number()]) %>% #Velocity
#     mutate(a_z_Loc = v_z_Loc[row_number()+1]- v_z_Loc[row_number()]) %>% #acceleration
#     
#     mutate(v_Abund = as.numeric(LogAbundance[row_number()+1]- LogAbundance[row_number()])) %>%
#     mutate(v_z_Abund = as.numeric(z_score_logAbund[row_number()+1]- z_score_logAbund[row_number()])) %>%
#     mutate(a_z_Abund = as.numeric(v_z_Abund[row_number()+1]- v_z_Abund[row_number()])) %>%
#     
#     mutate(hyp = as.numeric(sqrt((vx_axis_m)^2 + (vy_axis_m)^2))) %>% #hyp for scaling
#     mutate(v_y_sc_hyp = as.numeric(vy_axis_m/hyp)) %>% #normalize by hyp to get vectors of the same size
#     mutate(v_x_sc_hyp = as.numeric(vx_axis_m/hyp)) %>%
#     ungroup() %>%
#     drop_na()


#----  
load_dataset <- function(protein = NULL){
  df_for_facet_abundance <<-  read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>%
    select(LogAbundance, Loc_score, z_score_Loc, z_score_logAbund, Protein, Frames_post_treatment, Cell_Barcode, Unique_pos, Relocalized, Yet) %>%
    # filter(Frames_post_treatment %in% c(0,8,16,24,32)) %>%
    filter(Protein == protein) %>%
    filter(Frames_post_treatment >= -8, Frames_post_treatment <= 32) %>%
    #* This could be changed to have more than just the post
    collect()
  
  return(df_for_facet_abundance)
}

get_xy_vars <- function(norm){
  if (norm){
    x_axis_s <- 'z_score_Loc'
    y_axis_s <- 'z_score_logAbund'
    vx_axis_s <- 'v_z_Loc'
    vy_axis_s <- 'v_z_Abund'
  } else{
    x_axis_s <- 'Loc_score'
    y_axis_s <- 'LogAbundance'
    vx_axis_s <- 'v_Loc'
    vy_axis_s <- 'v_Abund'
  }
  
  x_axis <- rlang::ensym(x_axis_s)
  y_axis <- rlang::ensym(y_axis_s)
  vx_axis <- rlang::ensym(vx_axis_s)
  vy_axis <- rlang::ensym(vy_axis_s)
  
  return(x_axis_s, y_axis_s, vx_axis_s, vy_axis_s, x_axis, y_axis, vx_axis, vy_axis)
}

create_segment_spaces <- function(df_input,n_breaks){
  if (n_breaks < 1){
    stop('n_breaks must be >= 1')
  }
  if (n_breaks > 4){
    stop('n_breaks is >4. This is going to be too dense because this will result in >25 bounded areas with descrete vectors.')
  }
  
  if (n_breaks == 1){
    break_1_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*1))
    break_2_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*2))
    
    break_1_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*1))
    break_2_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*2))
    
    df_test_nest <- df_input %>% 
      mutuate(quick_group = case_when(({{x_axis}} <= break_1_x & {{y_axis}} <= break_1_y) ~ "x0y0",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <= break_2_x & {{y_axis}} <= break_1_y) ~ "x1y0",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x0y1",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x1y1",
                                      TRUE ~ "Other")
      )
    
  }else if (n_breaks == 2){
    break_1_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*1))
    break_2_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*2))
    break_3_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*3))
    
    break_1_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*1))
    break_2_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*2))
    break_3_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*3))
    
    df_test_nest <- df_input %>% 
      mutuate(quick_group = case_when(({{x_axis}} <= break_1_x & {{y_axis}} <= break_1_y) ~ "x0y0",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <= break_2_x & {{y_axis}} <= break_1_y) ~ "x1y0",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} <= break_1_y) ~ "x2y0",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x0y1",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x1y1",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x2y1",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x0y2",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x1y2",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x2y2",
                                      TRUE ~ "Other")
      )
  }else if (n_breaks == 3){
    break_1_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*1))
    break_2_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*2))
    break_3_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*3))
    break_4_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*4))
    
    break_1_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*1))
    break_2_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*2))
    break_3_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*3))
    break_4_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*4))
    
    df_test_nest <- df_input %>% 
      mutuate(quick_group = case_when(({{x_axis}} <= break_1_x & {{y_axis}} <= break_1_y) ~ "x0y0",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <= break_2_x & {{y_axis}} <= break_1_y) ~ "x1y0",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} <= break_1_y) ~ "x2y0",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x & {{y_axis}} <= break_1_y) ~ "x3y0",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x0y1",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x1y1",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x2y1",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x  & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x3y1",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x0y2",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x1y2",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x2y2",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x3y2",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_3_y & {{y_axis}} <=  break_4_y) ~ "x0y3",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_3_y & {{y_axis}} <= break_4_y) ~ "x1y3",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_3_y & {{y_axis}} <= break_4_y) ~ "x2y3",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x & {{y_axis}} > break_3_y & {{y_axis}} <= break_4_y) ~ "x3y3",
                                      TRUE ~ "Other")
      )
  }else if (n_breaks == 4){
    break_1_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*1))
    break_2_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*2))
    break_3_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*3))
    break_4_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*4))
    break_5_x <- quantile(df_input[[x_axis_s]], (1/(n_breaks + 1)*5))
    
    break_1_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*1))
    break_2_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*2))
    break_3_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*3))
    break_4_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*4))
    break_5_y <- quantile(df_input[[y_axis_s]], (1/(n_breaks + 1)*5))
    
    df_test_nest <- df_input %>% 
      mutuate(quick_group = case_when(({{x_axis}} <= break_1_x & {{y_axis}} <= break_1_y) ~ "x0y0",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <= break_2_x & {{y_axis}} <= break_1_y) ~ "x1y0",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} <= break_1_y) ~ "x2y0",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x & {{y_axis}} <= break_1_y) ~ "x3y0",
                                      ({{x_axis}} > break_4_x & {{x_axis}} <= break_5_x & {{y_axis}} <= break_1_y) ~ "x4y0",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x0y1",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x1y1",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x2y1",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x  & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x3y1",
                                      ({{x_axis}} > break_4_x & {{x_axis}} <= break_5_x & {{y_axis}} > break_1_y & {{y_axis}} <= break_2_y) ~ "x4y1",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x0y2",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x1y2",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x2y2",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x3y2",
                                      ({{x_axis}} > break_4_x & {{x_axis}} <= break_5_x & {{y_axis}} > break_2_y & {{y_axis}} <= break_3_y) ~ "x4y2",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_3_y & {{y_axis}} <=  break_4_y) ~ "x0y3",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_3_y & {{y_axis}} <= break_4_y) ~ "x1y3",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_3_y & {{y_axis}} <= break_4_y) ~ "x2y3",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x & {{y_axis}} > break_3_y & {{y_axis}} <= break_4_y) ~ "x3y3",
                                      ({{x_axis}} > break_4_x & {{x_axis}} <= break_5_x & {{y_axis}} > break_3_y & {{y_axis}} <= break_4_y) ~ "x4y3",
                                      
                                      ({{x_axis}} <= break_1_x & {{y_axis}} > break_4_y & {{y_axis}} <= break_5_y) ~ "x0y4",
                                      ({{x_axis}} > break_1_x & {{x_axis}} <=break_2_x & {{y_axis}} > break_4_y & {{y_axis}} <= break_5_y) ~ "x1y4",
                                      ({{x_axis}} > break_2_x & {{x_axis}} <= break_3_x & {{y_axis}} > break_4_y & {{y_axis}} <= break_5_y) ~ "x2y4",
                                      ({{x_axis}} > break_3_x & {{x_axis}} <= break_4_x & {{y_axis}} > break_4_y & {{y_axis}} <= break_5_y) ~ "x3y4",
                                      ({{x_axis}} > break_4_x & {{x_axis}} <= break_5_x & {{y_axis}} > break_4_y & {{y_axis}} <= break_5_y) ~ "x4y4",
                                      TRUE ~ "Other")) %>%
      nest() %>% 
      mutate(median_x = map_dbl(data, ~mean(.x[[x_axis_s]]))) %>% 
      mutate(median_y = map_dbl(data, ~mean(.x[[y_axis_s]]))) %>% 
      mutate(median_v_y = map_dbl(data, ~mean(.x[[vy_axis_s]]))) %>% 
      mutate(median_v_x = map_dbl(data, ~mean(.x[[vx_axis_s]]))) %>% 
      mutate(count = map_dbl(data, ~length(.x[[vy_axis_s]]))) %>% 
      filter(count > 5) %>% 
      mutate(count_s = log(count)/3)
  }
  
  return(df_test_nest)
}

get_velocity <- function(df_input, x_axis_s, y_axis_s, vx_axis_s, vy_axis_s, x_axis, y_axis, vx_axis, vy_axis){
  df <- df_input %>%
    filter(Protein == protein) %>%
    unique() %>%
    group_by(Cell_Barcode) %>%
    # mutate(Yet = as_factor(Yet)) %>% 
    # mutate(Relocalized = as_factor(Relocalized)) %>%
    arrange(Frames_post_treatment) %>%
    filter(Frames_post_treatment < (frame + space), Frames_post_treatment >= frame)%>%
    
    mutate(v_Loc = Loc_score[row_number()+1]- Loc_score[row_number()]) %>% #Velocity
    mutate(v_z_Loc = z_score_Loc[row_number()+1]- z_score_Loc[row_number()]) %>% #Velocity
    mutate(a_z_Loc = v_z_Loc[row_number()+1]- v_z_Loc[row_number()]) %>% #acceleration
    
    mutate(v_Abund = as.numeric(LogAbundance[row_number()+1]- LogAbundance[row_number()])) %>%
    mutate(v_z_Abund = as.numeric(z_score_logAbund[row_number()+1]- z_score_logAbund[row_number()])) %>%
    mutate(a_z_Abund = as.numeric(v_z_Abund[row_number()+1]- v_z_Abund[row_number()])) %>%
    
    mutate(vx_axis_m = mean({{vx_axis}}, na.rm = TRUE)) %>% 
    mutate(vy_axis_m = mean({{vy_axis}}, na.rm = TRUE)) %>% 
    
    mutate(hyp = as.numeric(sqrt((vx_axis_m)^2 + (vy_axis_m)^2))) %>% #hyp for scaling
    mutate(v_y_sc_hyp = as.numeric(vy_axis_m/hyp)) %>% #normalize by hyp to get vectors of the same size
    mutate(v_x_sc_hyp = as.numeric(vx_axis_m/hyp)) %>%
    ungroup() %>%
    drop_na()
  return(df)
}

complex_velocity <- function(protein, frame, space = 8, norm = TRUE, min_x = NULL, max_x = NULL, min_y = NULL, max_y = NULL){
  mypal <- pal_npg("nrc", alpha = 0.7)(3)
  
  
  if (exists('df_for_facet_abundance')){
    if (0 < df_for_facet_abundance %>%
        filter(Protein == protein) %>%
        nrow()){} else{
          load_dataset(protein = protein)
        }
  } else{
    load_dataset(protein = protein)
  }
  
  
  if (norm){
    x_axis_s <- 'z_score_Loc'
    y_axis_s <- 'z_score_logAbund'
    vx_axis_s <- 'v_z_Loc'
    vy_axis_s <- 'v_z_Abund'
  } else{
    x_axis_s <- 'Loc_score'
    y_axis_s <- 'LogAbundance'
    vx_axis_s <- 'v_Loc'
    vy_axis_s <- 'v_Abund'
  }
  
  x_axis <- rlang::ensym(x_axis_s)
  y_axis <- rlang::ensym(y_axis_s)
  vx_axis <- rlang::ensym(vx_axis_s)
  vy_axis <- rlang::ensym(vy_axis_s)
  
  ##### get max_mins ####
  if (is.null(min_x)){
  min_x <- df_for_facet_abundance %>%
    select({{x_axis}}) %>%
    min() %>%
    round(2)
  min_x <- (min_x - 0.1) %>%
    round(1)
  }
  
  if (is.null(max_x)){
  max_x <- df_for_facet_abundance %>%
    select({{x_axis}}) %>%
    max() %>%
    round(2)
  max_x <- max_x + 0.1 %>%
    round(1)
  }

  if (is.null(min_y)){
  min_y <- df_for_facet_abundance %>%
    select({{y_axis}}) %>%
    min() %>%
    round(2)
  min_y <- (min_y - 0.1) %>%
    round(1)
  }
  
  if (is.null(max_y)){
  max_y <- df_for_facet_abundance %>%
    select({{y_axis}}) %>%
    max() %>%
    round(2)
  max_y <- (max_y + 0.1) %>%
    round(1)
  }

  # df_for_facet_abundance %>%
  #   filter(Relocalized == 1) %>%
  #   select(Loc_score) %>%
  #   mutate(Loc_score = as.numeric(Loc_score)) %>%
  #   summarise(iqr = IQR(Loc_score)) * 3
  
  # max_x <- df_for_facet_abundance %>%
  #   filter(Relocalized == 1) %>%
  #   mutate(Loc_score = as.numeric(Loc_score)) %>%
  #   select(Loc_score) %>%
  #   summarise(qtile_m = quantile(Loc_score, 0.99)) %>%
  #   select(qtile_m) %>%
  #   round(1) %>%
  #   as.numeric()
  
  #####
  
  df <- df_for_facet_abundance %>%
    filter(Protein == protein) %>%
    unique() %>%
    group_by(Cell_Barcode) %>%
    # mutate(Yet = as_factor(Yet)) %>% 
    # mutate(Relocalized = as_factor(Relocalized)) %>%
    arrange(Frames_post_treatment) %>%
    filter(Frames_post_treatment < (frame + space), Frames_post_treatment >= frame)%>%
    
    mutate(v_Loc = Loc_score[row_number()+1]- Loc_score[row_number()]) %>% #Velocity
    mutate(v_z_Loc = z_score_Loc[row_number()+1]- z_score_Loc[row_number()]) %>% #Velocity
    mutate(a_z_Loc = v_z_Loc[row_number()+1]- v_z_Loc[row_number()]) %>% #acceleration
    
    mutate(v_Abund = as.numeric(LogAbundance[row_number()+1]- LogAbundance[row_number()])) %>%
    mutate(v_z_Abund = as.numeric(z_score_logAbund[row_number()+1]- z_score_logAbund[row_number()])) %>%
    mutate(a_z_Abund = as.numeric(v_z_Abund[row_number()+1]- v_z_Abund[row_number()])) %>%
    
    mutate(vx_axis_m = mean({{vx_axis}}, na.rm = TRUE)) %>% 
    mutate(vy_axis_m = mean({{vy_axis}}, na.rm = TRUE)) %>% 
    
    mutate(hyp = as.numeric(sqrt((vx_axis_m)^2 + (vy_axis_m)^2))) %>% #hyp for scaling
    mutate(v_y_sc_hyp = as.numeric(vy_axis_m/hyp)) %>% #normalize by hyp to get vectors of the same size
    mutate(v_x_sc_hyp = as.numeric(vx_axis_m/hyp)) %>%
    ungroup() %>%
    drop_na()
  
  
  df_temp_df <- df %>%
    filter(#{{y_axis}} <= 4, {{x_axis}} <= 4,
           Frames_post_treatment == frame) %>%
    select({{x_axis}}, {{y_axis}}, {{vx_axis}}, {{vy_axis}}, vx_axis_m, vy_axis_m, hyp, v_y_sc_hyp, v_x_sc_hyp, Relocalized, Yet)
  
  df_3_point <- df_temp_df %>%
    group_by(Relocalized) %>%
    summarise(lowest_Abund = min({{y_axis}}, na.rm = T), highest_Abund = max({{y_axis}}, na.rm = T),
              lowest_LocScore = min({{x_axis}}, na.rm = T), highest_LocScore = max({{x_axis}}, na.rm = T)) %>% 
    ungroup()
  
  mean_Loc <- mean(df_temp_df[[x_axis_s]])
  mean_Abund <- mean(df_temp_df[[y_axis_s]])
  mean_v_x <- mean(df_temp_df[[vx_axis_s]])*10
  mean_v_y <- mean(df_temp_df[[vy_axis_s]])*10
  quant_mid_hyp <- quantile(df_temp_df$hyp, 0.8)
  quant_20th_x <- quantile(df_temp_df[[x_axis_s]], 0.2)
  quant_40th_x <- quantile(df_temp_df[[x_axis_s]], 0.4)
  quant_60th_x <- quantile(df_temp_df[[x_axis_s]], 0.6)
  quant_80th_x <- quantile(df_temp_df[[x_axis_s]], 0.8)
  quant_100th_x <- quantile(df_temp_df[[x_axis_s]], 1)
  
  quant_20th_y <- quantile(df_temp_df[[y_axis_s]], 0.2)
  quant_40th_y <- quantile(df_temp_df[[y_axis_s]], 0.4)
  quant_60th_y <- quantile(df_temp_df[[y_axis_s]], 0.6)
  quant_80th_y <- quantile(df_temp_df[[y_axis_s]], 0.8)
  quant_100th_y <- quantile(df_temp_df[[y_axis_s]], 1)
  
  df_test_nest <- df_temp_df %>%
    mutate(quick_group = case_when(({{x_axis}} <= quant_20th_x & {{y_axis}} <= quant_20th_y) ~ "x0y0",
                                   ({{x_axis}} > quant_20th_x & {{x_axis}} <= quant_40th_x & {{y_axis}} <= quant_20th_y) ~ "x1y0",
                                   ({{x_axis}} > quant_40th_x & {{x_axis}} <= quant_60th_x & {{y_axis}} <= quant_20th_y) ~ "x2y0",
                                   ({{x_axis}} > quant_60th_x & {{x_axis}} <= quant_80th_x & {{y_axis}} <= quant_20th_y) ~ "x3y0",
                                   ({{x_axis}} > quant_80th_x & {{x_axis}} <= quant_100th_x & {{y_axis}} <= quant_20th_y) ~ "x4y0",
                                   
                                   ({{x_axis}} <= quant_20th_x & {{y_axis}} > quant_20th_y & {{y_axis}} <= quant_40th_y) ~ "x0y1",
                                   ({{x_axis}} > quant_20th_x & {{x_axis}} <=quant_40th_x & {{y_axis}} > quant_20th_y & {{y_axis}} <= quant_40th_y) ~ "x1y1",
                                   ({{x_axis}} > quant_40th_x & {{x_axis}} <= quant_60th_x & {{y_axis}} > quant_20th_y & {{y_axis}} <= quant_40th_y) ~ "x2y1",
                                   ({{x_axis}} > quant_60th_x & {{x_axis}} <= quant_80th_x  & {{y_axis}} > quant_20th_y & {{y_axis}} <= quant_40th_y) ~ "x3y1",
                                   ({{x_axis}} > quant_80th_x & {{x_axis}} <= quant_100th_x & {{y_axis}} > quant_20th_y & {{y_axis}} <= quant_40th_y) ~ "x4y1",
                                   
                                   ({{x_axis}} <= quant_20th_x & {{y_axis}} > quant_40th_y & {{y_axis}} <= quant_60th_y) ~ "x0y2",
                                   ({{x_axis}} > quant_20th_x & {{x_axis}} <=quant_40th_x & {{y_axis}} > quant_40th_y & {{y_axis}} <= quant_60th_y) ~ "x1y2",
                                   ({{x_axis}} > quant_40th_x & {{x_axis}} <= quant_60th_x & {{y_axis}} > quant_40th_y & {{y_axis}} <= quant_60th_y) ~ "x2y2",
                                   ({{x_axis}} > quant_60th_x & {{x_axis}} <= quant_80th_x & {{y_axis}} > quant_40th_y & {{y_axis}} <= quant_60th_y) ~ "x3y2",
                                   ({{x_axis}} > quant_80th_x & {{x_axis}} <= quant_100th_x & {{y_axis}} > quant_40th_y & {{y_axis}} <= quant_60th_y) ~ "x4y2",
                                   
                                   ({{x_axis}} <= quant_20th_x & {{y_axis}} > quant_60th_y & {{y_axis}} <=  quant_80th_y) ~ "x0y3",
                                   ({{x_axis}} > quant_20th_x & {{x_axis}} <=quant_40th_x & {{y_axis}} > quant_60th_y & {{y_axis}} <= quant_80th_y) ~ "x1y3",
                                   ({{x_axis}} > quant_40th_x & {{x_axis}} <= quant_60th_x & {{y_axis}} > quant_60th_y & {{y_axis}} <= quant_80th_y) ~ "x2y3",
                                   ({{x_axis}} > quant_60th_x & {{x_axis}} <= quant_80th_x & {{y_axis}} > quant_60th_y & {{y_axis}} <= quant_80th_y) ~ "x3y3",
                                   ({{x_axis}} > quant_80th_x & {{x_axis}} <= quant_100th_x & {{y_axis}} > quant_60th_y & {{y_axis}} <= quant_80th_y) ~ "x4y3",
                                   
                                   ({{x_axis}} <= quant_20th_x & {{y_axis}} > quant_80th_y & {{y_axis}} <= quant_100th_y) ~ "x0y4",
                                   ({{x_axis}} > quant_20th_x & {{x_axis}} <=quant_40th_x & {{y_axis}} > quant_80th_y & {{y_axis}} <= quant_100th_y) ~ "x1y4",
                                   ({{x_axis}} > quant_40th_x & {{x_axis}} <= quant_60th_x & {{y_axis}} > quant_80th_y & {{y_axis}} <= quant_100th_y) ~ "x2y4",
                                   ({{x_axis}} > quant_60th_x & {{x_axis}} <= quant_80th_x & {{y_axis}} > quant_80th_y & {{y_axis}} <= quant_100th_y) ~ "x3y4",
                                   ({{x_axis}} > quant_80th_x & {{x_axis}} <= quant_100th_x & {{y_axis}} > quant_80th_y & {{y_axis}} <= quant_100th_y) ~ "x4y4",
                                   TRUE ~ "Other")) %>%
    filter(quick_group != "Other") %>% 
    group_by(quick_group) %>% 
    nest() %>% 
    mutate(median_x = map_dbl(data, ~mean(.x[[x_axis_s]]))) %>% 
    mutate(median_y = map_dbl(data, ~mean(.x[[y_axis_s]]))) %>% 
    mutate(median_v_y = map_dbl(data, ~mean(.x[[vy_axis_s]]))) %>% 
    mutate(median_v_x = map_dbl(data, ~mean(.x[[vx_axis_s]]))) %>% 
    mutate(count = map_dbl(data, ~length(.x[[vy_axis_s]]))) %>% 
    filter(count > 5) %>% 
    mutate(count_s = log(count)/3)

  
    point1 <-  df_3_point %>%
      filter(Relocalized == 0) %>% 
      select(highest_Abund, lowest_LocScore) %>% 
      mutate(x_p = lowest_LocScore, y_p = highest_Abund)
    point2 <-  df_3_point %>%
      filter(Relocalized == 1) %>% 
      select(lowest_Abund, lowest_LocScore) %>% 
      mutate(x_p = lowest_LocScore, y_p = lowest_Abund)
    point3 <-  df_3_point %>%
      filter(Relocalized == 1) %>% 
      select(lowest_Abund, highest_LocScore) %>% 
      mutate(x_p = highest_LocScore, y_p = lowest_Abund)
    
    
  df_path <- df_for_facet_abundance %>% 
    filter(Frames_post_treatment >= 0, Frames_post_treatment > (frame - 8), Frames_post_treatment <= frame) %>% 
    arrange(Cell_Barcode, Frames_post_treatment) %>% 
    select({{x_axis}}, {{y_axis}}, Frames_post_treatment, Cell_Barcode)
  
p <- ggplot(data = df_temp_df, aes(x = {{x_axis}}, y = {{y_axis}}))+
    geom_point(mapping = aes({{x_axis}},{{y_axis}}), alpha = 0.7, color = '#8491B4B2')+
    # geom_line(data = filter(Frames_post_treatment >= 0, Frames_post_treatment > (frame - 8), Frames_post_treatment <= frame),
    #           mapping = aes(x = {{x_axis}}, y = {{y_axis}}, group = (Cell_Barcode)))
    # geom_line(data = df_path,
    #             mapping = aes(x = {{x_axis}}, y = {{y_axis}}, colour = Frames_post_treatment, group = Cell_Barcode), method = 'loess')+ #, method = 'lm', se = F)+
    # geom_smooth(data = df_path,
    #             mapping = aes(x = {{x_axis}}, y = {{y_axis}}, group = Frames_post_treatment), method = 'loess', se = F)+ #, method = 'lm', se = F)+
    # geom_path(data = df_path,
    #           mapping = aes(x = {{x_axis}}, y = {{y_axis}}, color = Cell_Barcode)) +
    # geom_segment(data = (df_temp_df %>% filter(hyp <=quant_mid_hyp)),
    #              aes(x = {{x_axis}}, y = {{y_axis}},
    #                  xend = ({{x_axis}} + v_x_sc_hyp/15),
    #                  yend = ({{y_axis}}+ v_y_sc_hyp/15),
    #                  colour = hyp,
    #                  alpha = 0.2),
    #              arrow = arrow(length = unit(0.5, "cm")), linewidth = 1)+
    geom_segment(data = (df_temp_df %>% filter(hyp > quant_mid_hyp)),
                 aes(x = {{x_axis}}, y = {{y_axis}},
                     xend = ({{x_axis}} + v_x_sc_hyp/15),
                     yend = ({{y_axis}}+ v_y_sc_hyp/15),
                     colour = hyp,
                     alpha = 0.7),
                 arrow = arrow(length = unit(0.5, "cm")), linewidth = 1)+
    scale_color_gradientn(colours = mypal)+
    geom_segment(data = df_test_nest,
                 aes(x = median_x,
                     y = median_y,
                     xend = (median_x + median_v_x/5),
                     yend = (median_y + median_v_y/5),
                     alpha = count_s),
                 colour = '#F39B7FB2',
                 arrow = arrow(length = unit(0.6, "cm")),
                 linewidth = 2)+
    
    # geom_segment(data = (df_temp_df %>% filter({{y_axis}} <= 2.5, {{x_axis}} <= 2.5)),
    #              aes(x = (mean({{x_axis}})),
    #                  y = (mean({{x_}})),
    #                  xend = (mean({{x_axis}}) + (mean(v_z_Loc)/5)),
    #                  yend= (mean({{x_axis}}) + (mean(v_z_Abund)/5)),
    #                  color = 'black'))+
    # 
    # annotate("segment",
    #          x = (mean_Loc),
    #          y = (mean_Abund),
    #          xend = (mean_Loc + mean_v_x),
    #          yend= (mean_Abund + mean_v_y),
    #          color = 'black',
    #          linewidth = 2,
    #          arrow = arrow(length = unit(0.6, "cm")))+
  
    annotate("segment",
             x = (point1$x_p),
             y = (point1$y_p),
             xend = (point2$x_p),
             yend= (point2$y_p),
             color = '#3C5488B2',
             linewidth = 2)+
    annotate("segment",
             x = (point2$x_p),
             y = (point2$y_p),
             xend = (point3$x_p),
             yend= (point3$y_p),
             color = '#3C5488B2',
             linewidth = 2)+
  geom_vline(xintercept = 1, alpha = 0.8, linetype = 'dotdash')+
  # ggside::ggside_layer(geom = "beeswarm",
  #                      data = df_temp_df,
  #                      mapping = aes(x = Yet, y = Loc_score, color = Relocalized, group = Yet),
  #              side = "x")
  geom_xsideboxplot(data = df_temp_df, mapping = aes(y = Yet, group = Yet, colour = Yet), orientation = 'y')+
  geom_ysideboxplot(data = df_temp_df, mapping = aes(x = Yet, group = Yet, colour = Yet), orientation = 'x')+
  scale_xsidey_discrete()+
  # # geom_ysidedensity(data = df_temp_df, mapping = aes(x = after_stat(density), group = Yet, yfill = Yet), position = "stack",) +
  scale_ysidex_discrete()+
  scale_xcolor_discrete()+
  scale_ycolor_discrete()+
  scale_xfill_discrete()+
  scale_yfill_discrete()+
  scale_x_continuous(limits = c(min_x, max_x))+
  scale_y_continuous(limits = c(min_y, max_y))+
  # # scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL)+ 
  # theme_ggside_minimal()+
  theme_minimal()+
  theme(legend.position = "none")
  
  # geom_ysideboxplot(data = df_temp_df, mapping = aes(x = ))
  # geom_xsidedensity(data = df_temp_df,
  #                     aes(x = {{x_axis}}, y = after_stat(density(x = {{y_axis}}))))
  return(p)
}

complex_heat <- function(protein, frame, space = 8, norm = TRUE, change_deg = 0, min_x = NULL, max_x = NULL, min_y = NULL, max_y = NULL){
  mypal <- pal_npg("nrc", alpha = 0.7)(9)
  
  if (exists('df_for_facet_abundance')){
    if (0 < df_for_facet_abundance %>%
        filter(Protein == protein) %>%
        nrow()){} else{
          load_dataset(protein = protein)
        }
  } else{
    load_dataset(protein = protein)
  }
  
  axis_variables <- get_xy_vars(norm = norm)
  
  #! The below could be replaced by modifying all calls to use "axis_variables[n]" notation
  x_axis_s = axis_variables[1]
  y_axis_s = axis_variables[2]
  vx_axis_s = axis_variables[3]
  vy_axis_s = axis_variables[4]
  x_axis = axis_variables[5]
  y_axis = axis_variables[6]
  vx_axis = axis_variables[7]
  vy_axis = axis_variables[8]
  
  ##### get max_mins ####
  if (is.null(min_x)){
    min_x <- df_for_facet_abundance %>%
      select({{x_axis}}) %>%
      min() %>%
      round(2)
    min_x <- (min_x - 0.1) %>%
      round(1)
  }
  
  if (is.null(max_x)){
    max_x <- df_for_facet_abundance %>%
      select({{x_axis}}) %>%
      max() %>%
      round(2)
    max_x <- max_x + 0.1 %>%
      round(1)
  }
  
  if (is.null(min_y)){
    min_y <- df_for_facet_abundance %>%
      select({{y_axis}}) %>%
      min() %>%
      round(2)
    min_y <- (min_y - 0.1) %>%
      round(1)
  }
  
  if (is.null(max_y)){
    max_y <- df_for_facet_abundance %>%
      select({{y_axis}}) %>%
      max() %>%
      round(2)
    max_y <- (max_y + 0.1) %>%
      round(1)
  }
  
  # df_for_facet_abundance %>%
  #   filter(Relocalized == 1) %>%
  #   select(Loc_score) %>%
  #   mutate(Loc_score = as.numeric(Loc_score)) %>%
  #   summarise(iqr = IQR(Loc_score)) * 3
  
  # max_x <- df_for_facet_abundance %>%
  #   filter(Relocalized == 1) %>%
  #   mutate(Loc_score = as.numeric(Loc_score)) %>%
  #   select(Loc_score) %>%
  #   summarise(qtile_m = quantile(Loc_score, 0.99)) %>%
  #   select(qtile_m) %>%
  #   round(1) %>%
  #   as.numeric()
  
  #####
  
  df <- df_for_facet_abundance %>% 
    get_velocity(.,x_axis_s = x_axis_s, y_axis_s = y_axis_s, vx_axis_s = vx_axis_s, vy_axis_s = vy_axis_s, x_axis = x_axis, y_axis = y_axis, vx_axis = vx_axis, vy_axis = vy_axis)
  
  df_temp_df <- df %>%
    filter(#{{y_axis}} <= 4, {{x_axis}} <= 4,
      Frames_post_treatment == frame) %>%
    select({{x_axis}}, {{y_axis}}, {{vx_axis}}, {{vy_axis}}, vx_axis_m, vy_axis_m, hyp, v_y_sc_hyp, v_x_sc_hyp, Relocalized, Yet)
  
  mean_Loc <- mean(df_temp_df[[x_axis_s]])
  mean_Abund <- mean(df_temp_df[[y_axis_s]])
  mean_v_x <- mean(df_temp_df[[vx_axis_s]])*10
  mean_v_y <- mean(df_temp_df[[vy_axis_s]])*10
  quant_mid_hyp <- quantile(df_temp_df$hyp, 0.8)
  
  df_test_nest <- create_segment_spaces(df_temp_df)
  
  
  if (change_deg == 0){
    
  } else if (change_deg == 0){
    prelude <- 'Velocity of'
  }else if (change_deg == 0){
    prelude <- 'Acceleration of'
  }
  
  p <- ggplot(data = df_temp_df, aes(x = {{x_axis}}, y = {{y_axis}}))+
    heat
    
    
    
    # geom_line(data = filter(Frames_post_treatment >= 0, Frames_post_treatment > (frame - 8), Frames_post_treatment <= frame),
    #           mapping = aes(x = {{x_axis}}, y = {{y_axis}}, group = (Cell_Barcode)))
    # geom_line(data = df_path,
    #             mapping = aes(x = {{x_axis}}, y = {{y_axis}}, colour = Frames_post_treatment, group = Cell_Barcode), method = 'loess')+ #, method = 'lm', se = F)+
    # geom_smooth(data = df_path,
    #             mapping = aes(x = {{x_axis}}, y = {{y_axis}}, group = Frames_post_treatment), method = 'loess', se = F)+ #, method = 'lm', se = F)+
    # geom_path(data = df_path,
    #           mapping = aes(x = {{x_axis}}, y = {{y_axis}}, color = Cell_Barcode)) +
    # geom_segment(data = (df_temp_df %>% filter(hyp <=quant_mid_hyp)),
    #              aes(x = {{x_axis}}, y = {{y_axis}},
    #                  xend = ({{x_axis}} + v_x_sc_hyp/15),
    #                  yend = ({{y_axis}}+ v_y_sc_hyp/15),
    #                  colour = hyp,
    #                  alpha = 0.2),
    #              arrow = arrow(length = unit(0.5, "cm")), linewidth = 1)+
    geom_segment(data = (df_temp_df %>% filter(hyp > quant_mid_hyp)),
                 aes(x = {{x_axis}}, y = {{y_axis}},
                     xend = ({{x_axis}} + v_x_sc_hyp/15),
                     yend = ({{y_axis}}+ v_y_sc_hyp/15),
                     colour = hyp,
                     alpha = 0.7),
                 arrow = arrow(length = unit(0.5, "cm")), linewidth = 1)+
    scale_color_gradientn(colours = mypal)+
    geom_segment(data = df_test_nest,
                 aes(x = median_x,
                     y = median_y,
                     xend = (median_x + median_v_x/5),
                     yend = (median_y + median_v_y/5),
                     alpha = count_s),
                 colour = '#F39B7FB2',
                 arrow = arrow(length = unit(0.6, "cm")),
                 linewidth = 2)+
    
    # geom_segment(data = (df_temp_df %>% filter({{y_axis}} <= 2.5, {{x_axis}} <= 2.5)),
    #              aes(x = (mean({{x_axis}})),
    #                  y = (mean({{x_}})),
    #                  xend = (mean({{x_axis}}) + (mean(v_z_Loc)/5)),
    #                  yend= (mean({{x_axis}}) + (mean(v_z_Abund)/5)),
    #                  color = 'black'))+
    # 
    # annotate("segment",
    #          x = (mean_Loc),
    #          y = (mean_Abund),
    #          xend = (mean_Loc + mean_v_x),
    #          yend= (mean_Abund + mean_v_y),
    #          color = 'black',
    #          linewidth = 2,
    #          arrow = arrow(length = unit(0.6, "cm")))+
    
    annotate("segment",
             x = (point1$x_p),
             y = (point1$y_p),
             xend = (point2$x_p),
             yend= (point2$y_p),
             color = '#3C5488B2',
             linewidth = 2)+
    annotate("segment",
             x = (point2$x_p),
             y = (point2$y_p),
             xend = (point3$x_p),
             yend= (point3$y_p),
             color = '#3C5488B2',
             linewidth = 2)+
    geom_vline(xintercept = 1, alpha = 0.8, linetype = 'dotdash')+
    # ggside::ggside_layer(geom = "beeswarm",
    #                      data = df_temp_df,
    #                      mapping = aes(x = Yet, y = Loc_score, color = Relocalized, group = Yet),
    #              side = "x")
    geom_xsideboxplot(data = df_temp_df, mapping = aes(y = Yet, group = Yet, colour = Yet), orientation = 'y')+
    geom_ysideboxplot(data = df_temp_df, mapping = aes(x = Yet, group = Yet, colour = Yet), orientation = 'x')+
    scale_xsidey_discrete()+
    # # geom_ysidedensity(data = df_temp_df, mapping = aes(x = after_stat(density), group = Yet, yfill = Yet), position = "stack",) +
    scale_ysidex_discrete()+
    scale_xcolor_discrete()+
    scale_ycolor_discrete()+
    scale_xfill_discrete()+
    scale_yfill_discrete()+
    scale_x_continuous(limits = c(min_x, max_x))+
    scale_y_continuous(limits = c(min_y, max_y))+
    # # scale_ysidex_continuous(guide = guide_axis(angle = 90), minor_breaks = NULL)+ 
    # theme_ggside_minimal()+
    theme_minimal()+
    theme(legend.position = "none")
  
  # geom_ysideboxplot(data = df_temp_df, mapping = aes(x = ))
  # geom_xsidedensity(data = df_temp_df,
  #                     aes(x = {{x_axis}}, y = after_stat(density(x = {{y_axis}}))))
  return(p)
}

# complex_velocity(df_for_facet_abundance, frame = 0, protein = 'RAD51', space = 8, norm = FALSE)

  # theme(legend.position = "none")


#------- Actual plotting ------------

plot_grid(complex_velocity(frame = -8, protein = 'RAD51', space = 8,norm = FALSE),
          complex_velocity(frame = 0, protein = 'RAD51', space = 8,norm = FALSE),
          complex_velocity(frame = 8, protein = 'RAD51', space = 8,norm = FALSE),
          complex_velocity(frame = 16,protein = 'RAD51', space = 8, norm = FALSE),
          complex_velocity(frame = 24,protein = 'RAD51', space = 8, norm = FALSE),
          complex_velocity(frame = 30,protein = 'RAD51', space = 8, norm = FALSE),
          labels = "AUTO", align = "h")


plot_grid(complex_heat(frame = -8, protein = 'RAD51', space = 8,norm = FALSE),
          complex_heat(frame = 0, protein = 'RAD51', space = 8,norm = FALSE),
          complex_heat(frame = 8, protein = 'RAD51', space = 8,norm = FALSE),
          complex_heat(frame = 16,protein = 'RAD51', space = 8, norm = FALSE),
          complex_heat(frame = 24,protein = 'RAD51', space = 8, norm = FALSE),
          complex_heat(frame = 30,protein = 'RAD51', space = 8, norm = FALSE),
          labels = "AUTO", align = "h")

#------------misc--------




keep_list <- df$Cell_Barcode %>%
  unique() %>%
  sample(50)

anim <- ggplot(df %>%
                 filter(z_score_logAbund <= 2.5, z_score_Loc <= 2.5) %>%
                 filter(Cell_Barcode %in% keep_list)) +
  transition_time(Frames_post_treatment)+
  # transition_reveal(Frames_post_treatment)
  geom_point(aes(x = Loc_score, y = LogAbundance,
                 colour = Cell_Barcode,
                 group = Cell_Barcode,
                 size = 2))+
  shadow_wake(wake_length = 0.1,
              alpha = FALSE)+
  theme(legend.position = "none")
anim



  # labs(title = 'Frame: {frame_time}', x = 'Z_ReLocScor', y = 'Z_Abundance') +
  # transition_states(Frames_post_treatment,#,
  #                   transition_length = 2,
  #                   state_length = 1)+
  shadow_wake(wake_length = 0.1,
              alpha = FALSE)+
  # shadow_mark(alpha = 0.6)
  # ease_aes('linear')


anim
#
animate(anim, fps = 24)

animate(
  anim +
    enter_fade()+
    exit_fade(),
  render = av_renderer(),
  width = 1000, height = 1000, res = 300
)

anim_save("test.gif", animation = anim)
`# anim_save("271-ggplot2-animated-gif-chart-with-gganimate2.gif")
