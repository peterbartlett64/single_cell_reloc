library(ggplot2)
library(wesanderson)
library(ggsci)
library(ggrepel)
library(tidymodels)
library(ggrepel)
library(gt)
library(data.table)
library(forcats)
library(arrow)
library(dplyr)


date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

#Load in the dataframe
# df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>%
#   filter(Frames_post_treatment >= 0)

# df_k <- read.csv("D:/ALL_FINAL/Combined_by_perc/median_kendall.csv") %>% 
#   filter(!(Protein %in% c("PPH22d0210", "LCD1", "RAD5", "RFA1d0210", "LSM3d0210r1", "RRP17d0210", "SLX8", "XRS2d0210", "Contam.")) | (str_like(Protein,'DDC2')))
# df_p <- read.csv("D:/ALL_FINAL/Combined_by_perc/median_pearson.csv") %>% 
#   filter(!(Protein %in% c("PPH22d0210", "LCD1", "RAD5", "RFA1d0210", "LSM3d0210r1", "RRP17d0210", "SLX8", "XRS2d0210", "Contam.")) | (str_like(Protein,'DDC2')))

# This is not currently being used
check_normality_me <- function(protein){
  df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>%
    filter(Frames_post_treatment >= 0, Protein == protein) %>% 
    collect() %>% 
    group_by(Frames_post_treatment) %>% 
    nest() %>%
    mutate(shapiro_A = map(data, ~shapiro.test(.x$z_score_logAbund)))
    
  return(df)
}
# 
# df_n <- df %>% 
#   group_by(Frames_post_treatment) %>% 
#   nest() %>%
#   mutate(shapiro_A = map(data, ~shapiro.test(.x$z_score_logAbund)))
# library("ggpubr")
# 
# ggscatter(df, x = "z_score_Loc", y ="z_score_logAbund", colour = "Unique_pos", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson")+
#   facet_wrap(~Frames_post_treatment)


#* Getting the correlation values again is not needed since did it in Python

# protein_corr <- function(df){
#   model <- cor.test(df$z_score_logAbund, df$z_score_Loc, method="pearson")
#   return(model)
# }
# 
# ### Confirmation of
# df_fitlered <- df %>% 
#   filter(Protein == 'SAE2') %>% 
#   filter(Frames_post_treatment %in% c('0', '8', '16', '24', '32')) %>%  # This is the consistent time
#   collect() %>% 
#   nest_by(Frames_post_treatment) %>% 
#   mutate(model = map(data, protein_corr, .progress = T),
#          model_estimate = model$estimate,
#          model_p = model$p.value,
#          model_tidy = map(model, tidy))


make_hose_spray <- function(corr_type, ll){
  tkach <- read_xls("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/Tkach_refined.xls") %>% 
    rename(Protein = "Standard Name") %>% 
    rename(MMS_Localization_change = "MMS localization change class") %>% 
    select(Protein, MMS_Localization_change)
  
  df <- read.csv(sprintf("D:/ALL_FINAL/Combined_by_perc/median_%s.csv", corr_type)) %>% 
    filter(!(Protein %in% c("PPH22d0210", "LCD1", "RAD5", "RFA1d0210", "LSM3d0210r1", "RRP17d0210", "SLX8", "XRS2d0210", "Contam.")))# %>% 
    # full_join(tkach, by = c(Protein = "Protein")) %>% 
    # drop_na()
  
  pal_g <- wes_palette("AsteroidCity2")
  pal <- wes_palette("AsteroidCity3", 20, type = "continuous")
  plot <- df %>% 
    # filter(-log10(pval) <= 50, corr <0.4, corr > -0.4) %>%
    {ggplot(data = ., aes(x=corr, y=-log10(pval), col=corr, label=Protein))+
    geom_point() + 
    geom_vline(xintercept=c(-0.3), col=pal_g[1]) +    
    geom_vline(xintercept=c(-0.1), col=pal_g[2]) +
    geom_vline(xintercept=c(0.1), col=pal_g[2]) +
    geom_vline(xintercept=c(0.3), col=pal_g[3]) +
    geom_vline(xintercept=c(0.7), col=pal_g[4]) +
    geom_vline(xintercept=c(0.7), col=pal_g[5]) +
        
    # geom_vline(xintercept=c(0.8), col=pal_g[6]) +
    geom_hline(yintercept=-log10(0.05), col=pal_g[2])+
    geom_text_repel(data = . %>% filter((pval < 0.05) &(Protein %in% c('SAE2', 'FLR1', 'EXO1', 'ECO1', 'RFA1d0213', 'HAC1')| corr > ll | corr < -ll)), box.padding = 0.5, max.overlaps = 40)+
    scale_colour_gradientn(colours = pal)+
    scale_x_continuous(breaks = c(-1, -0.7, -0.3, -0.1, 0, 0.1, 0.3, 0.7, 1))+
    scale_y_continuous(trans='log10')+
    theme_minimal()+
    labs(title=sprintf("Correlation by %s", corr_type), x="Correlation", y="-log10(p-value)")}
  
  return(plot)
}

make_hose_spray <- function(corr_type, ll){
  df <- read.csv(sprintf("D:/ALL_FINAL/Combined_by_perc/median_%s.csv", corr_type)) %>% 
    filter(!(Protein %in% c("PPH22d0210", "LCD1", "RAD5", "RFA1d0210", "LSM3d0210r1", "RRP17d0210", "SLX8", "XRS2d0210", "Contam.")))
  
  pal_g <- wes_palette("AsteroidCity2")
  pal <- wes_palette("AsteroidCity3", 20, type = "continuous")
  plot <- df %>% 
    # filter(-log10(pval) <= 50, corr <0.4, corr > -0.4) %>%
    {ggplot(data = ., aes(x=corr, y=-log10(pval), col=corr, label=Protein))+
    geom_point() + 
    geom_vline(xintercept=c(-0.3), col=pal_g[1]) +    
    geom_vline(xintercept=c(-0.1), col=pal_g[2]) +
    geom_vline(xintercept=c(0.1), col=pal_g[2]) +
    geom_vline(xintercept=c(0.3), col=pal_g[3]) +
    geom_vline(xintercept=c(0.7), col=pal_g[4]) +
    geom_vline(xintercept=c(0.7), col=pal_g[5]) +
        
    # geom_vline(xintercept=c(0.8), col=pal_g[6]) +
    geom_hline(yintercept=-log10(0.05), col=pal_g[2])+
    geom_text_repel(data = . %>% filter((pval < 0.05) &(Protein %in% c('SAE2', 'FLR1', 'EXO1', 'ECO1', 'RFA1d0213', 'HAC1')| corr > ll | corr < -ll)), box.padding = 0.5, max.overlaps = 40)+
    scale_colour_gradientn(colours = pal)+
    scale_x_continuous(breaks = c(-1, -0.7, -0.3, -0.1, 0, 0.1, 0.3, 0.7, 1))+
    scale_y_continuous(trans='log10')+
    theme_minimal()+
    labs(title=sprintf("Correlation by %s", corr_type), x="Correlation", y="-log10(p-value)")}
  
  return(plot)
}

make_hose_spray_location <- function(corr_type, ll){
  tkach <- read_xls("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/Tkach_refined.xls") %>% 
    rename(Protein = "Standard Name") %>% 
    rename(MMS_Localization_change = "MMS localization change class") %>% 
    mutate(MMS_Localization_change = as_factor(MMS_Localization_change)) %>% 
    select(Protein, MMS_Localization_change)
  
  df <- read.csv(sprintf("D:/ALL_FINAL/Combined_by_perc/median_%s.csv", corr_type)) %>% 
    filter(!(Protein %in% c("PPH22d0210", "LCD1", "RAD5", "RFA1d0210", "LSM3d0210r1", "RRP17d0210", "SLX8", "XRS2d0210", "Contam.")))%>% 
    full_join(tkach, by = c(Protein = "Protein")) %>%
    drop_na()

  pal_g <- wes_palette("AsteroidCity2")
  # pal <- wes_palette("AsteroidCity3", 10, "con")
  plot <- df %>% 
    # filter(-log10(pval) <= 50, corr <0.4, corr > -0.4) %>%
    {ggplot(data = ., aes(x=corr, y=-log10(pval), col=MMS_Localization_change, label=Protein))+
        geom_point() + 
        geom_vline(xintercept=c(-0.3), col=pal_g[1]) +    
        geom_vline(xintercept=c(-0.1), col=pal_g[2]) +
        geom_vline(xintercept=c(0.1), col=pal_g[2]) +
        geom_vline(xintercept=c(0.3), col=pal_g[3]) +
        geom_vline(xintercept=c(0.7), col=pal_g[4]) +
        geom_vline(xintercept=c(0.7), col=pal_g[5]) +
        
        # geom_vline(xintercept=c(0.8), col=pal_g[6]) +
        geom_hline(yintercept=-log10(0.05), col=pal_g[2])+
        # geom_text_repel(data = . %>% filter((pval < 0.05) &(Protein %in% c('SAE2', 'FLR1', 'EXO1', 'ECO1', 'RFA1d0213', 'HAC1')| corr > ll | corr < -ll)), box.padding = 0.5, max.overlaps = 40)+
        # scale_colour_gradientn(colours = pal)+
        scale_x_continuous(breaks = c(-1, -0.7, -0.3, -0.1, 0, 0.1, 0.3, 0.7, 1))+
        scale_y_continuous(trans='log10')+
        theme_minimal()+
        labs(title=sprintf("Correlation by %s", corr_type), x="Correlation", y="-log10(p-value)")}
  ggsave(plot = plot, filename = sprintf('%s_location_hose_spray_pearson.eps', date), width = 12, height = 12, device = cairo_pdf)
  ggsave(plot = plot, filename = sprintf('%s_location_hose_spray_pearson.pdf', date), width = 12, height = 12, device = cairo_pdf)
  ggsave(plot = plot, filename = sprintf('%s_location_hose_spray_pearson.png', date), width = 12, height = 12, device = cairo_pdf)
  return(plot)
}

hose <- make_hose_spray("kendall", 0.3)
hose
# hose <- make_hose_spray("spearman", 0.4)
# hose
hose <- make_hose_spray("pearson", 0.7)
hose <- make_hose_spray_location("pearson", 0.7)
hose

tkach <- read_xls("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/Tkach_refined.xls") %>% 
  rename(Protein = "Standard Name") %>% 
  rename(MMS_Localization_change = "MMS localization change class") %>% 
  select(Protein, MMS_Localization_change)

corr_type = 'pearson'
df_j <- read.csv(sprintf("D:/ALL_FINAL/Combined_by_perc/median_%s.csv", corr_type)) %>% 
  full_join(tkach, by = c(Protein = "Protein")) %>% 
  drop_na()
  
  
  
cross_join(df_j, tkach, join_by(Protein == Standard_Name)) 
  
  
ggsave(plot = hose, filename = sprintf('%s_hose_spray_pearson.eps', date), width = 12, height = 12, device = cairo_pdf)


#Generate the tables with information base

replace_by_p <- function(x, y){
  if(y < 0.05){
    return(x)
  }
  else{
    return("Negligible Correlation")
  }
}

create_corr_table <- function(corr_type){
  df <- read.csv(sprintf("D:/ALL_FINAL/Combined_by_perc/median_%s.csv", corr_type)) %>% 
    filter(!(Protein %in% c("PPH22d0210", "LCD1", "RAD5", "RFA1d0210", "LSM3d0210r1", "RRP17d0210", "SLX8", "XRS2d0210", "Contam."))) %>% 
    filter(!(Protein %like% 'DDC2')) %>% 
    filter(!(Protein %like% 'd0210')) %>%
    # filter((as.integer(diff_hex) < -1 )| (as.integer(diff_hex) > 1)) %>%
    mutate(Correlation_Group = cut(corr, breaks = c(-1, -0.7, -0.5, -0.3, 0.3, 0.5, 0.7, 1),
                                   labels = c("Strong Negative", "Moderate Negative", "Weak Negative","Negligible Correlation", "Weak Positive", "Moderate Postive","Stong Positive"))) %>% 
    as.data.table() %>%
    mutate(Correlation_Group = map2_vec(.x = Correlation_Group, .y = pval, replace_by_p)) %>% 
    mutate(Correlation_Group = fct_relevel(Correlation_Group, c("Strong Negative", "Moderate Negative", "Weak Negative","Negligible Correlation", "Weak Positive", "Moderate Postive","Stong Positive")))
  # filter((diff_hex_dup > 0 )| (diff_hex_dup < 0))
  # mutate(hexile_rank = fct_relevel(hexile_rank, sort(integer())))#Remove the proteins which did not move
  table <- df[, .("Protein Members" = list(unlist(unique(Protein))), N = .N), by = Correlation_Group] %>% 
    arrange(Correlation_Group) %>%
    gt()
  return(table)
}

t <- create_corr_table("pearson")
gtsave(t, sprintf('%s_Abundance_memmbers.pdf', date))
gtsave(t, sprintf('%s_Abundance_memmbers.html', date))

