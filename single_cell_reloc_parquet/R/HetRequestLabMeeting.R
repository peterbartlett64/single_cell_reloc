# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-02-02
#
# Script Name: Requested plots for lab meeting
#
# Script Description:
#
#
# SETUP ------------------------------------
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)
library(dplyr)
library(hrbrthemes)
library(ggthemes)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(scales)
library(stringr)
library(cowplot)
library(ggforce)
library(statsExpressions)
library(purrr)
# library(ggthemr)
library(extrafont)
library(forcats)
library(ggsignif)
library(ggrepel)
# library(showtext)
# showtext_auto()


# House Keeping -----------------------------------------
# This is the requested plot from Grant. The ggstatsplot will not work so had to re implement
n_samples = 5
date <- Sys.Date()

# Create the discrete color palette for use in this project
npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)
show_col(npg_clrs)

#Set working directory to save the plots.
#Todo: Change this so that it takes input from the global variables param
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

# Module Code --------------------------------------------

#Import fonts from system
# font_import(path="C:/Windows/Fonts", prompt=FALSE)
# choose_font("Arial")
# use_font <- "Arial"


# Load the data. This file is very large so bringing in with the arrow backend and will do filtering
data <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)

#Auto-define the end
# end <- data  %>%
#   group_by(Protein) %>%
#   summarise(max_frame = max(Frames_post_treatment)) %>%
#   summarise(mean_frame = mean(max_frame)) %>%
#   collect()
end = 32
len_end = 32

#Define functions



#, Paired t-test for variable protein. This takes a while to run, so not a great idea


#Unpaired comparison of the ReLocScores
unpaired_comparison  <- function(protein, len_end){
  # Filter the data to the protein of interest
  smaller<- data %>%
    filter(Frames_post_treatment >= 0 & Protein == protein) %>%
    filter(Frames_post_treatment %in% seq(0, len_end, n_samples)) %>%
    collect()
  # Set the variables to factors
  smaller$Relocalized <- as.factor(smaller$Relocalized)
  smaller$Yet <- as.factor(smaller$Yet)
  smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
  # Create the unpaired t-test plot
  unpaired_comparison <- ggbetweenstats(
    data = smaller,
    x = Yet,
    y = Loc_score,
    type = 'robust',
    violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
    point.args = list(size = 0)
  ) +
    geom_sina(smaller, mapping = aes(color = Relocalized, group = Yet), inherit.aes = T, size = 2, alpha = 0.6, scale = 'area')+
    ggtitle(protein)
  # Save the plot
  ggsave(paste0(protein, date, "_unpaired_comparison.pdf"), plot = unpaired_comparison, width = 12, height = 9, device = cairo_pdf)
  ggsave(paste0(protein, date, "_unpaired_comparison.eps"), plot = unpaired_comparison, width = 12, height = 9, device = cairo_pdf)
  # Return the plot
  return(unpaired_comparison) #. This will probably be discarded but still good to return
}


#Building scratch
regular_grouped_stats <- function(protein, len_end){
  smaller<- data %>%
    select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
    filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
    # This is workable right now
    collect() %>%
    filter(Frames_post_treatment %in% seq(0, len_end, length.out = n_samples)) %>%
    # Set the variables to factors
    mutate(Relocalized = as_factor(Relocalized)) %>%
    mutate(Yet = as_factor(Yet)) %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment))


  plot <- ggplot(smaller, aes(x = Frames_post_treatment, y = Loc_score, group = interaction(Frames_post_treatment, Yet))) +
    geom_sina(smaller, mapping = aes(color = Relocalized), inherit.aes = T, size = 2, alpha = 0.6, scale = 'width')+
    geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.6,
                 position = position_dodge(0.9), color = 'black')+
    geom_violin(fill = NA, color = "black", width = 0.6, linewidth = 0.4,
    position = position_dodge(0.9)) +
    geom_point(stat = "summary", size = 5, color = "#8a0f00",
               position = position_dodge(0.9), fun = median) +
    geom_label_repel(stat = "summary", fun = median, size = 3.5,
                     aes(label = paste0("hat(mu)*scriptstyle(median)==",
                                        round(after_stat(y), 2))),
                     parse = TRUE, position = position_dodge(0.9)) +
    geom_signif(y_position = 1.5, xmin = 1:5 - 0.2, xmax = 1:5 + 0.2,
                annotations = scales::pvalue(sapply(split(smaller, smaller$Frames_post_treatment),
                                                    \(x) wilcox.test(Loc_score~Yet, x)$p.value),
                                             add_p = TRUE)) +
    scale_color_npg() +
    scale_fill_npg()+
    theme_ipsum()+
    theme(axis.title = element_text(face = 2),
        legend.position = "bottom",
        axis.text.y.right = element_blank())+
    scale_y_continuous(name = "ReLocScore", labels = c(1, 1.2,1.4,1.6, 1.8))+
  ggtitle(protein)

  ggsave(paste0(protein, date, "violin_het_yet.pdf"), plot = plot, width = 12, height = 9, device = cairo_pdf)
  ggsave(paste0(protein, date, "violin_het_yet.eps"), plot = plot, width = 12, height = 9, device = cairo_pdf)
  
  return(plot)
}


# The below is similar but based on the ggbetweenstats backbone. It is ugly because the dots are in front
unpaired_group_comparison  <- function(protein, len_end){
  # Filter the data to the protein of interest
  smaller<- data %>%
    select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
    filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
    # This is workable right now
    collect() %>%
    filter(Frames_post_treatment %in% seq(0, 32, length.out = n_samples)) %>%
    # Set the variables to factors
    mutate(Relocalized = as_factor(Relocalized)) %>%
    mutate(Yet = as_factor(Yet)) %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment))

  # Create the grouped the plot
  group_violin <- grouped_ggbetweenstats(
    data = smaller,
    x = Yet,
    y = Loc_score,
    type = 'robust',
    grouping.var = Frames_post_treatment,
    violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
    point.args = list(size = 0)) +
    geom_sina(smaller, mapping = aes(color = Relocalized, group = Yet), inherit.aes = T, size = 2, alpha = 0.6, scale = 'area')+
    ggtitle(protein)
  # Save the plot
  ggsave(paste0(protein, date, "_grouped_comparison_ggbetweenBackbone.pdf"), plot = group_violin, width = 12, height = 9, device = cairo_pdf)
  ggsave(paste0(protein, date, "_grouped_comparison_ggbetweenBackbone.eps"), plot = group_violin, width = 12, height = 9, device = cairo_pdf)
  # Return the plot
  # return(unpaired_comparison) #. This will probably be discarded but still good to return
}


between_sina <- function(protein, len_end){
  smaller<- data %>%
    select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
    filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
    # This is workable right now
    collect() %>%
    filter(Frames_post_treatment %in% seq(0, 32, length.out = n_samples)) %>%
    # Set the variables to factors
    mutate(Relocalized = as_factor(Relocalized)) %>%
    mutate(Yet = as_factor(Yet)) %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment))

  between <- ggbetweenstats(
    data = smaller,
    x = Yet,
    y = Loc_score,
    type = 'robust',
    # points.color.palette = "Set1", # You can specify your desired color palette here
    # boxplot.args = list(fill = "white", width = 'area'),
    # violin.args = list(width = 0),
    violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
    point.args = list(size = 0),
    pairwise.display = TRUE
    # point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE)
    # title = "Group Comparison",
    # xlab = "Group",
    # ylab = "Variable"
  ) +
    # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
    # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
    geom_sina(smaller, mapping = aes(color = Relocalized, group = Yet), inherit.aes = T, size = 2, alpha = 0.2, scale = 'area')+
    theme_ipsum()+
    facet_grid(~Frames_post_treatment)+
    scale_color_aaas()+
    ggtitle(protein)
  ggsave(paste0(protein, date, "facet_sina_score.pdf"), plot = between, width = 12, height = 9, device = cairo_pdf)
  ggsave(paste0(protein, date, "facet_sina_score.eps"), plot = between, width = 12, height = 9, device = cairo_pdf)
}


paired_t_test <- function(protein, len_end){
  smaller <- data %>%
    select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
    filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
    # This is workable right now
    collect() %>%
    filter(Frames_post_treatment %in% seq(0, len_end, length.out = n_samples)) %>%
    # Set the variables to factors
    mutate(Relocalized = as_factor(Relocalized)) %>%
    mutate(Yet = as_factor(Yet)) %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment))
  
  smaller <- smaller %>% 
    filter(Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
           Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode)) %>%
    arrange(Frames_post_treatment, Cell_Barcode) %>% 
    mutate(Cell_Barcode = as_factor(Cell_Barcode)) %>%
    mutate(id = as.numeric(Cell_Barcode))
  
  paired_t_test <- ggwithinstats(smaller,
    #Drop duplicate frames_post_treatment
    # data = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
    #               Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode), )%>% 
    #   arrange(Frames_post_treatment, Cell_Barcode) %>% 
    #   mutate(Cell_Barcode = as_factor(Cell_Barcode)) %>%
    #   mutate(Cell_Barcode = as.numeric(Cell_Barcode)),
    x = Frames_post_treatment,
    y = Loc_score,
    type = "robust",
    withinVars = "id",
    pairwise_comparisons = TRUE,
    # points.color.palette = "Set1", # You can specify your desired color palette here
    # boxplot.args = list(fill = "white", width = 'area'),
    violin.args = list(width = 0),
    # violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
    point.args = list(size = 0)
    # point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE)
    # title = "Group Comparison",
    # xlab = "Group",
    # ylab = "Variable"
  )+
  # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
  # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
  geom_sina(data  = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
                           Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode)), mapping = aes(color = Yet, group = Frames_post_treatment), inherit.aes = T, size = 2, alpha = 0.6, scale = 'count')+
  scale_color_npg()+
  scale_fill_npg()+
  ggtitle(protein)

  ggsave(paste0(protein, date, "_paired_t_test.pdf"), plot = paired_t_test, width = 12, height = 9, device = cairo_pdf)
  ggsave(paste0(protein, date, "_paired_t_test.eps"), plot = paired_t_test, width = 12, height = 9, device = cairo_pdf)
  # Return the plot
  return(paired_t_test) #. This will probably be discarded but still good to return
}

#* This will take a long time to run
paired_swarm <- function(protein, len_end){
  smaller <- data %>%
    select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
    filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
    # This is workable right now
    collect() %>%
    filter(Frames_post_treatment %in% seq(0, len_end, length.out = n_samples)) %>%
    # Set the variables to factors
    mutate(Relocalized = as_factor(Relocalized)) %>%
    mutate(Yet = as_factor(Yet)) %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment))

  #* This is a new version with a sina plot
  paired_t_test_swarm <- ggwithinstats(
    #Drop duplicate frames_post_treatment
    data = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
                  Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode)),
    x = Frames_post_treatment,
    y = Loc_score,
    type = 'robust',
    # points.color.palette = "Set1", # You can specify your desired color palette here
    # boxplot.args = list(fill = "white", width = 'area'),
    violin.args = list(width = 0),
    # violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
    point.args = list(size = 0)
    # point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE)
    # title = "Group Comparison",
    # xlab = "Group",
    # ylab = "Variable"
    # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
  )+
    geom_sina(data  = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
                             Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode))
              , mapping = aes(color = Yet, group = Frames_post_treatment), inherit.aes = T, size = 2, alpha = 0.6, scale = 'count')+
    scale_color_npg()+
    scale_fill_npg()+
    ggtitle(protein)
  # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
  ggsave(paste0(protein, date, "paired_t_test_swarm.pdf"), plot = paired_t_test_swarm, width = 12, height = 9, device = cairo_pdf)
  ggsave(paste0(protein, date, "paired_t_test_swarm.eps"), plot = paired_t_test_swarm, width = 12, height = 9, device = cairo_pdf)
}

#####
facet_attempt <- function (protein, len_end){
  smaller <- data %>%
    select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
    filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
    # This is workable right now
    collect() %>%
    filter(Frames_post_treatment %in% seq(0, len_end, length.out = n_samples)) %>%
    # Set the variables to factors
    mutate(Relocalized = as_factor(Relocalized)) %>%
    mutate(Yet = as_factor(Yet)) %>%
    mutate(Frames_post_treatment = as_factor(Frames_post_treatment))

  facet_loc_graph <- smaller %>%
    group_by(Frames_post_treatment, Yet) %>%
    {ggplot(smaller, aes(x = Yet, y = Loc_score)) +
        # geom_violin(aes(fill = Yet), scale = 'width')+
        # geom_point(aes(fill = Relocalized), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.2)) +
        geom_boxplot(aes(group = Yet, fill = Yet), notch = T, alpha = 0.4)+
        geom_violin(aes(group = Yet), scale = 'width', color = npg_clrs[4], fill = NA) +
        geom_sina(aes(color = Relocalized, group = Yet), scale = 'width', jitter_y = F) +
        scale_fill_npg()+
        scale_color_npg()+
        # scale_y_log10(limits = c(0.9,2)) +
        # theme_ipsum() +
        # theme_minimal()+
        theme(legend.position = "bottom") +
        labs(title = "ECO1", y = "Loc_score")+
        scale_x_discrete("Frames_post_treatment", labels = c(
          "0" = "Has not",
          "1" = "Has relocalized far"
        ))+
        theme_light(base_family = "Arial")+
        # labels = c(paste0("Not Yet (", .$CurrNotYet, ")"), paste0("Yes (", .$CurrYet), ")"))+
        # subtitle = results_data$expression[[1]] +
        facet_grid(rows = NULL, cols = vars(Frames_post_treatment), scale = "free_x")+
        ggtitle(protein)}

  ggsave(paste0(protein, date, "_loc_score.pdf"), plot = facet_loc_graph, width = 12, height = 9, device = cairo_pdf)
  ggsave(paste0(protein, date, "_loc_score.eps"), plot = facet_loc_graph, width = 12, height = 9, device = cairo_pdf)
}



# paired_t_test(protein = 'EXO1', len_end = 32)
# paired_t_test(protein = 'FLR1', len_end = 32)
# paired_t_test(protein = 'ECO1', len_end = 32)

unpaired_comparison(protein = 'EXO1', len_end = 32)
unpaired_comparison(protein = 'FLR1', len_end = 32)
unpaired_comparison(protein = 'ECO1', len_end = 32)


unpaired_comparison(protein = 'SGS1', len_end = 32)
unpaired_comparison(protein = 'SLX8', len_end = 32)
unpaired_comparison(protein = 'MCM2', len_end = 32)
unpaired_comparison(protein = 'MRT4', len_end = 32)
unpaired_comparison(protein = 'RNR1d0216r1', len_end = 32)
unpaired_comparison(protein = 'RNR1d0222r1', len_end = 32)



regular_grouped_stats(protein = 'EXO1', len_end = 32)
regular_grouped_stats(protein = 'FLR1', len_end = 32)
regular_grouped_stats(protein = 'ECO1', len_end = 32)

unpaired_group_comparison(protein = 'EXO1', len_end = 32)
unpaired_group_comparison(protein = 'FLR1', len_end = 32)
unpaired_group_comparison(protein = 'ECO1', len_end = 32)


between_sina(protein = 'EXO1', len_end = 32)
between_sina(protein = 'FLR1', len_end = 32)
between_sina(protein = 'ECO1', len_end = 32)

between_sina(protein = 'SGS1', len_end = 32)
between_sina(protein = 'SLX8', len_end = 32)
between_sina(protein = 'MCM2', len_end = 32)
between_sina(protein = 'MRT4', len_end = 32)
between_sina(protein = 'RNR1d0216r1', len_end = 32)
between_sina(protein = 'RNR1d0222r1', len_end = 32)


paired_t_test(protein = 'EXO1', len_end = 32)
paired_t_test(protein = 'FLR1', len_end = 32)
paired_t_test(protein = 'ECO1', len_end = 32)

paired_t_test(protein = 'SGS1', len_end = 32)
paired_t_test(protein = 'SLX8', len_end = 32)
paired_t_test(protein = 'MCM2', len_end = 32)
paired_t_test(protein = 'MRT4', len_end = 32)
paired_t_test(protein = 'RNR1d0216r1', len_end = 32)
paired_t_test(protein = 'RNR1d0222r1', len_end = 32)


paired_swarm(protein = 'EXO1', len_end = 32)
paired_swarm(protein = 'FLR1', len_end = 32)
paired_swarm(protein = 'ECO1', len_end = 32)

paired_swarm(protein ='SGS1', len_end = 32)
paired_swarm(protein ='SLX8', len_end = 32)
paired_swarm(protein ='MCM2', len_end = 32)
paired_swarm(protein ='MRT4', len_end = 32)
paired_swarm(protein ='RNR1d0216r1', len_end = 32)
paired_swarm(protein ='RNR1d0222r1', len_end = 32)



#* No longer using these
# facet_attempt(protein = 'EXO1', len_end = 32)
# facet_attempt(protein = 'FLR1', len_end = 32)
# facet_attempt(protein = 'ECO1', len_end = 32)



# smaller<- data %>%
#   select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
#   filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
#   # This is workable right now
#   collect() %>%
#   filter(Frames_post_treatment %in% seq(0, len_end, length.out = n_samples)) %>%
#   # Set the variables to factors
#   mutate(Relocalized = as_factor(Relocalized)) %>%
#   mutate(Yet = as_factor(Yet)) %>%
#   mutate(Frames_post_treatment = as_factor(Frames_post_treatment)) %>%
#   # subset(Unique_Frame, Loc_score, Does, No_yet, Yes_yet, pres_end, CurrNot,
#   #               CurrYes, currProportion, RelocVelocity,RelocAcceleration, CurrNotYet,
#   #               CurrYet, currYetProportion, YetVelocity, YetAcceleration, countDoesNot,
#   #               countDoes, DoesProportion, DoesVelocity, Does_FinDiff) %>%
#   # filter(Frames_post_treatment %in% c(0, 8, 16, 32, 40)) %>%
#



  #, This is the simplified call with betweenstats defaults

#* This is the paired t-test between the values



#* This will take a long time to run
# paired_swarm <- function(protein, len_end){
#   smaller <- data %>%
#     select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
#     filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
#     # This is workable right now
#     collect() %>%
#     filter(Frames_post_treatment %in% seq(0, len_end, length.out = n_samples)) %>%
#     # Set the variables to factors
#     mutate(Relocalized = as_factor(Relocalized)) %>%
#     mutate(Yet = as_factor(Yet)) %>%
#     mutate(Frames_post_treatment = as_factor(Frames_post_treatment))
#
#   #* This is a new version with a sina plot
#   paired_t_test_swarm <- ggwithinstats(
#     #Drop duplicate frames_post_treatment
#     data = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
#                   Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode)),
#     x = Frames_post_treatment,
#     y = Loc_score,
#     # points.color.palette = "Set1", # You can specify your desired color palette here
#     # boxplot.args = list(fill = "white", width = 'area'),
#     violin.args = list(width = 0),
#     # violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
#     point.args = list(size = 0)
#     # point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE)
#     # title = "Group Comparison",
#     # xlab = "Group",
#     # ylab = "Variable"
#     # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
#   )+
#     geom_sina(data  = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
#                              Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode))
#               , mapping = aes(color = Yet, group = Frames_post_treatment), inherit.aes = T, size = 2, alpha = 0.6, scale = 'count')
#   # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
#   ggsave(paste0(protein, date, "paired_t_test_swarm.pdf"), plot = paired_t_test_swarm, width = 12, height = 9, device = cairo_pdf)
# }


#####
# facet_attempt <- function (proteien, len_end){
#   smaller <- data %>%
#     select(Cell_Barcode, Frames_post_treatment, Loc_score, MMS_localization_class, Protein, Relocalized, Yet) %>%
#     filter((Frames_post_treatment >= 0 ) & (Protein == protein)) %>%
#     # This is workable right now
#     collect() %>%
#     filter(Frames_post_treatment %in% seq(0, len_end, length.out = n_samples)) %>%
#     # Set the variables to factors
#     mutate(Relocalized = as_factor(Relocalized)) %>%
#     mutate(Yet = as_factor(Yet)) %>%
#     mutate(Frames_post_treatment = as_factor(Frames_post_treatment))
#
#   facet_loc_graph <- smaller %>%
#     group_by(Frames_post_treatment, Yet) %>%
#     {ggplot(smaller, aes(x = Yet, y = Loc_score)) +
#         # geom_violin(aes(fill = Yet), scale = 'width')+
#         # geom_point(aes(fill = Relocalized), position = position_jitterdodge(jitter.width = 2, dodge.width = 0.2)) +
#         geom_boxplot(aes(group = Yet, fill = Yet), notch = T, alpha = 0.4)+
#         geom_violin(aes(group = Yet), scale = 'width', color = npg_clrs[4], fill = NA) +
#         geom_sina(aes(color = Relocalized, group = Yet), scale = 'width', jitter_y = F) +
#         scale_fill_npg()+
#         scale_color_npg()+
#         # scale_y_log10(limits = c(0.9,2)) +
#         # theme_ipsum() +
#         # theme_minimal()+
#         theme(legend.position = "bottom") +
#         labs(title = "ECO1", y = "Loc_score")+
#         scale_x_discrete("Frames_post_treatment", labels = c(
#           "0" = "Has not",
#           "1" = "Has relocalized far"
#         ))+
#         theme_light(base_family = "Arial")+
#         # labels = c(paste0("Not Yet (", .$CurrNotYet, ")"), paste0("Yes (", .$CurrYet), ")"))+
#         # subtitle = results_data$expression[[1]] +
#         facet_grid(rows = NULL, cols = vars(Frames_post_treatment), scale = "free_x")}
#
#   ggsave(paste0(protein, date, "_loc_score.pdf"), plot = facet_loc_graph, width = 12, height = 9, device = cairo_pdf)
# }




### Testing something
# smaller %>%
#   {ggbetweenstats(
#     data = smaller,
#     x    = Yet,
#     y    = Loc_score,
#     type = "np",
#     var.equal = FALSE,
#     outlier.tagging = FALSE,
#     point.args = list(data = smaller, aes(color = smaller$Relocalized, group = smaller$Yet)))}




# # Plot using between and change point color by another variable
# ggbetweenstats(
#   data = data,
#   x = group,
#   y = variable,
#   # points.color.palette = "Set1", # You can specify your desired color palette here
#   # boxplot.args = list(fill = "white", width = 'area'),
#   # violin.args = list(width = 0),
#   violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
#   point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
#                       0.4, size = 3, stroke = 0, na.rm = TRUE),
#   title = "Group Comparison",
#   xlab = "Group",
#   ylab = "Variable"
# ) +
#   # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
#   # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
#   geom_sina(data, mapping = aes(color = color_variable, group = group), inherit.aes = T, size = 2, alpha = 0.6, scale = 'area')
# ###
#





# geom_sina(aes(color = Relocalized, group = Yet), scale = 'width') +
# scale_fill_npg()+
# scale_color_npg()+
# # scale_y_log10(limits = c(0.9,2)) +
# theme_ipsum() +
# theme(legend.position = "top") +
# labs(title = "ECO1", y = "Loc_score")
# npst


# ggwithinstats(
#   data = smaller,
#   x    = Frames_post_treatment,
#   y    = Loc_score,
#   type = "robust"
# )
#
#
#
# data_mods <- data %>%
#   group_by(Frames_post_treatment) %>%
#   mutate()
#   mutate(Yet_count = n(filter(., Yet == "Yes"))) %>%
#   mutate(Yet_display = paste0(, " (n = ", Yet_count, ")")) %>%
#   mutate(yet)
# data_mod
#
#
# rq_plot <- ggpt
#
#
#
# facet_grid(~Frames_post_treatment)
#
#
# paired_t_test <- function(protein, len_end){
#   # Filter the data to the protein of interest
#   smaller<- data%>%
#     filter(Frames_post_treatment >= 0 & Protein == protein) %>%
#     filter(Frames_post_treatment %in% seq(0, len_end, n_samples)) %>%
#     collect()
#   # Set the variables to factors
#   smaller$Relocalized <- as.factor(smaller$Relocalized)
#   smaller$Yet <- as.factor(smaller$Yet)
#   smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
#   # Create the paired t-test plot
#   ggwithinstats(
#     #Drop duplicate frames_post_treatment
#     data = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
#                   Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode)),
#     x = Frames_post_treatment,
#     y = Loc_score,
#     # points.color.palette = "Set1", # You can specify your desired color palette here
#     # boxplot.args = list(fill = "white", width = 'area'),
#     violin.args = list(width = 0),
#     # violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
#     point.args = list(size = 0)
#     # point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE)
#     # title = "Group Comparison",
#     # xlab = "Group",
#     # ylab = "Variable"
#     # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
#   )+
#     geom_sina(data  = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
#                              Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode)), mapping = aes(color = Yet, group = Frames_post_treatment), inherit.aes = T, size = 2, alpha = 0.6, scale = 'count')
#   ggsave(paste0(protein, date, "_paired_t_test.pdf"), plot = paired_t_test, width = 12, height = 9, device = cairo_pdf)
#   ggsave(paste0(protein, date, "_paired_t_test.eps"), plot = paired_t_test, device = cairo_pdf)
#   # Return the plot
#   return(paired_t_test) #. This will probably be discarded but still good to return
# }
