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
# library(showtext)
# showtext_auto()

# Module Code--------------------------------------------
# This is the requested plot from Grant. The ggstats plot will not

# I don't know why I have a number of samples here. I think it was for a different plot
# n_samples = 5

date <- Sys.Date()

#Set working directory to save the plots.
#Todo: Change this so that it takes input from the global variables param
setwd("D:/Second Mind/Academic/Project Stuff/Figures")
# font_import(path="C:/Windows/Fonts", prompt=FALSE)
# choose_font("Arial")
# use_font <- "Arial"


# Load the data. This file is very large so bringing in with the arrow backend and will do filtering
data <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F)

# Create the discrete color palette for use in this project
npg_clrs <-  pal_npg("nrc", alpha = 0.7)(4)
show_col(npg_clrs)

# rq_hetviolin <- function(protein){
#
#
#
# }

#, Paired t-test for variable protein
paired_t_test <- function(protein, count){
  # Filter the data to the protein of interest
  smaller<- data %>% filter(Frames_post_treatment >= 0 & Protein == protein) %>% collect()
  # Set the variables to factors
  smaller$Relocalized <- as.factor(smaller$Relocalized)
  smaller$Yet <- as.factor(smaller$Yet)
  smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
  # Create the paired t-test plot
  paired_t_test <- ggwithinstats(
    #Drop duplicate frames_post_treatment
    data = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
                  Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode), ),
    x = Frames_post_treatment,
    y = Loc_score
  )
  # Save the plot
  ggsave(paste0(protein, "_paired_t_test.pdf"), plot = paired_t_test, width = 4, height = 3, units = 'in', device = cairo_pdf)
  # Return the plot
  return(paired_t_test) #. This will probably be discarded but still good to return
}

unpaired_comparison  <- function(protein, count){
  # Filter the data to the protein of interest
  smaller<- data %>% filter(Frames_post_treatment >= 0 & Protein == protein) %>% collect()
  # Set the variables to factors
  smaller$Relocalized <- as.factor(smaller$Relocalized)
  smaller$Yet <- as.factor(smaller$Yet)
  smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
  # Create the unpaired t-test plot
  unpaired_comparison <- ggbetweenstats(
    data = smaller,
    x = Yet,
    y = Loc_score,
    violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
    point.args = list(size = 0)
  ) +
    geom_sina(smaller, mapping = aes(color = Relocalized, group = Yet), inherit.aes = T, size = 2, alpha = 0.6, scale = 'area')
  # Save the plot
  ggsave(paste0(protein, "_unpaired_comparison.pdf"), plot = unpaired_comparison, width = 4, height = 3, units = 'in', device = cairo_pdf)
  # Return the plot
  return(unpaired_comparison) #. This will probably be discarded but still good to return
}


#. In the future this can be set as an automated function subset
smaller<- data %>%
  filter(Frames_post_treatment >= 0 & Protein == 'FLR1') %>%
  # filter(Frames_post_treatment %in% seq(1, n_samples, 5)) %>%
  # subset(Unique_Frame, Loc_score, Does, No_yet, Yes_yet, pres_end, CurrNot,
  #               CurrYes, currProportion, RelocVelocity,RelocAcceleration, CurrNotYet,
  #               CurrYet, currYetProportion, YetVelocity, YetAcceleration, countDoesNot,
  #               countDoes, DoesProportion, DoesVelocity, Does_FinDiff) %>%
  # filter(Frames_post_treatment %in% c(0, 8, 16, 32, 40)) %>%
  collect()
  # smaller(ImageID, Loc_score, Does, No_yet, Yes_yet, pres_end, CurrNot,
  #        CurrYes, currProportion, RelocVelocity,RelocAcceleration, CurrNotYet,
  #        CurrYet, currYetProportion, YetVelocity, YetAcceleration, countDoesNot,
  #        countDoes, DoesProportion, DoesVelocity, Does_FinDiff)

smaller$Relocalized <- as.factor(smaller$Relocalized)
smaller$Yet <- as.factor(smaller$Yet)
smaller$Frames_post_treatment <- as.factor(smaller$Frames_post_treatment)
# smaller <- filter(smaller, Frames_post_treatment == "20")




#, This is the simplified call with betweenstats defaults
paired_t_test <- ggwithinstats(
  #Drop duplicate frames_post_treatment
  data = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
                Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode), ),
  x = Frames_post_treatment,
  y = Loc_score
  # points.color.palette = "Set1", # You can specify your desired color palette here
  # boxplot.args = list(fill = "white", width = 'area'),
  # violin.args = list(width = 0),
  # violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
  # point.args = list(size = 0)
  # point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE)
  # title = "Group Comparison",
  # xlab = "Group",
  # ylab = "Variable"
)
  # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
  # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
  # geom_sina(smaller, mapping = aes(color = Relocalized, group = Yet), inherit.aes = T, size = 2, alpha = 0.6, scale = 'area')
ggsave("paired_t_test.pdf", plot = paired_t_test, width = 4, height = 3, units = 'in', device = cairo_pdf)


#* This is a new version with a sina plot
paired_t_test_swarm <- ggwithinstats(
  #Drop duplicate frames_post_treatment
  data = filter(smaller, Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == tail(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode),
                Cell_Barcode %in% (smaller[smaller$Frames_post_treatment == head(smaller$Frames_post_treatment, n = 1),]$Cell_Barcode)),
  x = Frames_post_treatment,
  y = Loc_score,
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
          , mapping = aes(color = Yet, group = Frames_post_treatment), inherit.aes = T, size = 2, alpha = 0.6, scale = 'count')
# geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
ggsave("paired_t_test_swarm.pdf", plot = paired_t_test_swarm, width = 4, height = 3, units = 'in', device = cairo_pdf)




ggbetweenstats(
  data = smaller,
  x = Yet,
  y = Loc_score,
  # points.color.palette = "Set1", # You can specify your desired color palette here
  # boxplot.args = list(fill = "white", width = 'area'),
  # violin.args = list(width = 0),
  violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
  point.args = list(size = 0)
  # point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha = 0.4, size = 3, stroke = 0, na.rm = TRUE)
  # title = "Group Comparison",
  # xlab = "Group",
  # ylab = "Variable"
) +
  # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
  # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
  geom_sina(smaller, mapping = aes(color = Relocalized, group = Yet), inherit.aes = T, size = 2, alpha = 0.6, scale = 'area')






#####
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
    facet_grid(rows = NULL, cols = vars(Frames_post_treatment), scale = "free_x")}

facet_loc_graph

ggsave("ECO1_loc_score.pdf", plot = facet_loc_graph, width = 4, height = 3, units = 'in', device = cairo_pdf)


  # theme(axis.text.x = element_text(angle = 45, hjust = 1))

### Testing something
smaller %>%
  {ggbetweenstats(
    data = smaller,
    x    = Yet,
    y    = Loc_score,
    type = "np",
    var.equal = FALSE,
    outlier.tagging = TRUE,
    point.args = list(data = smaller, aes(color = smaller$Relocalized, group = smaller$Yet)))}


data <- data.frame(
  group = rep(c("A", "B", "C"), each = 30),
  variable = rnorm(90),
  color_variable = sample(c("X", "Y", "Z"), 90, replace = TRUE)
)

# Plot using ggbetween and change point color by another variable
ggbetweenstats(
  data = data,
  x = group,
  y = variable,
  # points.color.palette = "Set1", # You can specify your desired color palette here
  # boxplot.args = list(fill = "white", width = 'area'),
  # violin.args = list(width = 0),
  violin.args = list(alpha = 0.2, na.rm = TRUE, scale = 'area'),
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
                      0.4, size = 3, stroke = 0, na.rm = TRUE),
  title = "Group Comparison",
  xlab = "Group",
  ylab = "Variable"
) +
  # geom_point(data, mapping = aes(colour = color_variable), inherit.aes = TRUE, position = "jitter", show.legend = TRUE)+
  # geom_violin(data, mapping = aes(fill = group), inherit.aes = T, scale = 'count', alpha = 0.2)+
  geom_sina(data, mapping = aes(color = color_variable, group = group), inherit.aes = T, size = 2, alpha = 0.6, scale = 'area')
###


# geom_sina(aes(color = Relocalized, group = Yet), scale = 'width') +
  # scale_fill_npg()+
  # scale_color_npg()+
  # # scale_y_log10(limits = c(0.9,2)) +
  # theme_ipsum() +
  # theme(legend.position = "top") +
  # labs(title = "ECO1", y = "Loc_score")
npst


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
