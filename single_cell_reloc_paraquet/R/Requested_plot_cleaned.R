
library(httpgd)
library(forcats)
library(beeswarm)
require(pacman)
pacman:: p_load(ggdist, dplyr, GGally, ggplot2, ggthemes, ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr, ggpmisc, arrow)

#, Prepare the graphics to be output to html. Allows for larger figures to be plotted, sped up and keep ability to change the shape
hgd()
hgd_browse()


# wd = "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code"
wd = "D:/Microfluidics/RESULTS_ALL/Most_final_collected/Combined_by_perc/Col_with_Abund"
filename = "For_testing.parquet"


recloc_yet <- function(df, col_name, time_col_name, identity_col_name) {
  flag_col_name <- paste0(col_name, "_flag")
  df[[flag_col_name]] <- 0 #< The downstream function will only accept positive values so it must be 1 and 2 if true>

  for (id in unique(df[[identity_col_name]])){
    idx <- which(df[[identity_col_name]] == id)
    for (f in unique(df[[time_col_name]])){
      frm <- which(df[idx,time_col_name] == f)
      for (i in frm){
        if (sum(frm[1]:i, col_name)> 0){
          df[i, flag_col_name] <- 1
        }
      }
    }
  }
  # for (f in unique(df[[time_col_name]])) {
  #   frm <- which(df[[time_col_name]] == f)

  #   for (c in unique(frm[[identity_col_name]])){
  #     cell <- which(frm[[identity_col_name]] == c)
  #     for (i in cell) {
  #       if (sum(frm[cell[i], col_name]) > 0){
  #         df[i, flag_col_name] <- 1
  #       }
  #     }
  #   }

  # }
  # for (id in unique(df[[identity_col_name]])) {
  #   idx <- which(df[[identity_col_name]] == id)

  #   for (i in idx) {
  #     if (sum(df[idx[1]:i, col_name]) > 0) {
  #     # if (length(which(df[idx[1]:i, col_name] >= 1))>0){
  #       df[i, flag_col_name] <- 1
  #     }
  #   }

    # df[idx, flag_col_name] <- ave(df[idx, flag_col_name],
    #                               df[idx, time_col_name], FUN = cummax)
  return(df)
}


# output_het_violinSwarm <- fucntion (wd, filename){ #, This the function which crates a violin plot swarm to better descibe heterogeneity dynamics
#!This is temporary and will be set based on input from pipeline
setwd(wd) #* Test with a single file
#, Read in the file for testing/run
file = read_parquet(filename) # This should be made variable


reloc_yet_post <- file[which(file$Frames_post_treatment >= 0), ] #< Keep only the data which come after the treamtnet

# reloc_yet_post$refer <- reloc_yet_post$Frames_post_treatment #* Create a copy of the frame series before it is converted to a factor/string
# reloc_yet_post <- reloc_yet_post %>% arrange(Frames_post_treatment) #* Make sure that the Frames are sorted before passing to the _yet function

# file_with_flag <- recloc_yet(reloc_yet_post, "Relocalized", "Frames_post_treatment", "Cell_Barcode")
# file_with_flag <- file_with_flag[which((file_with_flag$refer)%%4 == 0), ]#* Keep every fourth time point - better for figures
# file_with_flag$Frames_post_treatment <- as.factor(file_with_flag$Frames_post_treatment)

reloc_yet_post <- reloc_yet_post[which((reloc_yet_post$refer)%%4 == 0), ]#* Keep every fourth time point - better for figures
reloc_yet_post$Frames_post_treatment <- as.factor(reloc_yet_post$Frames_post_treatment)


# binary_cond <- function(x) {
#   ifelse(x == 1, "Yes", "No")
# }

# ggplot(reloc_yet_post) +
# geom_violin(aes(x = Frames_post_treatment, y = Loc_score)) +
# geom_point(aes(x = Frames_post_treatment, y = Loc_score, shape = as.factor(Relocalized) ,col = binary_cond(Reloc_yet)), position = position_jitterdodge(jitter.width = 0.6, dodge.width = 1)) +
# # geom_point(aes(x = group, y = value2), position = position_jitterdodge(jitter.width = 0.2, dodge.width = -0.75)) +
# labs(x = "Frames After Treatment", y = "Localization Score") +
# geom_hline(yintercept = 1)+
# theme_clean()

reloc_yet_post <- reloc_yet_post %>%
	group_by(Frames_post_treatment, Relocalized)%>%
  	add_count(name = "CurReloc_count")

reloc_yet_post <- reloc_yet_post %>%
	group_by(Frames_post_treatment, Reloc_yet)%>%
  	add_count(name = "CurYet_count")


ggplot(reloc_yet_post) +
geom_violin(aes(x = Frames_post_treatment, y = Loc_score)) +
geom_point(aes(x = Frames_post_treatment, y = Loc_score, shape = as.factor(Relocalized) ,col = binary_cond(Reloc_yet)), position = position_jitterdodge(jitter.width = 0.6, dodge.width = 1)) +
geom_label(aes(x = Frames_post_treatment, y = 2, label = CurReloc_count, fill = Relocalized), position = position_dodge(width = 1))+
labs(x = "Frames After Treatment", y = "Localization Score") +
geom_hline(yintercept = 1)+
theme_clean()


file_with_flag %>%
  ggplot(aes(
    y = Frames_post_treatment,
    x = Loc_score)) +
  #xlim(0.95,1.125) +
  stat_eye(aes(fill = stat(1 < x))) +
  geom_dotsinterval(side = 'bottom', scale = 1.2, height = 0.6, fill = file_with_flag$Loc_score_flag) +
  scale_fill_manual(values = c("#d9d9d9fe", "skyblue")) +
  scale_color_brewer(palette = "Dark2")+
  theme_clean()

file_with_flag %>%
  ggplot(aes(
    y = Unique_Frame,
    x = Loc_score)) +
  #xlim(0.95,1.125) +
  stat_eye(aes(fill = stat(1 < x))) +
  geom_dotsinterval(side = 'bottom', scale = 1.2, height = 0.6, fill = file_with_flag$Loc_score_flag) +
  scale_fill_manual(values = c("#d9d9d9fe", "skyblue")) +
  scale_color_brewer(palette = "Dark2")+
  theme_clean()


beeswarm(file_with_flag$Loc_score, file_with_flag$Frames_post_treatment,
          method = "hex",
          col = file_with_flag$Loc_score_flag)
beeswarm()


# install.packages("beeswarm")
library(beeswarm)

# Data generation
set.seed(1995)
x <- rnorm(300)
g <- sample(c("G1", "G2", "G3"),
            size = 300, replace = TRUE)
z <- as.numeric(factor(sample(c("Yes", "No"),
                              size = 300, replace = TRUE)))

# Bee swarm plot by group
beeswarm(x ~ g,
          pch = 19,
          pwcol = as.numeric(z))

# Legend
legend("topright", legend = c("Yes", "No"),
        col = 1:2, pch = 19)




file_with_flag %>%
  ggplot(aes(
    y = Loc_score,
    x = Frames_post_treatment,
    z = Loc_score_flag
  )) +
  geom_contour()

file_with_flag %>%
  ggplot(aes(
    y = Frames_post_treatment,
    x = Loc_score)) +
  #xlim(0.95,1.125) +
  #stat_eye(aes(fill = stat(1 < x))) +
  geom_dotsinterval(side = 'bottom', scale = 1.2, height = 0.6, fill = file_with_flag$Loc_score_flag, ) +
  scale_fill_manual(values = c("#d9d9d9fe", "darkblue")) +
  scale_color_brewer(palette = "Dark2")+
  theme_clean()

  # return(file_with_flag)}

