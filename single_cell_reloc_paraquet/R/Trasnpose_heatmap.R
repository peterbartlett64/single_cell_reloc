require(pacman)
pacman:: p_load(here, dplyr, GGally, ggplot2, ggthemes, ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr, ggpmisc, data.table)
here() # Set



setwd("D:/Microfluidics/RESULTS_ALL/Most_final_collected/Combined_by_perc/Col_with_Abund")#!This is temporary
masters = fread("Subset_concat.parquet", sep = ",")
install.packages("here")
master_sub <- subset(masters, select = c(Protein, Loc_score, Abundance, z_score_Loc, z_score_Abund, Pearson_coeff_Loc_Abund))
master_sub_u <- master_sub[!duplicated(master_sub)]
rownames(master_sub_u) <- master_sub_u$Protein

heatmap(master_sub)


plot(df_corr)

dft_corr_lib = t(df_corr_lib)

#for_matrix = subset(df_corr_lib, select = c(Protein, Pearson_coeff_Loc_Abund))

Abundance_heatmap <- ggplot(data = master_sub, mapping = aes(x = Protein, y = Abundance,
  fill = Pearson_coeff_Loc_Abund)) +
  geom_tile() +
  xlab(label = "Protein") +
  ylab(label = "Frame")

Abundance_heatmap
