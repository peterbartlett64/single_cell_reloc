require(pacman)
pacman:: p_load(dplyr, GGally, ggplot2, ggthemes, ggvis, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr, ggpmisc)

setwd("/Users/pcnba/Desktop/Testing_myo/Bit_of_data") #!This is temporary

myo_l <- "x99thPercentile_Diff_background_mKate" #* Automate the search input. This is the other version if only a single is to be done. This could be set as an input in the future
myo_c <- "Myo1_mKa"
variable <- 'myo_smoothed_signal'

data = read.paraquet("entire_df.paraquet")

Peaks_frames <- df %>%
  group(Cell_Barcode) %>%
  findpeaks(.$variable, minpeakdistance = 9, sortstr = FALSE) %>% #This will return values which are labeled as TRUE or FALSE as peak state













Quant_d0222r2p1170300_primary %>%
  filter(Myo1Identity == myo_c) %>%
  ksmooth(x = .$Frame_x, y =.$x99thPercentile_Diff_background_mKate) #%>%

# Quant_d0222r2p1170300_primary$myo_smoothed_signal # This is an attempt to make a new column


Quant_d0222r2p1170300_primary$myo_smoothed_signal %>%
  group_by(Cell_Barcode) %>%
  ksmooth(., a)