library(tidyverse)
library(moments)
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
library(ggstatsplot)


get_max <- read_parquet("D:/ALL_FINAL/Combined_by_perc/new_percs.parquet", as_data_frame =T)%>% 
  group_by(Protein) %>% 
  summarise(max = max(Yet_perc)) %>%
  ungroup() %>% 
  write_csv('.csv')




test <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Quant_ALL.parquet", as_data_frame = F)%>% 
  select(Frames_post_treatment, Cell_Barcode, x95thPercentile_GFP_RAW, x99thPercentile_GFP_RAW, x95thPercentile_norm_OBJ_Median_GFP, x99thPercentile_norm_OBJ_Median_GFP) %>% 
  collect() %>% 
  group_by(Frames_post_treatment, Protein) %>%
  reframe(Upper = median(x99thPercentile_GFP), Lower = median(x95thPercentile_GFP))

UP <- test %>% 
  filter(Frames_post_treatment < 10, Frames_post_treatment > -20) %>% 
  select(Upper, Frames_post_treatment)
LOW <- test %>% 
  filter(Frames_post_treatment < 10, Frames_post_treatment > -20) %>% 
  select(Lower, Frames_post_treatment)

summary_data <- UP %>%
  group_by(Frames_post_treatment) %>%
  summarise(
    mean = mean(Upper),
    lower_ci = mean - 1.96 * sd(Upper) / sqrt(n()),
    upper_ci = mean + 1.96 * sd(Upper) / sqrt(n())
  )

# Create line plot with confidence intervals
ggplot(summary_data, aes(x = Frames_post_treatment, y = mean, color = series)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = series), alpha = 0.2) +
  labs(x = "X", y = "Y", color = "Series") +
  theme_minimal()


test_count <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Final_wAbund.parquet", as_data_frame = F)%>% 
  select(Frames_post_treatment, Cell_Barcode, Protein) %>% 
  filter(Protein == 'LSM12') %>% 
  filter(Frames_post_treatment == 0) %>% 
  collect() %>% 
  group_by(Protein) %>% 
  summarise(n = n()) %>% 
  arrange(n)

binom.test(99,100,p=1, alternative = 'less', conf.level = 0.95)




 testing <- test %>% 
  filter(Frames_post_treatment < 5)

# ggplot(data = testing, aes(x = Frames_post_treatment, y = Lower, group = Protein)) +
#   # geom_line()+
#   geom_smooth(aes(group = 1), linewidth = 2)


before_treat <- ggwithinstats(testing,
              x = Frames_post_treatment,
              y = Upper)



tested <- test %>% 
  select(Cell_Barcode, Unique_Frame,Upper, Lower) %>% 
  collect() %>% 
  summarise(Upper = mean())
  summarise(Upper = sum(Upper), Lower = sum(Lower))


df_library <- read_parquet("D:/ALL_FINAL/Combined_by_perc/penetrance_updated_trimmed.parquet", as_data_frame = T) %>%
  select(Percentage_reloc, updated_yet_perc)%>%
  summarise(PF_m = mean(Percentage_reloc, na.rm = T), TC_m = mean(updated_yet_perc, na.rm = T),
            PF_s = skewness(Percentage_reloc, na.rm = T), TC_s = skewness(updated_yet_perc, na.rm = T))


df_ddc2 <- read_parquet("D:/ALL_FINAL/Combined_by_perc/penetrance_updated_trimmed.parquet", as_data_frame = T) %>%
  select(Percentage_reloc, updated_yet_perc, Protein)%>%
  filter(Protein %like% 'DDC2')%>% 
  # filter(!(Protein %in% c("DDC2d0218r2p60KO", "DDC2d0224r1", "DDC2d0218r2p80KO", "DDC2d0218r2p20KO")))%>% #PF outliers
  # filter(!(Protein %in% c("DDC2d0223r1", "DDC2d0220r1", "DDC2d0218r2p100KO", "DDC2d0222r2")))%>% #TC outliers
  reframe(PF_m = median(Percentage_reloc, na.rm = T), PF_s = skewness(Percentage_reloc, na.rm = T), PF_var = var(Percentage_reloc, na.rm = T), PF_sd = sd(Percentage_reloc, na.rm = T),
            TC_m = median(updated_yet_perc, na.rm = T), TC_s = skewness(updated_yet_perc, na.rm = T), TC_var = var(updated_yet_perc, na.rm = T), TC_sd = sd(updated_yet_perc, na.rm = T),
            count = n(), tc_iqr = IQR(updated_yet_perc, na.rm = T), pf_iqr = IQR(Percentage_reloc, na.rm = T),
          PF_1stq = quantile(Percentage_reloc, probs = c(0.25), na.rm = T), PF_3rdq = quantile(Percentage_reloc, probs = c(0.75), na.rm = T),
          TC_1stq = quantile(updated_yet_perc, probs = c(0.25), na.rm = T), TC_3rdq = quantile(updated_yet_perc, probs = c(0.75), na.rm = T))
# df_reps <- read_parquet("D:/ALL_FINAL/Combined_by_perc/penetrance_updated_trimmed.parquet", as_data_frame = T) %>%
#   select(Percentage_reloc, updated_yet_perc, Protein)%>%
#   filter(Protein %like% 'XRS2')
# 
# copies_yet_pen_comparison <- ggdotplotstats(data = df_reps,
#                                           x = updated_yet_perc,
#                                           y = Protein,
#                                           # point.args = list(colour = npg_clrs[1]),
#                                           centrality.type = 'np',
#                                           title= "Confirmation of Timecourse  Penetrance consistency",
#                                           x_lab = "Percentage relocalization by Tracked Cell Series",
#                                           test.value = min(df_reps$updated_yet_perc)
# )+
#   scale_color_npg()
# copies_yet_pen_comparison

df_summ <- df %>% 
  filter(Frames_post_treatment == 'Start') %>% 
  summarise(TC_m = median(Percentage_reloc_less, na.rm = T),
            TC_s = skewness(Percentage_reloc_less, na.rm = T),
            TC_var = var(Percentage_reloc_less, na.rm = T),
            TC_sd = sd(Percentage_reloc_less, na.rm = T),
            tc_iqr = IQR(Percentage_reloc_less, na.rm = T))
  
  