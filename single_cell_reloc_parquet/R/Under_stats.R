pen_t <- pen_t %>% 
  mutate(Percentage_relocRn = percent_rank(Percentage_reloc),
         updated_yet_percRn = percent_rank(updated_yet_perc))


df_testing <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>% 
  filter(Protein == 'FLR1', Frames_post_treatment %in% c(0,8,16,24,32)) %>% 
  select("CurrNot", "CurrYes", "currProportion", "CurrNotYet", "CurrYet", "currYetProportion", "YetVelocity", "YetAcceleration") %>% 
  collect() %>% 
  head(5)

df_testing[,-1] <-round(df_testing[,-1],digits = 2)

under <- data.frame(cbind(names(df_testing), t(df_testing))) %>% 
  gt()
gtsave(under, filename = "FLR1_udnerstats.pdf")
gtsave(under, filename = "FLR1_udnerstats.html")



df_testing <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_everything.parquet", as_data_frame = F) %>% 
  filter(Protein == 'EXO1', Frames_post_treatment %in% c(0,8,16,24,32)) %>% 
  select("CurrNot", "CurrYes", "currProportion", "CurrNotYet", "CurrYet", "currYetProportion", "YetVelocity", "YetAcceleration") %>% 
  collect() %>% 
  head(5)

df_testing[,-1] <-round(df_testing[,-1],digits = 2)

under <- data.frame(cbind(names(df_testing), t(df_testing))) %>% 
  gt() %>% 
  gtsave(filename = "EXO1_udnerstats.pdf")
