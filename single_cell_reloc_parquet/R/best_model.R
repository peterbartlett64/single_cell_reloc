Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre-1.8')
library(car)        # extracts model results
library(MASS)       # provides "birthwt" dataset
library(ISLR)       # provides "Wage" dataset
library(tictoc)     # checks running time
library(sjPlot)     # visualizes model results
library(glmulti)    # finds the BEST model
library(flextable)  # beautifies tables
library(tidyverse)  # provides a lot of useful stuff !!! 
library(performance)# checks and compares quality of models
library(arrow)
library(strengejacke)
library(effects)


theme_set(theme_light(base_size = 12)) # beautifies plots
theme_update(panel.grid.minor = element_blank())
# prepare selection

df_all <- read_parquet("D:/ALL_FINAL/Combined_by_perc/merged_data_final.parquet")

df_all$Cell_Barcode <- as.factor(df_all$Cell_Barcode)
df_all$Date <- as.factor(df_all$Date)
df_all$Unique_Frame <- as.factor(df_all$Unique_Frame)
df_all$Is_treated <- as.factor(df_all$Is_treated)
df_all$ImageID <- as.factor(df_all$ImageID)
df_all$Unique_pos <- as.factor(df_all$Unique_pos)
df_all$Protein <- as.factor(df_all$Protein)

df_all$Relocalized <- as.logical(df_all$Relocalized)
df_all$Progen_bud <- as.logical(df_all$Progen_bud)
df_all$Yet <- as.logical(df_all$Yet)
df_all$Yes_yet <- as.logical(df_all$Yes_yet)
df_all$No_yet <- as.logical(df_all$No_yet)
df_all$Pos <- as.factor(df_all$Pos)

suteb <- filter(df_all, Protein == "RAD51") %>% as.data.frame()
suteb

nested_df <- df_all %>% 
    group_by(Protein, Unique_pos) %>%
    nest() #* This will group the data so that is separated and is not called every time a graph is created. This is in testing
nested_models <- nested_df %>% 
  mutate(models = map(data, ~ glmulti(Loc_score ~ log_Abundance + CoV_spos + cell_area + Is_treated,
                                        data   = ., 
                                        level  = 1,          # 2 with interactions, 1 without  
                                        method = "h",        # "d", or "h", or "g"
                                        fitfunction = glm,   # Type of model (LM, GLM etc.)
                                        confsetsize = 10)))





tic()
h_test <- glmulti(Loc_score ~ log_Abundance + CoV_spos + cell_area + Is_treated,
        data   = df_all, 
        level  = 1,          # 2 with interactions, 1 without  
        method = "h",        # "d", or "h", or "g"
        fitfunction = glm,   # Type of model (LM, GLM etc.)
        confsetsize = 10)   # Keep N best models
toc()

compare_performance(h_test@objects[[1]], h_test@objects[[2]])

plot(effects::allEffects(h_test@objects[[1]]),
     lines = list(multiline = T),
     confint = list(style = "auto"))

