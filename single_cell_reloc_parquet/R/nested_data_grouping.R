Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jre-1.8')
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(corrplot)
library(smplot2)
library(corrplot)
library(tidyverse)
# library(broom)
library(robustbase)
library(sjPlot)
library(olsrr)
library(sjstats)

#* Read in the full version of the file. This will be grouped into a nested df
df <- read_parquet("D:/ALL_FINAL/Combined_by_perc/Final_wAbund.parquet") #This should actually call the global file with all the information

sub <- df %>% 
  filter(Protein == 'SAE2', Frames_post_treatment == 30)
m <- lm(Loc_score ~ Abundance, data = sub)
ols_plot_resid_lev(m) # Run residuals diagnostics, to find obvious outlier
plot_residuals(m)

rm <- lmrob(z_score_Loc ~ z_score_Abund, data = sub)
summary(rm)
# plot_model(rm, type = 'pred')


# nested_df <- df %>% 
#   group_by(Protein, Unique_pos) %>%  
#   nest() #* This will group the data so that is separated and is not called every time a graph is created. This is in testing
# remove(df) #* Free up some space because this is a pretty big dataframe
# 
# nested_models <- nested_df %>% 
#   mutate(models = map(data, ~ lm(z_score_Loc ~ z_score_Abund, data = .)),
#          coeff = map(models, tidy, conf.int = TRUE),
#          quality= map(models, glance),
#          preds = map(models, augment))
# library(performance)
# nested_models$models[[2]] %>% performance::check_model()


# m <- lm(outcome ~ predictor, data = df) #Create linear model

# rm <- lmrob(outcome ~ predictor, data = df) #* gives different residual measurements to each measurements 
                                            #*based on how close to the linear model. Chances the influence of the points
# rm <- lmrob(Loc_score ~ Abundance, data = df)


#Visualize both models
# plot_model(m, type = 'pred', ci_style = 'bca')
plot_model(m, type = 'pred', show.data = T)
plot_model(rm, type = 'pred', show.data = T)

#Compare coefficients and r^2 of both models
tab_model(m)
tab_model(rm)
