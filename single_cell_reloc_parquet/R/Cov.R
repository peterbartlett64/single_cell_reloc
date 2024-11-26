#This is the plot for creating facet plots of the abundance versus the localization.
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)

ggplot(aes(x=Log_loc, y=Log_abund)
       + geom_point )

ggplot(df, aes(z_score_logLoc, z_score_logAbund)) +
  #geom_pointdensity() +
  geom_point() +
  sm_statCorr() + 
  #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
  geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  geom_vline(xintercept=0, size=1.5, color="red") + 
  geom_vline(xintercept=median(x), size=1.5, color="blue", linetype = 3) + 
  geom_hline(yintercept=0, size=1.5, color="red") +
  geom_vline(xintercept=median(y), size=1.5, color="blue", linetype = 3) + 
  facet_wrap(vars(Frames_post_treatment))

for (prots in )