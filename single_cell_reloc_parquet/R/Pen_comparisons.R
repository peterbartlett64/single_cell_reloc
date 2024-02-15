# HEADER --------------------------------------------
#
# Author:     Peter Bartlett
# Copyright     Copyright 2024 - Peter Bartlett
# Email:      p.bartlett@mail.utoronto.ca
#
# Instance:     2024-02-01
#
# Script Name:    
#
# Script Description:
#
#
# SETUP ------------------------------------
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(ggsci)
library(dplyr)
library(hrbrthemes)
library(wesanderson)
library(ggsci)
library(gridExtra)
library(scales)
library(stringr)
library(cowplot)
library(ggstatsplot)
library(tidylog)

# Module Code--------------------------------------------
date <- Sys.Date()
setwd("D:/Second Mind/Academic/Project Stuff/Figures")




ggplot(filter(df, Protein %in% Protein %in% c("RFA1", "HTA2", "RDH54", "POL30",
                                              "NPL4", "CHK1", "MRC1", "MSH3", "MMS21",
                                              "RAD51", "RAD24", "RTT107", "RRD1", "SLD2",
                                              "DOA1", "ECO1", "RPN4", "SUB2", "RAD57",
                                              "DBF4", "PPH3", "RAD55", "CDC1", "RAD9",
                                              "XRS2", "LRS4", "SLD3", "INO80", "RAD54",
                                              "SAE2", "ZIP2", "ARP4", "DPB11", "SRS2",
                                              "CDC6", "RFC2", "SLX4", "TOP3", "NEJ1",
                                              "RAD33", "RAD52", "NAM7", "YKU80", "PSO2",
                                              "SGS1", "MRE11", "YKU70", "MGS1", "RAD50",
                                              "RFC3", "RFC4", "RTS1", "EXO1", "ULS1",
                                              "REV1", "RAD53", "DDC1", "RIM1")),
       aes(x = Frames_post_treatment, y = Loc_score, group = 'Protein'))+
  geom_line(aes(x = Frames_post_treatment, ))

