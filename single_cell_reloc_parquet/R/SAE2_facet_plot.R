#* This script is for creating a facet plot comparison of abundance and localization.
#* Included below is a hard coded script for SAE2 along with variable for other proteins
library(ggplot2)
library(ggExtra)
library(arrow)
library(ggpointdensity)
library(dplyr)
library(ggstatsplot)
library(corrplot)
library(smplot2)
library(corrplot)

# setwd("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc_parquet/R")
setwd("D:/Second Mind/Academic/Project Stuff/Figures")

#* The below is for the SAE2 plot
df = read_parquet("D:/ALL_FINAL/Combined_by_perc/MODSAE2_selected.parquet")
#df_filter = filter(df, Frames_post_treatment == 20)
numeric_columns <- df %>% select_if(is.numeric)
numeric_columns <- df %>% select(z_score_logLoc, z_score_logAbund, z_score_Loc, z_score_Abund)


# Simple Facet 

gg_temp <-ggplot(df, aes(z_score_logLoc, z_score_logAbund)) +
  #geom_pointdensity() +
  geom_point() +
  sm_statCorr() + 
  #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
  # geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
  # geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  facet_wrap(vars(Frames_post_treatment))+
  ggtitle("Corr_facet_logAbund-logLoc_Pearson_facet")
ggsave("Corr_facet_logAbund-logLoc_Pearson_facet.png", gg_temp, width = 30, height = 18)

gg_temp <-ggplot(df, aes(z_score_Loc, z_score_logAbund)) +
  #geom_pointdensity() +
  geom_point() +
  sm_statCorr() + 
  #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
  # geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
  # geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  facet_wrap(vars(Frames_post_treatment))+
  ggtitle("Corr_facet_logAbund-regLoc_Pearson_facet")
ggsave("Corr_facet_logAbund-regLoc_Pearson_facet.png", gg_temp, width = 30, height = 18)

gg_temp <-ggplot(df, aes(z_score_logLoc, z_score_Abund)) +
  #geom_pointdensity() +
  geom_point() +
  sm_statCorr() + 
  #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
  # geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
  # geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  facet_wrap(vars(Frames_post_treatment))+
  ggtitle("Corr_facet_regAbund-logLoc_Pearson_facet")
ggsave("Corr_facet_regAbund-logLoc_Pearson_facet.png", gg_temp, width = 30, height = 18)

gg_temp <-ggplot(df, aes(z_score_Loc, z_score_Abund)) +
  #geom_pointdensity() +
  geom_point() +
  sm_statCorr() + 
  #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
  # geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
  # geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  facet_wrap(vars(Frames_post_treatment))+
  ggtitle("Corr_facet_regAbund-regLoc_Pearson_facet")
ggsave("Corr_facet_regAbund-regLoc_Pearson_facet.png", gg_temp, width = 30, height = 18)



# M = cor(numeric_columns)
# corrplot(M, method = 'number')
# corrplot.mixed(M, order = 'alphabet')
# corrplot(M, addCoef.col = 'black', tl.pos = 'd',
#          cl.pos = 'n', col = COL2('PiYG'))


list_frames <- unique(df[["Frame"]])

for (c in list_frames){
  df_frame <- filter(df, Frame == c)
  numeric_columns <- df_frame %>% select(z_score_logLoc, z_score_logAbund, z_score_Loc, z_score_Abund)
  M = cor(numeric_columns)
  png(file=sprintf("Corr_frame%s_SAE2.png", c), width=600, height=350)
  corr_graph <- corrplot(M, order = 'AOE', addCoef.col = 'black', type = 'lower', diag = FALSE)
  # ggsave(sprintf("Corr_frame%s_SAE2.png", c), plot = corr_graph)
  dev.off()
}


#The below if for making a facet plot for all the proteins
proteins = list("ZPR1", "ZIP2", "YTA8", "YPR174C", "YOX1", "YOF1", "YMR291W", "YMR160W", "YMR061C", "YML108W", "YML011C", "YLR363W-A", "YLR297W", "YLR126C", "YLR108C", "YKU80", "YKU70", "YKL060W", "YJR056C", "YIL108W", "YGR151C", "YGR042W", "YER064C", "YDR348C", "YDR170W-A", "YDR132C", "YDR115W", "YDR089W", "YDL156W", "YDL129W", "YDL111C", "YBR259W", "YBR197C", "YAR009C", "YAP1", "XRS2d0215", "XRS2d0210", "XBP1", "VPS1", "VPH1", "ULS1", "ULP2", "ULP1", "UFD4", "UBC9", "TUB1", "TSR3", "TSR1", "TSC13", "TSA1", "TRM112", "TOS4", "TOP3", "TOF2", "TIS11", "TDR3", "SVL3", "SUT1", "SUB2", "STB4", "STB2", "SRS2", "SRP68", "SQS1", "SPT21", "SNT2", "SLX8", "SLX4", "SLD3", "SLD2", "SKG3", "SIP5", "SGT2", "SGS1", "SFH5", "SEC3", "SEC11", "SCH9", "SCD6", "SAE2", "SAC6", "RTR2d0222r1", "RTR2d0215r2", "RSF2", "RRP5", "RRP17d0217", "RRP17d0210", "RRB1", "RQC2", "RPN4", "RPL40A", "RPL15B", "RPC10", "RNR4", "RNR1d0222r1", "RNR1d0216r1", "RMT2", "RMI1", "RME1", "RIM1", "RGA1", "RFC4", "RFC3", "RFA2", "RFA1d0213", "RFA1d0210", "REV1", "RDH54", "RBD2", "RAS1", "RAD9", "RAD5", "RAD5d0223r1", "RAD57", "RAD55", "RAD54", "RAD53", "RAD51", "RAD24", "QCR6", "PXL1", "PSY1", "PSO2", "PRS5", "PRE3", "PPN1", "PPH3", "PPH22d0214", "PPH22d0210", "PPH21", "POL30", "PNC1", "PKP2", "PHO81", "PEX29", "PEX21", "PDR3", "PBP4", "PBP2", "PBP1", "PAT1", "OPY2", "NSG1", "NOP58", "NOP13", "NMD4", "NEJ1", "NAM7", "MTR10", "MSN2d0222r2", "MSN2d0218r2", "MSH3", "MSD1", "MSB3", "MSB1", "MRT4", "MRS8", "MRE11", "MRC1", "MODSAE2", "MOB1", "MMS21", "MKT1", "MGS1", "MCM2", "LST8", "LSM7", "LSM4", "LSM3d0217", "LSM3d0214r1", "LSM3d0210r1", "LSM2", "LSM1", "LSM12", "LRS4", "LOC1", "LCD1", "LAP4", "LAG1", "KTR3", "KTR1", "IZH4", "IWR1", "ITR1", "IQG1", "IPL1", "INO80", "HTA2", "HSP42", "HSP26", "HSL7", "HOS2", "HNT3", "HMG2", "HMG1", "HGH1", "HAC1", "HAA1", "GYP5", "GYL1", "GTB1", "GSY2", "GSY1", "GLN1", "GLC3", "GIC1", "GCD8", "FUI1", "FPR2", "FLR1", "FIG4", "FGV2", "FAR1", "FAA4", "FAA1", "EXO70", "EXO1", "ESP2", "ENT1", "EDC3", "EDC2", "EDC1", "ECO1", "ECM3", "ECM29", "DUS3", "DSF2", "DSE3", "DSE1", "DPB11", "DOT6", "DOA1", "DMA2", "DIP5", "DHH1", "DDC2d0224r1", "DDC2d0223r1", "DDC2d0222r2", "DDC2d0222r1", "DDC2d0221r1", "DDC2d0220r2", "DDC2d0220r1", "DDC2d0219r1", "DDC2d0218r2p90KO", "DDC2d0218r2p80KO", "DDC2d0218r2p70KO", "DDC2d0218r2p60KO", "DDC2d0218r2p50KO", "DDC2d0218r2p40KO", "DDC2d0218r2p30KO", "DDC2d0218r2p20KO", "DDC2d0218r2p110KO", "DDC2d0218r2p100KO", "DDC2d0218r1", "DDC2d0217r2", "DDC2d0217r1", "DDC2d0216r2", "DDC2d0216r1", "DDC2d0215r2", "DDC2d0215r1", "DDC2d0214r2", "DDC2d0214r1", "DDC2d013r1", "DDC1", "DCP2", "DCP1", "DBF4", "CYK3", "CTR86", "CTR1", "CSM1", "CRM11", "Contam.", "CMS1", "CLB3", "CHS7", "CHK1", "CGR1", "CDC6", "CDC48", "CDC40", "CDC27", "CDC24", "CDC20", "CDC15", "CDC14", "CCC1", "CBP2", "CAT8", "CAP2", "CAP1", "CAF8", "BUD4", "BTN2", "BMH2", "BMH1", "BIR1", "BCH1", "AVL9", "ATG3", "ATG29", "ATG18", "ATC1", "ASE1", "ARP4", "ARO4", "APJ1", "APC5", "APC4", "APC1", "AME1", "AIP1", "AIM20", "AGR2", "AFT1", "ACE2", "AAC1")
for (prot in proteins){
  file = sprintf('D:/ALL_FINAL/Combined_by_perc/%s_selected.parquet', prot[1])
  df <- read_parquet(file)
  gg_temp <-ggplot(df, aes(z_score_Loc, z_score_logAbund)) +
    #geom_pointdensity() +
    geom_point() +
    sm_statCorr() + 
    #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
    # geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
    # geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
    facet_wrap(vars(Frames_post_treatment))+
    ggtitle(sprintf("%s_Corr_facet_logAbund-regLoc_Pearson_facet.png", prot[1]))+
    theme(plot.title = element_text(size=14))
  ggsave(sprintf("%s_Corr_facet_logAbund-regLoc_Pearson_facet.png", prot[1]), gg_temp, width = 30, height = 18)
}
       
