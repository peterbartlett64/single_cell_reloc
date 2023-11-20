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

setwd("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc_parquet/R")

df = read_parquet("D:/ALL_FINAL/Combined_by_perc/SAE2_selected.parquet")
#df_filter = filter(df, Frames_post_treatment == 20)

# Simple Facet 
ggplot(df, aes(z_score_logLoc, z_score_logAbund)) +
  #geom_pointdensity() +
  geom_point() +
  sm_statCorr() + 
  #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
  geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
  facet_wrap(vars(Frames_post_treatment))



#The below if for making a facet plot for all the proteins
proteins = list("AAC1", "ACE2" , "AFT1" , "AGR2" , "AIM20" , "AIP1" , "AME1" , "APC1" , "APC4" , "APC5" , "APJ1" , "ARO4" , "ARP4" , "ASE1" , "ATC1" , "ATG3" , "ATG18" , "ATG29" , "AVL9" , "BCH1" , "BIR1" , "BMH1" , "BMH2" , "BTN2" , "BUD4" , "CAF8" , "CAP1" , "CAP2" , "CAT8" , "CBP2" , "CCC1" , "CDC6" , "CDC14" , "CDC15" , "CDC20" , "CDC24" , "CDC27" , "CDC40" , "CDC48" , "CGR1" , "CHK1" , "CHS7" , "CLB3" , "CMS1" , "Contam." , "CRM11" , "CSM1" , "CTR1" , "CTR86" , "CYK3" , "DBF4" , "DCP1" , "DCP2" , "DDC1" , "DDC2d013r1" , "DDC2d0214r1" , "DDC2d0214r2" , "DDC2d0215r1" , "DDC2d0215r2" , "DDC2d0216r1" , "DDC2d0216r2" , "DDC2d0217r1" , "DDC2d0217r2" , "DDC2d0218r1" , "DDC2d0218r2p20KO" , "DDC2d0218r2p30KO" , "DDC2d0218r2p40KO" , "DDC2d0218r2p50KO" , "DDC2d0218r2p60KO" , "DDC2d0218r2p70KO" , "DDC2d0218r2p80KO" , "DDC2d0218r2p90KO" , "DDC2d0218r2p100KO" , "DDC2d0218r2p110KO" , "DDC2d0219r1" , "DDC2d0220r1" , "DDC2d0220r2" , "DDC2d0221r1" , "DDC2d0222r1" , "DDC2d0222r2" , "DDC2d0223r1" , "DDC2d0224r1" , "DHH1" , "DIP5" , "DMA2" , "DOA1" , "DOT6" , "DPB11" , "DSE1" , "DSE3" , "DSF2" , "DUS3" , "ECM3" , "ECM29" , "ECO1" , "EDC1" , "EDC2" , "EDC3" , "ENT1" , "ESP2" , "EXO1" , "EXO70" , "FAA1" , "FAA4" , "FAR1" , "FGV2" , "FIG4" , "FLR1" , "FPR2" , "FUI1" , "GCD8" , "GIC1" , "GLC3" , "GLN1" , "GSY1" , "GSY2" , "GTB1" , "GYL1" , "GYP5" , "HAA1" , "HAC1" , "HGH1" , "HMG1" , "HMG2" , "HNT3" , "HOS2" , "HSL7" , "HSP26" , "HSP42" , "HTA2" , "INO80" , "IPL1" , "IQG1" , "ITR1" , "IWR1" , "IZH4" , "KTR1" , "KTR3" , "LAG1" , "LAP4" , "LCD1" , "LOC1" , "LRS4" , "LSM1" , "LSM2" , "LSM3d0210r1" , "LSM3d0214r1" , "LSM3d0217" , "LSM4" , "LSM7" , "LSM12" , "LST8" , "MCM2" , "MGS1" , "MKT1" , "MMS21" , "MOB1" , "MRC1" , "MRE11" , "MRS8" , "MRT4" , "MSB1" , "MSB3" , "MSD1" , "MSH3" , "MSN2d0218r2" , "MSN2d0222r2" , "MTR10" , "NAM7" , "NEJ1" , "NMD4" , "NOP13" , "NOP58" , "NSG1" , "OPY2" , "output_selection.csv" , "parquet_test.R" , "PAT1" , "PBP1" , "PBP2" , "PBP4" , "PDR3" , "PEX21" , "PEX29" , "PHO81" , "PKP2" , "PNC1" , "POL30" , "PPH3" , "PPH21" , "PPH22d0210" , "PPH22d0214" , "PPN1" , "PRE3" , "PRS5" , "PSO2" , "PSY1" , "PXL1" , "QCR6" , "RAD5" , "RAD5d0223r1" , "RAD9" , "RAD24" , "RAD51" , "RAD53" , "RAD54" , "RAD55" , "RAD57" , "RAS1" , "RBD2" , "RDH54" , "REV1" , "RFA1d0210" , "RFA1d0213" , "RFA2" , "RFC3" , "RFC4" , "RGA1" , "RIM1" , "RME1" , "RMI1" , "RMT2" , "RNR1d0216r1" , "RNR1d0222r1" , "RNR4" , "RPC10" , "RPL15B" , "RPL40A" , "RPN4" , "RQC2" , "RRB1" , "RRP5" , "RRP17d0210" , "RRP17d0217" , "RSF2" , "RTR2d0215r2" , "RTR2d0222r1" , "SAC6" , "SAE2" , "SCD6" , "SCH9" , "SEC3" , "SEC11" , "SFH5" , "SGS1" , "SGT2" , "SIP5" , "SKG3" , "SLD2" , "SLD3" , "SLX4" , "SLX8" , "SNT2" , "SPT21" , "SQS1" , "SRP68" , "SRS2" , "STB2" , "STB4" , "SUB2" , "SUT1" , "SVL3" , "TDR3" , "TIS11" , "TOF2" , "TOP3" , "TOS4" , "TRM112" , "TSA1" , "TSC13" , "TSR1" , "TSR3" , "TUB1" , "UBC9" , "UFD4" , "ULP1" , "ULP2" , "ULS1" , "VPH1" , "VPS1" , "XBP1" , "XRS2d0210" , "XRS2d0215" , "YAP1" , "YAR009C" , "YBR197C" , "YBR259W" , "YDL111C" , "YDL129W" , "YDL156W" , "YDR089W" , "YDR115W" , "YDR132C" , "YDR170W-A" , "YDR348C" , "YER064C" , "YGR042W" , "YGR151C" , "YIL108W" , "YJR056C" , "YKL060W" , "YKU70" , "YKU80" , "YLR108C" , "YLR126C" , "YLR297W" , "YLR363W-A" , "YML011C" , "YML108W" , "YMR061C" , "YMR160W" , "YMR291W" , "YOF1" , "YOX1" , "YPR174C" , "YTA8" , "ZIP2" , "ZPR1")

for (prot in proteins){
  file = sprintf('D:/ALL_FINAL/Combined_by_perc/%s_selected.parquet', prot[1])
  df <- read_parquet(file)
  gg <- ggplot(df, aes(z_score_logLoc, z_score_logAbund)) +
    geom_pointdensity() +
    #geom_point() +
    sm_statCorr() + 
    #geom_pointdensity(x = 'z_score_logLoc', y = 'z_score_logAbund')+
    #geom_xsidedensity(aes(y = after_stat(density)), position = "stack") +
    #geom_ysidedensity(aes(x = after_stat(density)), position = "stack") +
    facet_wrap(vars(Frames_post_treatment))
  prot = "SAE2"
  ggsave(sprintf("Corr_facet%s_Pearson_facet.png", prot[1]), gg, width = 30, height = 18)
}
       
