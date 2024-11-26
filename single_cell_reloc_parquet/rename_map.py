#%%
import pandas as pd
import os
#%%
rename_dict = {"YBL111C": "YDL111C",
"YKL069W": "YKL060W",
"YML018C": "YML011C",
"YNR061C": "YMR061C",
"ESP1": "NSP1",
"GCD6": "SCD6",
"LSB1": "MSB1",
"RFA1": "HFA1",
"RTN2": "BTN2",
"VTC5": "MTC5",
"XRS2": "SRS2",
"RAD52": "Rad5",
"RAD59": "Rad5",
"AIM44": "AIM4",
"BUD14": "BUD4",
"PPH22": "PPH21",
"APE1": "ASE1",
"ATG2": "ATG29",
"CBF2": "CBP2",
"CMR1": "CMS1",
"CRM1": "CRM11",
"CSM3": "CSM1",
"DNA2": "DMA2",
"ECM5": "ECM3",
"MRS6": "MRS8",
"MSH2": "MSH3",
"MSN2": "MSN1",
"MSO1": "MSB1",
"PAL1": "PIL1",
"PBY1": "PSY1",
"RFC5": "RFC3",
"RGD2": "RBD2",
"RNR1": "RNR4",
"SEC8": "SEC3",
"SIZ1": "Siz2",
"STB3": "STB4",
"STP4": "STB4",
"TDA1": "TSA1",
"TOP1": "TOP3",
"YOP1": "YAP1",
"YTA6": "YTA8",
"YMR111C": "YMR031C"}

percs = pd.read_csv("C:/Users/pcnba/OneDrive - University of Toronto/Desktop/Percs.csv")
percs['Protein_cpy'] = percs['Protein'].map(rename_dict)
percs["Protein_cpy"].fillna(percs["Protein"], inplace=True)
percs.drop("Protein", axis=1, inplace=True)
percs.rename(columns={"Protein_cpy": "Protein"}, inplace=True)
percs.to_csv("C:/Users/pcnba/OneDrive - University of Toronto/Desktop/Percs_rep.csv", index=False)
#%%
pd.read_parquet
#%%
