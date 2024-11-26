#%%
import pandas as pd
import os
import statistics as stats
import numpy as np
import math
import numpy as np
#Todo:convert this file into function form so that it is callable
#%%
#Todo: complete functions
def relate_compare(cell_value, pop_value, mol_value):
	ratio = mol_value/pop_value
	#* There may need to be a log tranformation
	cell_mol = cell_value * ratio

def convert_time_frames(df, hrs, frame_gap = 7.5):
	frames_post = (hrs*60)/(frame_gap)
	return(df.loc['Frames_post_treatment'] == frames_post)

def convert_frames_time(x, frame_gap = 7.5):
	time_post_post = (x*frame_gap) #. Decide if this should be in hours or minutes
	return(time_post_post)
#%%
# def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
#     new = path.replace(os.sep, '/')
#     return (new)

# Directory = input("directory?")
# Directory = slash_switch(Directory)
# os.chdir(Directory)

#%%
# File_n = input("file_name")
# # File = pd.read_csv(f"{File_n}.csv")
# File = pd.read_parquet(f"{File_n}.parquet")
# print("File has been read in")

# File.to_parquet(f"{File_n}.parquet")
# print("File has been saved to paraquet")
#%%

# File_n2 = input("file_name2")
# File2 = pd.read_csv(File_n2)

# l = int(input("files"))
# for n in range(l):
#     File = input("file name?")
#     file = pd.read_csv(File + ".csv", sep = " ")
#     file_np = file.copy()
#     file = pd.merge(file_n2, file, left_index=True, right_index= True)


#%%
# Abundance_information = pd.read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Abundance_Mod_Stress.xlsx", sheet_name="Table S8")
#? UserWarning: Unknown extension is not supported and will be removed for idx, row in parser.parse():
Abundance_information = pd.read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/ProteinAbundance-main/supplementalTables/TableS8_Final_Mod.xlsx", engine='openpyxl')
Localization_information = pd.read_parquet("D:/Microfluidics/Most_final_collected/Most_final_collected/Combined_by_perc/Col_with_Abund/Merged/Final_wAbundLocal.parquet")

Abundance_information.drop(Abundance_information.columns[Abundance_information.columns.str.contains('Unnamed',case = False)],axis = 1, inplace= True)
Localization_information.drop(Localization_information.columns[Localization_information.columns.str.contains('Unnamed',case = False)],axis = 1, inplace= True)


#%% Get the population abundances
#* Group cells by the strain frame
Localization_information['Strain_frame'] = Localization_information['Protein'] + Localization_information['ImageID'].apply(lambda x: x[-5:])
#%%
try: Localization_information['Time_post_treatment']
except KeyError:
	Localization_information['Time_post_treatment'] = pd.Series(Localization_information['Frames_post_treatment']).apply(lambda x: convert_frames_time(x, frame_gap = 7.5))
#%% Group by the protein time, so that can aggregate mean and get a global measure of cell abundance

#* Step 1: Get rid of the extra columns and convert from A.U to mols per cell
Localization_information_cpy = Localization_information.loc[:,['Protein', 'Frame_x', 'Time_post_treatment','Cell_Barcode', 'ImageID', 'Abundance', 'Strain_frame']].copy()
Localization_information_cpy['Mol_Abundance'] = pd.Series(Localization_information_cpy['Abundance']).apply(lambda x: math.log(x))

#* Step 2: Group by protein time, create independent dataframes and get the means (to deal with different cell counts)
global_pop_Abund =  Localization_information.loc[:,['Protein', 'Frame_x', 'Time_post_treatment','Cell_Barcode', 'ImageID', 'Abundance', 'Strain_frame']].copy()
global_pop_Abund['Median_Mol_Abundance'] = Localization_information_cpy.groupby('Strain_frame')['Mol_Abundance'].transform('median')

# zero_hrs = Localization_information_cpy.loc[Localization_information_cpy['Time_post_treatment'] == 0].reset_index(drop=True)
# zero_hrs_median = zero_hrs.groupby('Protein')['Mol_Abundance'].agg('median')
# two_hrs = Localization_information_cpy.loc[Localization_information_cpy['Time_post_treatment'] == 120].reset_index(drop=True)
# two_hrs_median = two_hrs.groupby('Protein')['Mol_Abundance'].agg('median')
# four_hrs = Localization_information_cpy.loc[Localization_information_cpy['Time_post_treatment'] == 240].reset_index(drop=True)
# four_hrs_median =four_hrs.groupby('Protein')['Mol_Abundance'].agg('median')

#* Step 3: Get modes for each time dataframe
zero_hrs_mode = stats.mode(global_pop_Abund.loc[global_pop_Abund['Time_post_treatment'] == 0].drop(columns=['Cell_Barcode', 'Abundance']).drop_duplicates())
two_hrs_mode = stats.mode(global_pop_Abund.loc[global_pop_Abund['Time_post_treatment'] == 120].drop(columns=['Cell_Barcode', 'Abundance']).drop_duplicates())
four_hrs_mode = stats.mode(global_pop_Abund.loc[global_pop_Abund['Time_post_treatment'] == 240].drop(columns=['Cell_Barcode', 'Abundance']).drop_duplicates())

global_pop_Abund.groupby('Time_post_treatment')['Mean_Mol_Abundance'].agg(pd.Series.mode)

zero_hrs_mode = stats.mode(zero_hrs_median)
two_hrs_mode = stats.mode(two_hrs_median)
four_hrs_mode = stats.mode(four_hrs_median)

#* Step 4: Relate the original measurements to the time modes
zero_hrs["shifted_Mol_Abundance"] = pd.Series(zero_hrs.loc[:,"Mol_Abundance"]).apply(lambda x: x*(100/zero_hrs_mode))
two_hrs["shifted_Mol_Abundance"] = pd.Series(two_hrs.loc[:,"Mol_Abundance"]).apply(lambda x: x*(100/two_hrs_mode))
four_hrs["shifted_Mol_Abundance"] = pd.Series(four_hrs.loc[:,"Mol_Abundance"]).apply(lambda x: x*(100/two_hrs_mode))

#%%
Abundance_information['unt_ref'] = (Abundance_information['DEN_Untreated'] + Abundance_information['MAZ_Untreated'] + Abundance_information['LEE_Untreated'] + Abundance_information['TKA_Untreated'])/4 #* This was easier than using mean function
Abundance_information['two_hrs_ref'] = (Abundance_information['LEE_MMS'] + Abundance_information['TKA_MMS'])/2
Abundance_information['four_hrs_ref'] = Abundance_information[ 'MAZ_MMS']


ref = Abundance_information.loc[:, ['Standard Name', 'unt_ref', 'two_hrs_ref', 'four_hrs_ref']]
del Abundance_information

#%% Compare the abundance data generated to Brandon's aggregated data
#* Step 1: Merge the data from this study and the MMS study times
zero_hrs = pd.merge(zero_hrs, ref, left_on = 'Protein', right_on = 'Standard Name', how = 'left')
two_hrs = pd.merge(two_hrs, ref, left_on = 'Protein', right_on = 'Standard Name', how = 'left')
four_hrs = pd.merge(four_hrs, ref, left_on = 'Protein', right_on = 'Standard Name', how = 'left')


#%% Relate the time means to reference datasets via least-squares regression
alpha_zero = np.linalg.lstsq(zero_hrs['zero_hrs_ref'], zero_hrs["shifted_Mol_Abundance"])
alpha_two = np.linalg.lstsq(two_hrs['two_hrs_ref'], two_hrs["shifted_Mol_Abundance"])
alpha_four = np.linalg.lstsq(four_hrs['four_hrs_ref'], four_hrs["shifted_Mol_Abundance"])

#%% Merge Data
#%%
# merged = pd.merge(Localization_information, Abundance_information,left_on= "Protein",right_on= "Standard Name", how = 'left')
#%%
# merged.dropna(subset=['MMS_Mean'], inplace= True)
# merged.drop(merged.columns[merged.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)

#%%
merged['Unt_mean'] = merged.groupby('Strain_frame')['Mol_Abundance'].agg('median')
#* 2hours is frame 16 with 7.5 min per
merged['2hrs_mean'] =
#* 4hours is frame 32 with 7.5 min per
merged['4hrs_mean'] =

#%%

#%%

#. This is the hard coded fast version





#Todo: This should be included in more places to reference the stain frames. In many places, there are single cell aggs based on frame groupings.
df['Strain_frame'] = df['ImageID'] + df['Protein']
# %%
