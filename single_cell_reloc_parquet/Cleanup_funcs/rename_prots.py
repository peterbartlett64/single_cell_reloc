#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
rename_dict = pd.read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/RenameTable.xlsx", usecols=['Reference Array (Universe)', 'Table8.Map_COPU']).dropna().set_index('Table8.Map_COPU')
# %%
test = pd.Series(rename_dict.index.value_counts())
test[test > 1]
#%%
#, Take the index (old name) and convert as the key in the dictionary and the value as the value (new name)
dictionary = rename_dict.drop_duplicates().loc[~rename_dict.index.duplicated()].squeeze().to_dict() #* It is safer to just drop duplicates in the rename dictionary, rather than trying to merge them or keep the first or last one.


#%%
# Percs = pd.read_csv("C:/Users/pcnba/OneDrive - University of Toronto/Desktop/Percs_rep.csv") #. This is the old one
Percs = pd.read_parquet("D:/ALL_FINAL/Combined_by_perc/penetrance_updated_trimmed.parquet")
Percs['Unmod_Name'] = Percs['Protein']
Percs['Protein'] = Percs['Protein'].map(dictionary).fillna(Percs['Unmod_Name'])
Percs.set_index('Protein', inplace=True)
Percs.to_csv("C:/Users/pcnba/OneDrive - University of Toronto/Desktop/Percs_rep_MODNAMES.csv")
#%%
Schemas = pd.read_csv("D:/Second Mind/Academic/Project Stuff/Figures/sankey_data_move_TopSplit.csv")
Schemas = Schemas.loc[Schemas['Frames_post_treatment'] == 'Start']
Schemas['Unmod_Name'] = Schemas['Protein']
Schemas['Protein'] = Schemas['Protein'].map(dictionary).fillna(Schemas['Unmod_Name'])
Schemas.to_csv("D:/Second Mind/Academic/Project Stuff/Figures/sankey_data_move_TopSplit_MODNAMES.csv")
#%%

#, THis is the recursive code to replace all protien names for the library based on the map. HIGHLY recommend making a copy of the files before running

import glob
files = glob.glob("D:/ALL_FINAL/Combined_by_perc/*") #. Find all dataframe files of either parquet or csv format
files = pd.DataFrame(files)
files = files.iloc[:,0].str.split(".", expand=True).drop(columns=2)
files.rename(columns={0:'path', 1:'ext'}, inplace=True)

for f in files.itertuples():
	if f[2] == 'parquet':
		temp_f = pd.read_parquet(f[1] + '.' + f[2])
	elif f[2] == 'csv':
		temp_f = pd.read_csv(f[1] + '.' + f[2])
	else:
		continue

	#. This maps the new names to the old names and fills in the old names if the new names are not found
	temp_f['Unmod_Name'] = temp_f['Protein']
	temp_f['Protein'] = temp_f['Protein'].map(dictionary).fillna(temp_f['Unmod_Name'])

	if f[2] == 'parquet':
		temp_f.to_parquet(f[1] + '.parquet')
	elif f[2] == 'csv':
		temp_f.to_csv(f[1] + '.csv')
