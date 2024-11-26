#Cleanup_merged

#, This script is just to confirm which proteins are annotated to move in MMS
#%%
if __name__ == '__main__':
	import pandas as pd
	import os
	Global_Variables = {
	# 	"analyze": "F:/Microfluidics/Missing_Analyze2",
		"microfluidics_results": "F:/Microfluidics/RES_N_ULTS",
		"information_path": "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging",
		"post_path": "D:/ALL_FINAL"}
	os.chdir(Global_Variables['information_path'])

	movers_mms = pd.read_csv('library_prots_moveMMS.csv')
	proteins = pd.read_csv("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc_parquet/R/new_temp.csv")

	matched = pd.merge(proteins, movers_mms, left_on='Protein',right_on= 'Standard Name', how='outer')
	