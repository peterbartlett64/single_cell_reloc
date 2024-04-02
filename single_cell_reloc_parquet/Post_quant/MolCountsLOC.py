
#, This file is workable merged the ReLoc information (LocScore) with the true molecular count
def read_mol_counts(file):
	df = pd.read_excel(file, sheet_name = 'Table S8')
	all_cols = df.columns

	#Get the columns for the untreated
	untreated_cols = [col for col in all_cols if 'Untreated' in col]

	# list_untreated = ["BRE_Untreated","CHO_Untreated","DAV_Untreated" "DEN_Untreated", "DEN_MMS" "MAZ_Untreated", "MAZ_MMS",	"LEE_Untreated", "LEE_MMS",	"TKA_Untreated", "TKA_MMS"]

	df['Untreated_mean'] = df[untreated_cols].mean(axis=1)
	df['Untreated_std'] = df[untreated_cols].std(axis=1)

	return(df)

def read_loc_file(file):
	df = pd.read_parquet(file, columns=['Gene_Standard_Name', 'Loc_score', 'Cell_Barcode', 'Unique_Frame', 'Frames_post_treatment', 'Abundance', 'log_Abundance', 'z_score_Abund'])
	return(df)

def main():
	mol_counts = read_mol_counts(file = "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/TableS8_Final.xlsx")# gv.slash_switch(input('What is the mol counts file?')))
	loc_file = read_loc_file(file = 'D:\ALL_FINAL\Combined_by_perc\Loc_data_comp_merged_Tkach.parquet') #gv.slash_switch(input('What is the loc file?')))

	merged = pd.merge(loc_file, mol_counts, left_on= 'Gene_Standard_Name', right_on= 'Standard Name', how = 'left')
	return(merged)



if __name__ == '__main__':
	import pandas as pd
	import arrow
	import pyarrow
	import pyarrow as pa
	import pyarrow.parquet as pq
	import os
	import single_cell_reloc_parquet.global_functions.global_variables as gv

	#Todo: finish the variable read here. Because the name of the loc file is variable and the mol counts file is variable, decided to just do a simple import
	# Global_variables = gv.Global_variables()
	# Global_variables['post_path'] = 'D:/ALL_FINAL'
	# os.join(Global_variables['post_path'], 'Combined_by_perc')
	# os.chdir(Global_variables['post_path'])

	table = pa.Table.from_pandas(main())
	pq.write_table(table, 'D:/ALL_FINAL/Combined_by_perc/Loc_data_comp_merged_mol_counts.parquet')