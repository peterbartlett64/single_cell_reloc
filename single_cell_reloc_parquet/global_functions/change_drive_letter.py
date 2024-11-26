#%%
import pandas as pd
import pathlib
import os
from single_cell_reloc_parquet.global_functions.global_variables import slash_switch

def change_letter(file,orig_letter, change_letter):
	# temp_file = file.replace(f"{orig_letter}:", f"{change_letter}:" #? This line was for when reading based on the old path. This is unlikely as the user would have no reason to do so. Left here to be uncommented just in case though

	temp_file = file
	file_ext = pathlib.Path(temp_file).suffix

	if file_ext == '.parquet':
		temp = pd.read_parquet(temp_file)
	elif file_ext == '.csv':
		temp = pd.read_csv(temp_file)


	try:
		temp['Path'] = pd.Series(temp['Path']).apply(lambda x: x.replace(f"{orig_letter}:", f"{change_letter}:"))
	except:
		return(f"Not and index file: {file}")

	if file_ext == '.parquet':
		temp.to_parquet(temp_file)
	elif file_ext =='.csv':
		temp.to_csv(temp_file)

#%%
if __name__ == '__main__':
	# microfluidics_results = slash_switch(input('microfluidics_results')) #* This could have been written without the slash switching, but it should be more stable between systems
	# df_list = input("Files to be mod")
	# df_list = df_list.split(', ')
	df_list = ["Allmasks.csv", "Allmasks.parquet", "info_index.csv", "info_index.parquet", "Allmasks_exp.csv", "Allmasks_exp.parquet", "orgAllmasks.csv", "orgAllmasks.parquet", "imgIndex.csv", "imgIndex.parquet"]
	change_f = "D" #input("change from")
	change_t = "E" #input("change to")
	for d in df_list:
		df = os.path.join(microfluidics_results, d)
		change_letter(file= df, orig_letter= change_f , change_letter= change_t)
		print(d)

# change_letter("D:\\Microfluidics\\RESULTS\\d0214r1p010200\\GFP_mKO_mKa\\2023613220212_TrackingResults_d0214r1p010200_2023-06-13.txt", "D", "E")