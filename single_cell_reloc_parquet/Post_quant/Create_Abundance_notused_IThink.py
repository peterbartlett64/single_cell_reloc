#Create Abundance
#%%
from joblib import Parallel, delayed
import pandas as pd
import math
from scipy import stats
from joblib import Parallel, delayed
from datetime import date
import os
import psutil as p
import math
from glob import glob
import single_cell_reloc_parquet.global_functions.global_variables as gv

#%%
def size(byte): #, this the function to convert bytes into more suitable reading format.
  #* I do not know why I have this
  size_GB = byte/(1024)**4
  return(size_GB)

slash_switch = gv.slash_switch


# def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
#     new = path.replace(os.sep, '/')
#     return (new)
#, Capture system information to decide run parameters
mem = p.virtual_memory()
mem_avail = mem.available
RAM_GB = int(mem_avail)/((1024)**3) #* Assume that col file is ~1GB so check how much memory is available
print(RAM_GB)
pn = os.cpu_count()
if RAM_GB >= pn:
	pr = pn
elif RAM_GB <= pn:
	pr = int(RAM_GB)

today = str(date.today())

#%%
if __name__ == '__main__':
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': 'E:/Microfluidics/MOST FINAL', #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': os.cpu_count(),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True}

	post_quant = Global_variables["post_path"]
	percentage_combined_path = os.path.join(post_quant, "Combined_by_perc")
	os.chdir(post_quant)

	All_protein_files = sorted(glob('*.parquet'))
	All_proteins_df = pd.concat((pd.read_parquet(file) for file in All_protein_files), ignore_index= True)
	All_proteins_df.to_parquet("Final_w_Abund.parquet")


#%%
# path = "C:/Users/pcnba/Desktop/Test/Chamber_Col_d0213_r1.0_Ch4.0_2022-05-27.csv"
# col_df = pd.read_csv(path) #! This is functional and probably faster because the f string can be avoided
# #! Remove the date information from the column output for ease of use
# #* Calculate the abundance for every cell at every timepoint. This is done with such density so that associations can be made to the current abundance and the Loc_score
# col_df["Abundance"] = pd.Series(col_df["factor_median_OBJ_GFP"]).apply(lambda x: math.log(x)) #* Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset
# #* purposely overwrit the old version. Because the dataframes are quite large, it is excessive to keep the unmodified df in memory for a couple non-mutative lines
# #? Perfroming with the lambda function might be slightly faster because the function does not need to be declared for each core/run

# group_prot_frame = col_df.groupby(["Protein", "Frame_x"]) #. This may need to be changed to jsut regular "Frame"

# #* This seems like a good plac to grab the z_scores for the Loc_score and Abundance
# col_df["z_score_Loc"] = group_prot_frame["Loc_score"].transform(stats.zscore)
# col_df["z_score_Abund"] = group_prot_frame["Abundance"].transform(stats.zscore)

#%%
#, NEW version of Abundance genration.
#! This was moved into 'Intensity_select_percentage_combine.py'
def Abundance_calc_manager(local_prot_df:pd.DataFrame): #* The abundacne could be incorporated in the post-quant function, but it could slow down. It also makes the process less modular and therefore makes timing of computer location a restraint
	#, Calculate the abundance for every cell at every timepoint. This is done with such density so that associations can be made to the current abundance and the Loc_score
	local_prot_df["Abundance"] = pd.Series(local_prot_df["factor_median_OBJ_GFP"]).apply(lambda x: math.log(x)) #* Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset
	#* purposely overwrit the old version. Because the dataframes are quite large, it is excessive to keep the unmodified df in memory for a couple non-mutative lines
	#? Perfroming with the lambda function might be slightly faster because the function does not need to be declared for each core/run

	group_prot_frame = local_prot_df.groupby(["Protein", "Frame_x"]) #. This may need to be changed to jsut regular "Frame"

	# #* This seems like a good plac to grab the z_scores for the Loc_score and Abundance
	local_prot_df["z_score_Loc"] = group_prot_frame["Loc_score"].transform(stats.zscore)
	local_prot_df["z_score_Abund"] = group_prot_frame["Abundance"].transform(stats.zscore)
	return(local_prot_df)

#%%
# #, OLD version of running Abundance on the column level.
# def old_Abundance_calc_manager(path): #* The abundacne could be incorporated in the post-quant function, but it could slow down. It also makes the process less modular and therefore makes timing of computer location a restraint
# 	try:
# 		col_df = pd.read_csv(path) #! This is functional and probably faster because the f string can be avoided
# 		#! Remove the date information from the column output for ease of use
# 		#* Calculate the abundance for every cell at every timepoint. This is done with such density so that associations can be made to the current abundance and the Loc_score
# 		col_df["Abundance"] = pd.Series(col_df["factor_median_OBJ_GFP"]).apply(lambda x: math.log(x)) #* Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset
# 		#* purposely overwrit the old version. Because the dataframes are quite large, it is excessive to keep the unmodified df in memory for a couple non-mutative lines
# 		#? Perfroming with the lambda function might be slightly faster because the function does not need to be declared for each core/run

# 		group_prot_frame = col_df.groupby(["Protein", "Frame_x"]) #. This may need to be changed to jsut regular "Frame"

# 		# #* This seems like a good plac to grab the z_scores for the Loc_score and Abundance
# 		col_df["z_score_Loc"] = group_prot_frame["Loc_score"].transform(stats.zscore)
# 		col_df["z_score_Abund"] = group_prot_frame["Abundance"].transform(stats.zscore)

# 		col_df.to_csv(path)
# 		return(f"{path} complete")
# 	except:
# 		return(f"FAILED ON {path}")


# try:
# 	path = os.path.join(post_quant, "Col_w_abundance") #* Added this to make the system non-destructive
# 	os.mkdir(path)
# except FileExistsError:
# 	path = os.path.join(post_quant, "Col_w_abundance")
# 	pass

# x = Parallel(n_jobs=pr, verbose= 100)(delayed(Abundance_calc_manager)(path = p) for p in path_cols)
# print(x) # This is a good way of showing the log in the terminal
# x = pd.Series(x) #* This is the log of positions which have been calculated
# x.to_csv(f'Abundance_func_output_{today}.csv')
#%%

#, Calcualte the Frame with maximum Loc_score. This would be far easier to apply to the ALL_pos_df

post_treatment_df = col_df.loc[col_df["Frames_post_treatment"] > 0]
Agg_loc = post_treatment_df.groupby(["Protein", "Frame_x"])["Loc_score"].median() #* First get the median/mean of Loc_scores in each protein frame
Agg_loc = pd.DataFrame(Agg_loc)
Agg_loc.reset_index(drop = False, inplace = True)
max_agg_Loc_score = Agg_loc.loc[Agg_loc.groupby(["Protein"])["Loc_score"].idxmax()][["Protein", "Frame_x"]] #*Create a table of maximum aggregated Loc_scores for each protein

for p in max_agg_Loc_score["Protein"]:
	pearson_r = col_df.apply(lambda x: stats.pearsonr(x["z_score_Abund"], x["z_score_Loc"]))


#%%
#, Running Abundance on the prim index/pos level - This is currenly not in use.
# Quant_prim_index = pd.read_csv("Quant_prim_index.csv")
# run_list = input("What positinos to run on?").split(", ")
# def Abundance_calc(Quant_prim_name) #* The abundacne could be incorporated in the post-quant function, but it could slow down. It also makes the process less modular and therefore makes timing of computer location a restraint
# 	def Abundance_conv_f(cell_median):
# 			Abundance_prot = math.log(cell_median) #Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset
# 			return(Abundance_prot)

# 	pos = Quant_prim_name[6:-8]# Get position name. Quant_|d0218r1p870300|_ALL.csv

# 	Quant_prim = pd.read_csv(Quant_prim_name)
# 	# Quant_prim["Abundance"] = pd.Series(Quant_prim["factor_median_OBJ_GFP"]).apply(Abundance_conv_f)

# 	Quant_prim["Abundance"] = pd.Sereis(Quant_prim["factor_median_OBJ_GFP"]).apply(lambda math.log(x)) #* Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset
# 	#? THis might be slightly faster because the function does not need to be declared for each core/run

# 	Quant_prim.to_csv(Quant_prim_name)
# 	return(f"{pos} complete")

# Quant_movement["Abundance"] = pd.Series(Quant_movement["factor_median_OBJ_GFP"]).apply(Abundance_calc)

# # def secondary_calc_manger(pos)-> None:
# # 	#! if goint this route, must add in the Abundance_calc funtion so that it is accessible to all paralllel cores
# # 	pos_quant_mov = pd.read_csv(f"Quant_movement_{pos}.csv")
# # 	pos_quant_mov["Abundance"] = pd.Series(pos_quant_mov["factor_median_OBJ_GFP"]).apply(Abundance_calc)
# # 	return(f"{pos} abundane calculation is complete")

# Parallel(n_jobs=pr)(delayed(secondary_calc_manger)(p) for p in positions)

# TODO: Modify to allow for mutiple chanels to be Abundance checked. It will be a fairly simple change, but should make decision when to be running the function
#, This requires modification to apply to multiple fluorophores
#* This version could be applied to all fluorphores for the next expirimentParallel(n_jobs=pr)(delayed(tdim_metric_auto)(c) for c in cols)
# fluor = input("What flourophore is this for?")
# def variable_look(fluor):  #TODO: incorporate this into another function like the movement.
# 	Quant_movement = pd.read_csv("")
# 	Quant_movement[f"Abundance_{fluor}"] = pd.Series(Quant_movement[f"factor_median_OBJ_{fluor}"]).apply(Abundance_calc)
# 	return(Quant_movement)
# #* could also run through the fluors in the library and calculate

# for fluor in flourescent:
# 	Quant_movement = pd.read_csv("")
# 	Quant_movement[f"Abundance_{fluor}"] = pd.Series(Quant_movement[f"factor_median_OBJ_{fluor}"]).apply(Abundance_calc)










