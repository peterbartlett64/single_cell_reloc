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
if __name__ == '__main__':
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze', #* This is the location of the images
	'microfluidics_results': 'E:/Microfluidics/RESULTS', #* This is the location of the segmentation,tracking and pre-quant index files
	'post_path': 'D:/ALL_FINAL', #gv.slash_switch(input('What is the post quant directory')) , #Todo: This needs to be changed to a input call
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

	# All_protein_files = sorted(glob('*.parquet'))
	# All_proteins_df = pd.concat((pd.read_parquet(file) for file in All_protein_files), ignore_index= True)

#%%
#, NEW version of Abundance genration.
#! This was moved into 'Intensity_select_percentage_combine.py' before
#. Decided that it would be better to call this file as a module instead

def Abundance_log_manager(local_prot_df:pd.DataFrame): #* The abundacne could be incorporated in the post-quant function, but it could slow down. It also makes the process less modular and therefore makes timing of computer location a restraint
	#, Calculate the abundance for every cell at every timepoint. This is done with such density so that associations can be made to the current abundance and the Loc_scores
	local_prot_df["Abundance"] = local_prot_df['factor_median_OBJ_GFP']
	local_prot_df["log_Abundance"] = pd.Series(local_prot_df["factor_median_OBJ_GFP"]).apply(lambda x: math.log(x)) #* Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset so it is still an arbitraty unit value

	local_prot_df['log_Loc_score'] = pd.Series(local_prot_df["Loc_score"]).apply(lambda x: math.log(x))
	group_prot_frame = local_prot_df.groupby(["Protein", "Frame"])

	#* The z.score is now a comparison between the both log transformed values. Abundance does not scale linearly with intensity so it is required for conversion of the median factor. Obviouly cannot do linear regression between the a reg vs log transformed scale
	local_prot_df["z_score_logLoc"] = group_prot_frame["log_Loc_score"].transform(stats.zscore)
	local_prot_df["z_score_logAbund"] = group_prot_frame["log_Abundance"].transform(stats.zscore)

	#* Added these lines below for comparison to the original values
	local_prot_df["log_Abundance"] = local_prot_df['factor_median_OBJ_GFP']
	local_prot_df["z_score_Loc"] = group_prot_frame["Loc_score"].transform(stats.zscore)
	local_prot_df["z_score_Abund"] = group_prot_frame["Abundance"].transform(stats.zscore)
	return(local_prot_df)

c = 1
for filename in os.listdir(percentage_combined_path):
	c += 1
	if filename.endswith('.parquet'):
		# Load the Parquet file into a DataFrame
		file_path = os.path.join(percentage_combined_path, filename)
		df = pd.read_parquet(file_path).reset_index(drop = False)
		df = Abundance_log_manager(df)
		df.to_parquet(file_path)

		# Append the data to the merged_data DataFrame
		# merged_data = pd.concat([merged_data, df], ignore_index=True)
		print(c)


# #. The below function is just for getting back the regular values/= This is not required as the factor is already stored under a different name and Loc is in the df
# def Abundance_fix_manager(local_prot_df:pd.DataFrame): #* The abundacne could be incorporated in the post-quant function, but it could slow down. It also makes the process less modular and therefore makes timing of computer location a restraint
# 	#, Calculate the abundance for every cell at every timepoint. This is done with such density so that associations can be made to the current abundance and the Loc_scores
# 	local_prot_df["Abundance"] = pd.Series(local_prot_df["factor_median_OBJ_GFP"]).apply(lambda x: math.exp(x)) #* Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset
# 	#* purposely overwrit the old version. Because the dataframes are quite large, it is excessive to keep the unmodified df in memory for a couple non-mutative lines
# 	#? Perfroming with the lambda function might be slightly faster because the function does not need to be declared for each core/run

# 	group_prot_frame = local_prot_df.groupby(["Protein", "Frame"])


# 	# #* This seems like a good plac to grab the z_scores for the Loc_score and Abundance
# 	local_prot_df["z_score_Loc"] = group_prot_frame["Loc_score"].transform(stats.zscore)
# 	local_prot_df["z_score_Abund"] = group_prot_frame["Abundance"].transform(stats.zscore)
# 	return(local_prot_df)