#%%
import datetime
from scipy import stats
from joblib import Parallel, delayed
from glob import glob
from scipy.stats import variation
# import single_cell_reloc.Post_quant.post_quant_JAN20_2023 as post_quant_funcs
import single_cell_reloc.Post_quant.Strain_ID as Strain_ID_funcs
import single_cell_reloc.Post_quant.File_indexers as pq_files
import single_cell_reloc.Post_quant.Move as Move
import os
import pandas as pd

def Post_quant_strain_manager(percentiles, mulitplex, cores) -> None:
	Quantification_index = pq_files.Quantification_index_er()
	positions = Quantification_index["PositionID"].unique()

	max_set_cores = Global_variables["cpu_se"] #* This is the maximum number of cores 'available' for use as set in the global variables
	l = len(positions)
	if l < Global_variables:
		use_cores_len = l
	else:
		use_cores_len = max_set_cores

	Parallel(n_jobs=use_cores_len, verbose = 100, prefer='threads')(delayed(pq_files.combine_pos)(p) for p in positions) #* This should prefer threads as it largely IO limited.

	Quant_ALL_index = pq_files.Quant_ALL_index_er()


	os.chdir(Global_variables["microfluidics_results"]) #* This should already be globally defined as a variable. It will either be created in
	condition_information = pq_files.read_condition_informaiton()

	if mulitplex == True:
		Strain_ID = Strain_ID_funcs.Strain_ID_multiplex() #Create a variable synonym
	elif mulitplex == False:
		Strain_ID = Strain_ID_funcs.Strain_ID_single()

	Parallel(n_jobs=cores, verbose = 100, prefer='threads')(delayed(Strain_ID)(p) for p in range(len(Quant_ALL_index))) #? Determine if this is mostly IO limited or computation

	for p in percentiles:
		Move.I_move(p)
		

	return()#* This is to automate running both percentages






if __name__ == "__main__": #* This should be included in all major folders
	import single_cell_reloc.global_functions.global_variables as glbl_vars
	Global_variables = glbl_vars.global_manager() #* This includes all the checks to make sure that all params passed are permissible
	Post_quant_manager(percentiles=Global_variables["percentiles"], mulitplex= Global_variables["mulitplex"] )














# %%
