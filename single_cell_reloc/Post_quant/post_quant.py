#%%
import datetime
from scipy import stats
from joblib import Parallel, delayed
from glob import glob
from scipy.stats import variation
import single_cell_reloc.Post_quant.post_quant_JAN20_2023 as post_quant_funcs
import single_cell_reloc.Post_quant.Strain_ID as Strain_ID_funcs
import os
import pandas as pd

def Post_quant_manager(percentiles, mulitplex) -> None:
	Quant_ALL_index = post_quant_funcs.Quantification_index_er()

	if mulitplex == True:
		Strain_ID = Strain_ID_funcs.Strain_ID_multiplex()
	elif mulitplex == False:
		Strain_ID = Strain_ID_funcs.Strain_ID_single()

	l = len(Quant_ALL_index)
	if l < Global_variables:
		use_cores = l
	else:
		use_cores = Global_variables["cpu_se"]

	Parallel(n_jobs=use_cores, verbose = 100)(delayed(Strain_ID)(p) for p in range(len(Quant_ALL_index)))

	os.chdir(Global_variables["Post_path"])

	Quant_prim_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith("mary.csv") and name.startswith("Quant"): # fix naming
				Quant_prim_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		break #This makes the program run non recursively and not decend into daughter folders

	Quant_prim_index = pd.DataFrame(Quant_prim_index)
	Quant_prim_index["PositionID"] = pd.Series(Quant_prim_index.iloc[:,0]).apply(f_Position_ID_qALLi)
	# Quant_prim_index["Position"] = pd.Series(Quant_prim_index["Path"]).apply(f_Position_ID_qALLi)
	# Quant_prim_index["Run"] = pd.Series(Quant_prim_index["Position"]).apply(f_run)
	# Quant_prim_index["Col"] = pd.Series(Quant_prim_index["Position"]).apply(f_col)
	# Quant_prim_index["Date"] = pd.Series(Quant_prim_index["Position"]).apply(f_expdate)

	Quant_prim_index.sort_values(by = "PositionID", inplace = True)
	Quant_prim_index.to_csv("Quant_prim_index.csv")# , index = False)

	# for p in percentiles:



if __name__ == "__main__": #* This should be included in all major folders
	import single_cell_reloc.global_functions.global_variables as glbl_vars
	Global_variables = glbl_vars.global_vars()
	Post_quant_manager(percentiles=Global_variables["percentiles"], mulitplex= Global_variables["mulitplex"] )














# %%
