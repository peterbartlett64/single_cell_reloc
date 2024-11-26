#
#%% Read in the modulesX
import pandas as pd
import os
import single_cell_reloc_paraquet.global_functions.global_variables as gv
# import single_cell_reloc_paraquet.Post_quant.Myo_ID_different as mi
# import single_cell_reloc_paraquet.Post_quant.Move as mv

#%% Read in files
# mplx_state = input("Are this sample multiplexed?")
# if mplx_state == "yes" or mplx_state == "Yes" or mplx_state == "YES" or mplx_state == "y"
# 	mplx_state = True
# else:
# 	mplx_state = False

#. Each of these files should be a version of Quant_ALL (format)
try:
	ref_path #* Check to see if the variable has been defined by calling it
except NameError:
	ref_path = gv.slash_switch(input("Path to reference file")) #* This is just a simpified call to os.path(input, sep = '\')
	compare_path = gv.slash_switch(input("Path to comparative file"))

reference_normal = pd.read_parquet(ref_path).set_index(['Cell_Barcode', "Unique_Frame"])
reference_compare = pd.read_parquet(compare_path).set_index(['Cell_Barcode', "Unique_Frame"])
#%% Myo_id
# if mplx_state == True:
# 	reference_normal = mi(reference_normal)
# 	reference_compare = mi(reference_compare)
# else:
# 	pass

# reference_normal = mv(reference_normal) #* Calculte the movement for the normal file
# reference_compare = mv(reference_compare) #* Calculte the movement for the compare file

#%% Create comparison and apply logic. Main part of script
comparison = pd.merge(reference_normal, reference_compare, left_index= True, right_index= True, suffixes=(None,'_compare'))
comparison = comparison.reindex(sorted(comparison.columns), axis=1)
# comparison.sort_index(axis=1, inplace=True)

#%% Apply check_funcs

def check_dif(x, col):
	try:
		y = x[col]- x[f"{col}_compare"]
	except TypeError:
		if comparison[x] == comparison[f"{x}_compare"]:
			y = False
		else:
			y = True
	return(y)

cols_list = reference_normal.columns


for x in cols_list:
	comparison[f"{x}_diff"] = comparison[x].transform(check_dif(col = x)) #* This is the version for a function call
	# try:
	# 	comparison[f"{x}_diff"] = comparison[x] - comparison[f"{x}_compare"] #* This is the version for a vector subtraction of rows
	# except TypeError:
	# 	if comparison[x] == comparison[f"{x}_compare"]:
	# 		comparison[f"{x}_diff"] = False
	# 	else:
	# 		comparison[f"{x}_diff"] = True

#* Confirm that the reloc results are the same
#. It is probably easier to run the post-loc individually right now
col