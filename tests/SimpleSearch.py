
#%%
#* Run through the indexing for images, original masks, trackmasks, and expanded masks, and quantification
from single_cell_reloc_paraquet.Pre_track import *
from single_cell_reloc_paraquet.Post_quant import *
from single_cell_reloc_paraquet.Pre_seg import *
import single_cell_reloc_paraquet.global_functions.global_variables as gv
import pandas as pd
import pickle as pk
import json
import os

#%%
#, Define the global comparison logic
#. The below version is optimized to run on a single folder.
def comparison_logical(ref1_df, ref2_df, col_n): #* This is the refercen search for all versions
	#* Call from known:
	#< ref1_df = pd.read_csv(f"{ref1}.csv")
	#< ref2_df = pd.read_csv(f"{ref2}.csv")

	#* Call respective indexer based on defined global variables
	#< exec("ref1_df = {ref1}_er()")
	#< exec("ref2_df = {ref2}_er()")

	merged = pd.merge(ref1_df, ref2_df, left_index= True, right_index= True, how = 'left') #* This may not be best as a left-side comparison
	subset_temp = merged.loc[(merged[col_n]).isna()]
	subset_temp.to_parquet(f'{ref1}_{ref2}_comparison.parquet')
	return(subset_temp, ref1_df, ref2_df)

#* The below were removed, to make the file more concise.
#< def img_seg(ref1, ref2):
#< def seg_track(ref1, ref2):
#< def track_exp(ref1, ref2):
#< def exp_quant(ref1, ref2):

#, Fill Segmentation: (identified by the missing original masks relative to the imgIndex)
#! Pass the missing to the snakemake object to fill the missing
# def seg_fill(microfluidcs_results, subset = False):
	missing_seg_fill = comparison_logical(imgIndex, orgAllmasks)

	#todo: complete the subset logic. Not needed for current run.
	#< if subset == False:
	#< 	pass
	#< else:
	#< 	missing_seg_fill = missing_seg_fill.loc[missing_seg_fill[Global_variables["Subset_by"]].isin(Global_variables["Subset"])]

	seg_fill = 1
	while seg_fill < 3:

		#> Call to snakemake.
		seg_fill += 1
		print(f"Run {seg_fill} of seg_fill complete")

	#* Remove any positions which are not completed in this part. There is no point trying to force and carry through bad original data which can't be segmented

#, Fill tracking masks (identified as the missing Allmasks relative to the orgAllmasks)
#! Pass the missing to the trackmanager as a subset. Run twice to make sure that all spots have been fillled. This is in case there was some memory overflow or other unkown cause of failure
#%%
def convert_fluor_to_folderName(fluors):
	fluors = ["GFP", "mKO", "mKa"]
	for n,f in enumerate(fluors):
		if n == 0:
			concat = f
		else:
			concat += "_" + f
	return (concat)

def seg_fill(microfluidcs_results, subset = False, fluors_folder = "GFP_mKO_mKa"):
	tm_fill = 1
	while tm_fill < 3:
		imgIndex = imgIndex_er() #* This is the index of the images
		orgAllmaks = orgAllmak_er() #* This is the index of the original masks

		subset_temp = comparison_logical(orgAllmasks, orgAllmasks)

		#! The functin will the called again as the looper has only be increased by 1
		track_manager(temp_subset)
		tm_fill += 1
		print(f"Run {tm_fill} of org_fill complete")


#, Fill the expansion mask (identified as the missing expmasks relative to the Allmasks)
def exp_fill(microfluidcs_results, subset = False):
	em_fill = 1
	while em_fill < 3:
		Allmasks = Allmasks() #* This is the index of the Allmasks
		expAllmasks = expAllmasks() #* This is the index of the expanded masks
		comparison_logical(Allmasks, expAllmasks)
		#! The functin will the called again as the looper has only be increased by 1
		expansion_manager(temp_subset)
		em_fill += 1
		print(f"Run {em_fill} of exp_fill complete")


#, Fill the quantification (identified as missing quantification relative to the expanded masks)
def quant_fill(microfluidcs_results, subset = False):
	q_fill = 1
	while q_fill < 3:
		expAllmasks = expAllmasks() #* This is the index of the expanded masks
		quantification_index = quantification_index() #* This is the index of the quantification

		expAllmasks.set_index("Unique_frame")
		quantification_index.set_index("Unique_frame")

		comparison_logical(expAllmasks, quantificatoin_index)
		q_fill += 1
		print(f"Run {q_fill} of quant_fill complete")

#%%
#Todo: Automate the simple failures (non seg issues) to complete the dataset
if __name__ == "__main__":
	try: #* This is to load a json file. May change to load a pickle instead
		Global_variables = json.load("Global_variables.json", "r")
	except:
		print("Global_variables.json does not exist. Please create one.")
		Global_variables = gv.global_manager()

	temp_output = gv.slash_switch(input("What should the temp output be called?"))
	temp_output_path = os.path.join(Global_variables["microfluidcs_results"], temp_output)
	os.makedirs(temp_output_path, exist_ok=True)
	os.chdir(temp_output_path) #* The functions should be changed to the

	# seg_state = input("Run segmentation? True or False?")
	# if seg_state == "True":
	# 	seg_state = True
	# if seg_state == "False":
	# 	seg_state = False

	track_state = input("Run tracking? True or False?")
	if track_state == "True":
		track_state = True
	if track_state == "False":
		track_state = False

	quant_state = input("Run quantification? True or False?")
	if quant_state == "True":
		quant_state = True
	if quant_state == "False":
		quant_state = False

	# if seg_state == True:
	# 	seg_fill(Global_variables["microfluidcs_results"], subset = False, fluors_folder = convert_fluor_to_folderName(Global_variables["Fluor_names"]))
	if track_state == True:
		exp_fill(Global_variables["microfluidcs_results"])
	if quant_fill == True:
		quant_fill(Global_variables["microfluidcs_results"])