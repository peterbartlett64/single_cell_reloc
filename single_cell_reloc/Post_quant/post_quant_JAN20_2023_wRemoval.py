#%%
#, This document was modified from the AUG2 tree. This is a major modifiacation to make try and deal with the DivideByZero error which removes some cells
#? These are a lot of imports. I doubt that I need them all
from datetime import datetime
from decimal import DivisionByZero
import math
import pandas as pd
import numpy as np
import glob
import os
import datetime
from scipy import stats
from joblib import Parallel, delayed
from glob import glob
from scipy.stats import variation

#%%

# post_path =  "C:/Users/Peter/Desktop/Subset_Images/JAN2/2022-01-03"
microfluidic_results = Global_variables["microfluidic_results"]
post_path= Global_variables["post_path"]
pr = Global_variables["cpu_se"]

percentile = int(input("What percentile would you like to run on"))

#Change the paths to the corrected version becuase of the windows input
microfluidic_results = slash_switch(microfluidic_results)
post_path = slash_switch(post_path)
os.chdir(post_path) # This would normally be the path set in the quantification scipt

#%%
Experimental_info = {
	"timepoint_space": 7.5,
	"cell_tracking": True,
	"flourescent": True,
	"data_subset": False
	}

if Experimental_info["data_subset"] == True:
	subdat = str(input("Enter position barcodes to be analyzed (comma seperated and in the form dMMDDpPosition)"))

	instances = subdat.split(", ")
	print(f"Analysis will be performed for positions which correspond to the following pos_barcodes: {instances}")
else:
	pass

try: subdat
except NameError: pass
else: del subdat

#%% # Below was just copied from Quantificaiton.py becuase it was working at last check.
def Quantification_index_er():
	#! commented out because only Quant_ALL files are present
	Quantification_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith(".csv") and name.startswith("Quantification_d"):
				Quantification_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		# break #This makes the program run non-recursively and not decend into daughter folders

	Quantification_index = pd.DataFrame(Quantification_index)

	year = str(datetime.datetime.now().year)
	decade = year[:3] # This assumes that the analysis is done within the same decade as starting
	del year

	def f_Position_ID(z):
		start = z.find('tion_')+5 #Note: The shift forward is dependant upon how you write out position
		end = z.find("_"+ decade)-5   #Assume that this pipeline will only be used this century! Make sure that the 'n' is present in the array to confirm no place 0 has been lost
		return(z[start:end])

	def f_Frame(z):
		start = z.find('tion_')+5
		end = z.find('_' + "20")
		return(z[start:end])

	def f_expdate(x):
		dend = x.find('r')
		expdate = x[0:dend]
		return(expdate)

	def Mod_epoch(f):
		t = os.path.getmtime(f)
		return(t)
		# return(time.ctime(t))

	def f_non_sync(p):
		substring = ".sync"
		if substring in p:
			return(None)
		else:
			return(p)

	def epoch_convert(timestamp):
		date_time = datetime.datetime.fromtimestamp(timestamp)
		d = date_time.strftime("%m/%d/%Y")
		return(d)

	# Quantification_index["Pad"] =
	Quantification_index["PositionID"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Position_ID)
	Quantification_index["Date"] = pd.Series(Quantification_index["PositionID"]).apply(f_expdate)
	Quantification_index["Frame"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Frame)

	#* NEWLY ADDED
	Quantification_index["Mod_epoch"] = pd.Series(Quantification_index["Path"]).apply(Mod_epoch)
	Quantification_index["Mod_date"] = pd.Series(Quantification_index["Mod_epoch"]).apply(epoch_convert)

	Quantification_index["Max_epoch_frame"] = Quantification_index.groupby("Frame")["Mod_epoch"].transform('max') #? This isn't sylish but works
	Quantification_index = Quantification_index.loc[Quantification_index["Mod_epoch"] == Quantification_index["Max_epoch_frame"]]
	Quantification_index.drop(columns=["Max_epoch_frame"])


	if Experimental_info["data_subset"] == True:
		Quantification_index = Quantification_index[Quantification_index["PositionID"].isin(instances)]
		print(f"Post-quant will be performed for positions which correspond to the following pos_barcodes: {instances}")
	else:
		pass
	Quantification_index.sort_values(by = "PositionID", inplace = True)
	Quantification_index.to_csv("Quantification_index.csv")
	return(Quantification_index)

####### Function run
# Quantification_index = Quantification_index_er()
#######

# Quantification_index = pd.read_csv("Quantification_index.csv") # Temp added because want to read the version whiich was created pre-quant today

#%% Combine frames into posision
from joblib import Parallel, delayed

def combine_pos(pos):
	Quant_frame_comb = pd.DataFrame([])
	subset = Quantification_index.loc[Quantification_index["PositionID"] == pos]
	subset_s = subset.sort_values(by = "Frame")
	for qf in range(len(subset)):
		try:
			q = pd.read_csv(subset_s.iloc[qf,0])
		except:
			continue
		Quant_frame_comb = pd.concat([Quant_frame_comb, q])
	Quant_frame_comb.to_csv(f"Quant_{pos}_ALL.csv")
	return(f"{pos} frame merging is complete")

positions = Quantification_index["PositionID"].unique()


l = len(positions)
if  l > pn:
	pr = pn
else:
	pr = l

# for pos in positions:
# 	Quant_frame_comb = pd.DataFrame([])
# 	subset = Quantification_index.loc[Quantification_index["PositionID"] == pos]
# 	subset.sort_values(by = "Frame", inplace = True)
# 	for qf in range(len(subset)):
# 		q = pd.read_csv(subset.iloc[qf,0])
# 		Quant_frame_comb = pd.concat([Quant_frame_comb, q])
# 	Quant_frame_comb.to_csv(f"Quant_{pos}_ALL.csv")
# 	print(f"{pos} merging is complete")
####### Function run
Parallel(n_jobs=pr, verbose = 100)(delayed(combine_pos)(p) for p in positions)
#######


#%% Generate an index of Quant_ALL
# # Create an index of Quantification dfs
############################################

year = str(datetime.datetime.now().year) #this assumes the analysis is completed within the same decade of starting. A reasonable expectation, but something to check if running aroudn the new year
decade = year[:3]
del year

def f_Position_ID_qALLi(z):
	start = z.find('ant_')+4
	end = z.find("_ALL")
	return(z[start:end])

# def Quant_ALL_index_er():
Quant_ALL_index  = []
count = 0
for root, dirs, files, in os.walk(os.getcwd()):
	for name in files:
		if name.endswith("_ALL.csv") and name.startswith("Quant"): # fix naming
			Quant_ALL_index.append({'Path': os.path.join(root, name)})
			count = count + 1
			print(count, end="\r")
		else:
			pass
	break #This makes the program run non-recursively and not decend into daughter folders

Quant_ALL_index = pd.DataFrame(Quant_ALL_index)
Quant_ALL_index["PositionID"] = pd.Series(Quant_ALL_index.iloc[:,0]).apply(f_Position_ID_qALLi)
Quant_ALL_index.sort_values(by = "PositionID", inplace = True)
Quant_ALL_index.to_csv("Quant_ALL_index.csv")

###### Function run
Quant_ALL_index = Quant_ALL_index_er()
######

#%%
# Cell_index_timing = info_simple[["track_index", "track_start_frame"]]
# os.chdir(microfluidics_results)
Quant_ALL_index = pd.read_csv("Quant_ALL_index.csv") #*#* Does the file need to read in to work?
Quant_ALL_index.drop(Quant_ALL_index.columns[Quant_ALL_index.columns.str.contains('Unnamed',case = False)],axis = 1, inplace= True)
#%%
os.chdir(microfluidic_results)
Condition_information = pd.read_excel("MicrofluidicsMap_wCol.xlsx", sheet_name='ProteinLocations', dtype={'Date' : str})
Condition_information = Condition_information.drop(columns=["Notes", "Predicted localization Change"]).dropna()
Condition_information.drop(Condition_information.columns[Condition_information.columns.str.contains('Unnamed',case = False)],axis = 1, inplace= True)
#%%
# year = str(datetime.date.today().year)
time_per_frame = Experimental_info['timepoint_space'] ###### This is set and could be placed in a more convient place.
year = str(datetime.datetime.now().year) #this assumes the analysis is completed in the same decade as the aquisition
decade = year[:3]
del year

def f_convert_date(d):
	d = str(d)
	s = d.find(decade)+5
	return("d" + d[s:])

Condition_information["Date"] = pd.Series(Condition_information["Date"]).apply(f_convert_date)
# Condition_information["Time"] = Condition_information["Time.1"] # This is required if there is a "Time" for aquisition, and timing for treatment
# Condition_information = Condition_information[["Date", "Run Number", "Time", "MapID (Col_Range)", "Myo1Marker", "Protein", "Time"]]
#%%
## Below is a far more complex implementation of Myo1 determination
#Instead of being just based on the max and min of each channel, this version also takes into account the start frame and also the dynamic range
os.chdir(post_path)


#%%
def f_Position_ID_qALLi(z):
	start = z.find('ant_')+4
	end = z.find("_prim")
	return(z[start:end])

# def f_col(p):
# 	Col = int(p[p.find("p")+1:-4])
# 	return(Col)

# def f_expdate(x):
# 	dend = x.find('r')
# 	expdate = x[0:dend]
# 	return(expdate)

# def f_run(x):
# 	s = x.find("r")+1
# 	e = s+1
# 	r = x[s:e]
# 	return(r)


##### MOVED THE FRAME JUSTIFICATION UPSTREAM
#%%
####### THIS IS ALL CELLS USED. THE ONLY VERSION THAT DOES NOT IS DEC16_SUBSET
####### This is the two sides Locp = 1_score test
# percentile = str(input("What is the the parameter to check Loc Change on?"))

def I_move(p, percentile = percentile):
	try: #TODO: Rename Frame_x to just Frame
		Quant_prim = pd.read_csv(Quant_prim_index.iloc[p,0], usecols = ['Cell_Barcode', 'ImageID', 'Date', 'Frame_x', 'Unique_Frame', 'factor_median_OBJ_GFP', 'factor_mean_OBJ_GFP', 'x90thPercentile_norm_OBJ_Median_GFP', 'x95thPercentile_norm_OBJ_Median_GFP', 'x99thPercentile_norm_OBJ_Median_GFP', 'Progen_bud', 'x90thPercentile_GFP_RAW', 'x95thPercentile_GFP_RAW', 'x99thPercentile_GFP_RAW', 'max_GFP_RAW', 'max_mKa_RAW', 'max_mKO_RAW', 'TrackID_valid', 'Myo1Identity', 'mKO_foldChange', 'mKO_direction', 'mKA_foldChange', 'mKa_direction', 'byProgen_bud', 'byRange', 'Col_info', 'Protein', 'Is_treated', 'Frames_post_treatment', #* From here on are newly added columns to determine the cell stage. Could do a pd.merge at the end instead.
		"x80thPercentile_Diff_background_mKate", "x90thPercentile_Diff_background_mKate", "x99thPercentile_Diff_background_mKate", "x80thPercentile_Diff_background_mKO", "x90thPercentile_Diff_background_mKO", "x99thPercentile_Diff_background_mKO", "averageIntensity_mKO_Frame", "averageIntesntiy_mKO_Background", "averageIntensity_mKO_Object", "mKO_spread", "averageIntensity_mKate_Frame", "averageIntesntiy_mKate_Background", "averageIntensity_mKate_Object", "mKate_spread", "factor_median _OBJ_KO", "factor_mean_OBJ_KO", "factor_total_OBJ_KO", "factor_mKO_background_Med", "factor_mKO_background_Avg", "factor_mKO_background_Tot", "factor_median_OBJ_mKate", "factor_mean_OBJ_mKate", "factor_total_OBJ_mKate", "factor_mKate_background_Med", "factor_mKate_background_Avg", "factor_mKate_background_Tot", "x60thPercentile_mKa_RAW", "x80thPercentile_mKa_RAW", "x90thPercentile_mKa_RAW", "x95thPercentile_mKa_RAW", "x99thPercentile_mKa_RAW", "max_mKa_RAW", "x60thPercentile_mKO_RAW", "x80thPercentile_mKO_RAW", "x90thPercentile_mKO_RAW", "x95thPercentile_mKO_RAW", "x99thPercentile_mKO_RAW", "max_mKO_RAW"])

		Quant_prim["Unique_pos"] = pd.Series(Quant_prim["Unique_Frame"]).apply(lambda x: x[:x.find('f')]) #* This was added on AUG31,2022

		Pos = Quant_prim_index.iloc[p,1]
		Pos = str(Pos)
		Quant_prim.drop(Quant_prim.columns[Quant_prim.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
		# try:
		MyoSeen = Quant_prim["Myo1Identity"].unique()
		if "Myo1_mKO" in MyoSeen:
			mKO_skip = False
		else:
			mKO_skip = True

		if "Myo1_mKa" in MyoSeen:
			mKa_skip = False
		else:
			mKa_skip = True

		#This is for the removal of dead cells based on high mKa and low GFP relative to the population.  The code is not complete and will remove to many cells
		Quant_f__twe = Quant_prim.loc[(Quant_prim["Frame_x"] <= 20) & (Quant_prim["Is_treated"] == 0)].copy()
		Median_pop = Quant_f__twe.groupby(by = "Cell_Barcode").agg('min').groupby(by = "Myo1Identity").agg('median')
		Std_pop = Quant_f__twe.groupby(by = "Cell_Barcode").agg('min').groupby(by = "Myo1Identity").agg('std')

		if mKa_skip == False:
			MedGFPm2std_pop_mKa = Median_pop.loc["Myo1_mKa", "max_GFP_RAW"] + 2*Std_pop.loc["Myo1_mKa", "max_GFP_RAW"]
			MedmKap2std_pop_mKa = Median_pop.loc["Myo1_mKa", "max_mKa_RAW"] + 2*Std_pop.loc["Myo1_mKa", "max_mKa_RAW"]
			MedmKOp2std_pop_mKa = Median_pop.loc["Myo1_mKa", "max_mKO_RAW"] + 2*Std_pop.loc["Myo1_mKa", "max_mKO_RAW"]

			Low_GFP_pop_mKa = Quant_f__twe[(Quant_f__twe["max_GFP_RAW"] > MedGFPm2std_pop_mKa) & (Quant_f__twe["Myo1Identity"] == "Myo1_mKa")]["Cell_Barcode"]
			High_mKa_pop_mKa = Quant_f__twe[(Quant_f__twe["max_mKa_RAW"] > MedmKap2std_pop_mKa) & (Quant_f__twe["Myo1Identity"] == "Myo1_mKa")]["Cell_Barcode"]
			High_mKO_pop_mKa = Quant_f__twe[(Quant_f__twe["max_mKO_RAW"] > MedmKOp2std_pop_mKa) & (Quant_f__twe["Myo1Identity"] == "Myo1_mKa")]["Cell_Barcode"]
		else:
			pass

		if mKO_skip == False:
			MedGFPm2std_pop_mKO = Median_pop.loc["Myo1_mKO", "max_GFP_RAW"] + 2*Std_pop.loc["Myo1_mKO", "max_GFP_RAW"]
			MedmKap2std_pop_mKO = Median_pop.loc["Myo1_mKO", "max_mKa_RAW"] + 2*Std_pop.loc["Myo1_mKO", "max_mKa_RAW"]
			MedmKOp2std_pop_mKO = Median_pop.loc["Myo1_mKO", "max_mKO_RAW"] + 2*Std_pop.loc["Myo1_mKO", "max_mKO_RAW"]
			Low_GFP_pop_mKO = Quant_f__twe[(Quant_f__twe["max_GFP_RAW"] > MedGFPm2std_pop_mKO) & (Quant_f__twe["Myo1Identity"] == "Myo1_mKO")]["Cell_Barcode"]
			High_mKa_pop_mKO = Quant_f__twe[(Quant_f__twe["max_mKa_RAW"] > MedmKap2std_pop_mKO) & (Quant_f__twe["Myo1Identity"] == "Myo1_mKO")]["Cell_Barcode"]
			High_mKO_pop_mKO = Quant_f__twe[(Quant_f__twe["max_mKO_RAW"] > MedmKOp2std_pop_mKO) & (Quant_f__twe["Myo1Identity"] == "Myo1_mKO")]["Cell_Barcode"]
		else:
			pass

		# if MedGFPm2std_pop_mKO < MedGFPm2std_pop_mKa:
		# 	Low_GFP_pop_Cannot = Quant_f__twe[(Quant_f__twe["max_GFP_RAW"] > MedGFPm2std_pop_mKO) & (Quant_f__twe["Myo1Identity"] == "Cannot Be Determined")]["Cell_Barcode"]
		# 	High_mKa_pop_Cannot = Quant_f__twe[(Quant_f__twe["max_mKa_RAW"] > MedmKap2std_pop_mKO) & (Quant_f__twe["Myo1Identity"] == "Cannot Be Determined")]["Cell_Barcode"]
		# elif MedGFPm2std_pop_mKO > MedGFPm2std_pop_mKa:
		# 	Low_GFP_pop_Cannot = Quant_f__twe[(Quant_f__twe["max_GFP_RAW"] > MedGFPm2std_pop_mKa) & (Quant_f__twe["Myo1Identity"] == "Cannot Be Determined")]["Cell_Barcode"]
		# 	High_mKa_pop_Cannot = Quant_f__twe[(Quant_f__twe["max_mKa_RAW"] > MedmKap2std_pop_mKa) & (Quant_f__twe["Myo1Identity"] == "Cannot Be Determined")]["Cell_Barcode"]

		if mKa_skip == False and mKO_skip == False:
			if MedGFPm2std_pop_mKa > MedGFPm2std_pop_mKO:
				gGFP_sel = MedGFPm2std_pop_mKa
			elif MedGFPm2std_pop_mKO > MedGFPm2std_pop_mKa:
				gGFP_sel = MedGFPm2std_pop_mKO
			else:
				gGFP_sel = MedGFPm2std_pop_mKa

			if MedmKap2std_pop_mKa > MedmKap2std_pop_mKO:
				gKa_sel = MedmKap2std_pop_mKa
			elif MedmKap2std_pop_mKO > MedmKap2std_pop_mKa:
				gKa_sel = MedmKap2std_pop_mKO
			else:
				gKa_sel = MedmKap2std_pop_mKa

		elif mKa_skip == False and mKO_skip == True:
			gGFP_sel = MedGFPm2std_pop_mKa
			gKa_sel = MedmKap2std_pop_mKa
		elif mKa_skip == True and mKO_skip == False:
			gGFP_sel = MedGFPm2std_pop_mKO
			gKa_sel = MedmKap2std_pop_mKO
		else:
			return(f"Failure to detect/store either population")

		Quant_f__twe["min_max_GFP_Object"] = Quant_f__twe.groupby(by = "Cell_Barcode").max_GFP_RAW.transform('min')
		Quant_f__twe["min_max_mKate_Object"] = Quant_f__twe.groupby(by = "Cell_Barcode").max_mKa_RAW.transform('min')
		# Quant_yes = Quant_f__twe[(Quant_f__twe["min_max_GFP_Object"] < gGFP_sel) & (Quant_f__twe["min_max_mKate_Object"] <gKa_sel)]["Cell_Barcode"]
		Quant_no = Quant_f__twe[(Quant_f__twe["min_max_GFP_Object"] > gGFP_sel) & (Quant_f__twe["min_max_mKate_Object"] > gKa_sel)]["Cell_Barcode"] # Creat a list of cells which are above the GFP and mKa thresholds.

		Quant_prim = Quant_prim.loc[~Quant_prim["Cell_Barcode"].isin(Quant_no)] # This is functional but messy. Remove the cells marked as too flourescent to be considered viable. Removal is only performed BASED ON the first 10 frames
		Quant_no.to_csv(f"Dropped_dead_{Pos}.csv")

		Quant_prim.set_index(["Cell_Barcode", "ImageID"], inplace = True)

		fluor_perc = (f"x{percentile}thPercentile_norm_OBJ_Median_GFP")
		fluor_perc_r = (f"x{percentile}thPercentile_GFP_RAW")

		UNT_mKa_all = Quant_prim[(Quant_prim["Myo1Identity"] == "Myo1_mKa") & (Quant_prim["Is_treated"] == 0)]
		UNT_mKa = UNT_mKa_all[UNT_mKa_all["Frames_post_treatment"] >= -10] # This is the subset of pre-treatment, 3 frames before treatment. May be able to expand but not clear that it is needed #! TESTING THE EFFECT OF USING FRAMES MORE THAN -10 INSTEAD OF -3

		UNT_mKa_GFP = UNT_mKa[fluor_perc].dropna()
		if len(UNT_mKa_GFP) == 0:
			return(f"{Pos} has no mKa")

		# UNT_GFP_mKa_factor_upper = np.median(UNT_mKa_GFP) + 2*np.std(UNT_mKa_GFP) ## This was the median plus the regular sd. To be more inline with
		# UNT_GFP_mKa_factor_lower = np.median(UNT_mKa_GFP) - 2*np.std(UNT_mKa_GFP)

		UNT_GFP_mKa_factor_upper = np.median(UNT_mKa_GFP) + 2*stats.median_abs_deviation(UNT_mKa_GFP)
		UNT_GFP_mKa_factor_lower = np.median(UNT_mKa_GFP) - 2*stats.median_abs_deviation(UNT_mKa_GFP)


		UNT_mKO_all = Quant_prim[(Quant_prim["Myo1Identity"] == "Myo1_mKO") & (Quant_prim["Is_treated"] == 0)] # The second part seems unescessary
		UNT_mKO = UNT_mKO_all[UNT_mKO_all["Frames_post_treatment"] >= -10] # This may seem wrong  but have already selected for cells pre-treatment above #! TESTING THE EFFECT OF USING FRAMES MORE THAN -10 INSTEAD OF -3 # AUG29 2022


		UNT_mKO_GFP = UNT_mKO[fluor_perc].dropna()
		if len(UNT_mKO_GFP) == 0:
			return(f"{Pos} has no mKO")

		# UNT_GFP_mKO_factor_upper = np.median(UNT_mKO_GFP) + 2*np.std(UNT_mKO_GFP)
		# UNT_GFP_mKO_factor_lower = np.median(UNT_mKO_GFP) - 2*np.std(UNT_mKO_GFP)

		UNT_GFP_mKO_factor_upper = np.median(UNT_mKO_GFP) + 2*stats.median_abs_deviation(UNT_mKO_GFP)
		UNT_GFP_mKO_factor_lower = np.median(UNT_mKO_GFP) - 2*stats.median_abs_deviation(UNT_mKO_GFP)

		MMS_mKa = Quant_prim[(Quant_prim["Myo1Identity"] == "Myo1_mKa") & (Quant_prim["Is_treated"] == 1)]
		MMS_mKa["mKa_factor_upper"] = UNT_GFP_mKa_factor_upper
		MMS_mKa["mKa_factor_lower"] = UNT_GFP_mKa_factor_lower

		MMS_mKO = Quant_prim[(Quant_prim["Myo1Identity"] == "Myo1_mKO") & (Quant_prim["Is_treated"] == 1)]
		MMS_mKO["mKO_factor_upper"] = UNT_GFP_mKO_factor_upper
		MMS_mKO["mKO_factor_lower"]= UNT_GFP_mKO_factor_lower

		Raw_factorUpper_mKa = np.median(UNT_mKa[fluor_perc_r]) + 2*stats.median_abs_deviation(UNT_mKa[fluor_perc]) #This calculated the parameters based on the variable percentile. By writing in this way, I can test different percentile for differnt loc types
		Raw_factorLower_mKa = np.median(UNT_mKa[fluor_perc_r]) - 2*stats.median_abs_deviation(UNT_mKa[fluor_perc])
		Raw_factorUpper_mKO = np.median(UNT_mKO[fluor_perc_r]) + 2*stats.median_abs_deviation(UNT_mKO[fluor_perc])
		Raw_factorLower_mKO = np.median(UNT_mKO[fluor_perc_r]) - 2*stats.median_abs_deviation(UNT_mKO[fluor_perc])

		#Calculate raw threshold values. These values are helpful for checking whether reloc has occured. This how to derive from the value can be directly calculated as above
		# GFP_normal_avg_mKa = np.mean(UNT_mKa['factor_median_OBJ_GFP'])
		# GFP_normal_avg_mKO = np.mean(UNT_mKO['factor_median_OBJ_GFP'])
		# MMS_mKa["mKa_factor_upper"] = UNT_GFP_mKa_factor_upper * GFP_normal_avg_mKa
		# MMS_mKa["mKa_factor_lower"] = UNT_GFP_mKa_factor_lower * GFP_normal_avg_mKa
		# MMS_mKO["mKO_factor_upper"] = UNT_GFP_mKO_factor_upper * GFP_normal_avg_mKO
		# MMS_mKO["mKO_factor_lower"] = UNT_GFP_mKO_factor_lower * GFP_normal_avg_mKO

		UNT_mKa.to_csv(f"UNT_{Pos}mKa.csv") # Output the cells and information on which parameters are derived
		UNT_mKO.to_csv(f"UNT_{Pos}mKO.csv")

		def reloc_score_mKa(row):
			metric = row[fluor_perc]
			loc_score_upper = metric/UNT_GFP_mKa_factor_upper
			loc_score_lower = 1/(metric/UNT_GFP_mKa_factor_lower) # The inverse. Eq to UNT_GFP_mKa_factor_lower/metric   # This makes it so the threshold is always 1
			if loc_score_upper > loc_score_lower: # store the information on which threshold is being crossed, in case there is a mix of deloc and reloc.
												# The information is stored based on proximity and crossing the threshold so that it doesn't suddenly switch
				Upper  = 1
				Lower = 0
				return loc_score_upper, Upper, Lower
			if loc_score_lower > loc_score_upper:
				Upper = 0
				Lower = 1
				return loc_score_lower, Upper, Lower

		def reloc_score_mKO(row):
			metric = row[fluor_perc]
			loc_score_upper = metric/UNT_GFP_mKO_factor_upper
			loc_score_lower = 1/(metric/UNT_GFP_mKO_factor_lower) # The inverse. Eq to UNT_GFP_mKO_factor_lower/metric
			if loc_score_upper >= loc_score_lower:
				Upper  = 1
				Lower = 0
				return loc_score_upper, Upper, Lower
			if loc_score_lower >= loc_score_upper:
				Upper = 0
				Lower = 1
				return loc_score_lower, Upper, Lower

		def loc_score_to_binary(loc_score): # Conver the Loc_score value into a binary of whether the cell is currenlty dispalying relocalization
			if loc_score > 1:
				return(1)
			elif loc_score <=1:
				return(0)

		unif_mKa = Quant_prim.loc[Quant_prim["Myo1Identity"] == "Myo1_mKa"]
		unif_mKO = Quant_prim.loc[Quant_prim["Myo1Identity"] == "Myo1_mKO"]

		unif_mKa[["Loc_score", "Upper", "Lower"]] = unif_mKa.apply(reloc_score_mKa, axis = 1, result_type = "expand")

		unif_mKO[["Loc_score",  "Upper", "Lower"]] = unif_mKO.apply(reloc_score_mKO, axis = 1, result_type = 'expand')

		unif_mKa["Relocalized"] = pd.Series(unif_mKa["Loc_score"]).apply(loc_score_to_binary)

		unif_mKO["Relocalized"] = pd.Series(unif_mKO["Loc_score"]).apply(loc_score_to_binary)

		unif_mKO["mKO_factor_upper"]= UNT_GFP_mKO_factor_upper
		unif_mKO["mKO_factor_lower"]= UNT_GFP_mKO_factor_lower
		unif_mKO["mKO_RAWfactor_upper"] = Raw_factorUpper_mKO
		unif_mKO["mKO_RAWfactor_lower"] = Raw_factorLower_mKO

		unif_mKa["mKa_factor_upper"]= UNT_GFP_mKa_factor_upper
		unif_mKa["mKa_factor_lower"]= UNT_GFP_mKa_factor_lower
		unif_mKa["mKa_RAWfactor_upper"] = Raw_factorUpper_mKa
		unif_mKa["mKa_RAWfactor_lower"] = Raw_factorLower_mKa

		MMS_mKa[["Loc_score", "Upper", "Lower"]] = MMS_mKa.apply(reloc_score_mKa, axis = 1, result_type = "expand")
		MMS_mKO[["Loc_score", "Upper", "Lower"]] = MMS_mKO.apply(reloc_score_mKO, axis = 1, result_type = 'expand')

		MMS_mKa["Relocalized"] = pd.Series(MMS_mKa["Loc_score"]).apply(loc_score_to_binary)
		MMS_mKO["Relocalized"] = pd.Series(MMS_mKO["Loc_score"]).apply(loc_score_to_binary)

		MMS_mKa.to_csv(f"MMS_{Pos}mKa_wReloc.csv") #These are just the post_treatment observations
		MMS_mKO.to_csv(f"MMS_{Pos}mKO_wReloc.csv")

		# unif_mKa = pd.concat([UNT_mKa_all, MMS_mKa])
		# unif_mKO = pd.concat([UNT_mKO_all, MMS_mKO])

		CV_table_mKa = unif_mKa.groupby('Frames_post_treatment').Loc_score.agg([lambda x: variation(x), 'count'])
		CV_table_mKa.rename(columns  = {"lambda_x": "CV_sp", "count" : "cell_count"}, inplace = True) # 02/24/21 made this inplace
		unif_mKa = pd.merge(unif_mKa, CV_table_mKa, left_on= "Frames_post_treatment", right_index=True)

		CV_table_mKO  = unif_mKO.groupby('Frames_post_treatment').Loc_score.agg([lambda x: variation(x), 'count'])
		CV_table_mKO.rename(columns  = {"lambda_x": "CV_sp", "count" : "cell_count"})
		unif_mKO = pd.merge(unif_mKO, CV_table_mKa, left_on= "Frames_post_treatment", right_index=True)

		unif_mKa.to_csv(f"unif_mKa_{Pos}.csv")
		unif_mKO.to_csv(f"unif_mKO_{Pos}.csv")

		Movement_treat_course = pd.concat([unif_mKa, unif_mKO]) # was previously names Movement_pseudo
		Movement_treat_course.to_csv(f"Movement_treat_course_{Pos}.csv")

		return(f"The movement calculation and identification for position {Pos} is complete. Based on {percentile}th percentile")
	except:
		return(f"The movement calculation and identification for position {Pos} HAS FAILED")


pr = pn # pn # 8 #8 #pn - 1

# for p in range(1): #len(Quant_prim_index)):
	# I_move(p)

# #, Testing this new code to make sure that all cells presnet at treatment are
# Tracked_subset_pres_end = Tracked_subset[Tracked_subset["Cell_Barcode"].isin(Tracked_subset.iloc[-1,:]["Cell_Barcode"])]
# unif_mKa["Max_frame_pos"] = unif_mKa.groupby(["Unique_Pos"])["Frame_x"].transform('max') #* This is the max frame for a given position

# #!WhyTF is Unique_frame being ouput as an obbject and not a string?!
# df.groupby("Cell_Barcode")["Frame_x"].transform('max') #* This is the max frame for a given cell barcode
# df["Max_frame_pos"] = df.groupby(["Unique_Pos"])["Frame_x"].transform('max')


# def compare:
# 	if max(ser["Frame"]):
# 		return(0)
# 	if max(ser["Frame"]):
# 		return(1)


# df["Pres_end"] = df.apply(lambda x: x["Frame_x"] - x["Max_frame_pos"])

# df = df.loc[df["Pres_end"] == 1]
# dropped = df.loc[df["Pres_end"] != 1]

# dropped =
# df.drop(columns = "Max_frame_pos")

# df.groupby(["Cell_Barcode", "Unique_Frame"])




Parallel(n_jobs=pr, verbose = 100)(delayed(I_move)(p) for p in range(len(Quant_prim_index)))
#%%
# This is the postion combining code based on the column infromation.
os.chdir(post_path)

today = str(datetime.date.today())
year = str(datetime.date.today().year)
decade = year[:3]
#Write another indexing loop for movement

movement_unif_index = []
count = 0
for root, dirs, files, in os.walk(os.getcwd()):
	for name in files:
		if name.startswith("Move") and name.endswith(".csv") : # fix naming
			movement_unif_index.append({'Path': os.path.join(root, name)})
			count = count + 1
			print(count, end="\r")
		else:
			pass
	break #This makes the program run non recursively and not decend into daughter folders
movement_unif_index = pd.DataFrame(movement_unif_index)

def f_col(p):
	Col = int(p[p.find("p")+1:-4])
	return(Col)

def f_Position_move(z):
	start = z.find('course_')+7 #Note: The shift forward is dependant upon how you write out position
	end = z.find(".csv")   #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
	return(z[start:end])

def f_expdate(x):
	dend = x.find('r')
	expdate = x[0:dend] ### d****|r
	return(expdate)

def f_run(x):
	s = x.find("r")+1
	e = s+1
	r = x[s:e] ### d****r*|
	return(r)


movement_unif_index["Position"] = pd.Series(movement_unif_index["Path"]).apply(f_Position_move)
movement_unif_index["Run"] = pd.Series(movement_unif_index["Position"]).apply(f_run)
movement_unif_index["Col"] = pd.Series(movement_unif_index["Position"]).apply(f_col)
movement_unif_index["Date"] = pd.Series(movement_unif_index["Position"]).apply(f_expdate)
movement_unif_index.sort_values(by = ["Col"], inplace= True)  # the sorting is not required but makes it easier to track the progress of the jobs
movement_unif_index.reset_index(inplace= True, drop = True) #even though it is dropped when saving, leaving this here in case I do not want to read in later. Would be faster but don't want to change anything right now

movement_unif_index.to_csv("index_Movement.csv", index = False)

#%%
Col_list = Condition_information[["Date", "MapID (Col_Range)", "Run Number"]]
Col_list = Col_list.drop_duplicates().reset_index(drop = True) # This is the set of combination for the above information components
movement_unif_index = pd.read_csv("index_Movement.csv")
positions = movement_unif_index["Position"].unique()
#%%
os.chdir(post_path)
def f_Prot_combine_df(col_i): #This function runs based on the information file and searched for corresponding data. Done this way becasue it allows for parrallel joining without interferance
	try:
		col_inf = Col_list.iloc[(col_i), :] # col_inf here represents the condtion information search value NOT the name derived
		col_inf_minus = Col_list.iloc[(col_i - 1), :]
		matches = movement_unif_index[(movement_unif_index["Date"] == col_inf["Date"])  & (movement_unif_index["Run"] == col_inf["Run Number"]) & (col_inf_minus["Date"] == col_inf["Date"]) & (movement_unif_index["Col"] > col_inf_minus["MapID (Col_Range)"]) & (movement_unif_index["Col"] <= col_inf["MapID (Col_Range)"])]
		if len(matches) == 0: #Test to see the position being tested corresponds to the first column of the run being tested
			smallmap = Col_list[["Date", "Run Number", "MapID (Col_Range)"]]
			firstof_this_run = smallmap[(smallmap["Run Number"] == col_inf["Run Number"]) & (smallmap["Date"] == col_inf["Date"]) ].iloc[0]# This is a bad way of doing it but was quick to write
			matches = movement_unif_index[(movement_unif_index["Date"] == firstof_this_run["Date"]) & (movement_unif_index["Run"] == firstof_this_run["Run Number"]) & (movement_unif_index["Col"] <= firstof_this_run["MapID (Col_Range)"])] # Take the first row because it is a single ended comparison
		else:
			pass

		if len(matches) == 0: # If there are still no matched, it is assumed there is a mislabeling or the data has not being created yet
			return("No matches")
		else:
			Col = col_inf["MapID (Col_Range)"] #Return the col_range which the postional file falls within
			Date = col_inf["Date"] # Return the date which the postional file falls within
			Run = col_inf["Run Number"] #Return the run number which the postional file falls within
			concated = pd.DataFrame([])
			for pos in matches["Path"]: #Sequentially read in each matching postion data files and add it to the Column dataframe
				p = pd.read_csv(pos)
				concated = pd.concat([concated, p])


			CV_table= concated.groupby(['Protein', 'Frames_post_treatment']).Loc_score.agg([lambda x: variation(x), 'count']) #Calculated the varation and the count of cells for each protein in each postion file
			CV_table.rename(columns  = {"lambda_0": "CV_apos", "count" : "cell_count"}) # Should rename, but does not seem to be working
			concated = pd.merge(concated, CV_table, left_on= ['Protein', 'Frames_post_treatment'], right_index=True)
			concated.to_csv(f"Chamber_Col_{Date}_r{Run}_Ch{Col}_{today}.csv")
			return("Complete")
	except: #This is bad practice but should be fine for now
		return("Failed")

if len(movement_unif_index) > 20: #Arbitrarily set
	Parallel(n_jobs=pr, verbose = 100)(delayed(f_Prot_combine_df)(col_i) for col_i in range(len(Col_list))) #The parallel version should only be used when the number of positions being combined is large
else:
	Col_len = len(Col_list)
	for col_i in range(Col_len):
		res = f_Prot_combine_df(col_i)
		print(f"{col_i+1}/{Col_len}")
		print(res)

#%%
year = str(datetime.date.today().year)

## Create a list of combined column files
Chamber_index = []
count = 0
for root, dirs, files, in os.walk(os.getcwd()):
	for name in files:
		if name.startswith("Chamber") and name.endswith(f".csv"): # fix naming
			if name.endswith("index.csv"):
				pass
			else:
				Chamber_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
		else:
			pass
	break #This makes the program run non recursively and not decend into daughter folders
Chamber_index = pd.DataFrame(Chamber_index)
def f_chamber(col):
	s = col.find("_Col")+1
	e = col.find("_" + year)
	Chamber = col[s:e]
	return(Chamber)

Chamber_index["Chamber"] = pd.Series(Chamber_index.iloc[:,0]).apply(f_chamber)
Chamber_index.sort_values(by = "Chamber", inplace = True)
Chamber_index.to_csv("Chamber_index.csv")

#%%

#%%
# This is the percentage localization. Below, the percentages are calculated both including and excluding the pre-treatment values
try: Chamber_index
except NameError:
	os.chdir(post_path)
	Chamber_index = pd.read_csv("Chamber_index.csv")
today = str(datetime.date.today())
#%%

# This is a modification of Percent Trajectory concept which was calculated in Brandon't pipleine. Percent response is by nature a population metric but this is calculated and stored at the single cell level.
#The following series of functions produces the variabbles
#	1. Gobal percent repsonse for each protein (The number of cells which pass the threshold in total)
#	2. The time specific percentage of cells which are past the threshold
#	3. The percent trajectory of the max percentage response

def new_percentages_post_t(chamber, Chamber_df): # This function calculated the percentages with just the post-treatment values
	# Function passed the chamber number and the corresponding df from the master function
	#Subset the data to just the post treatment
	#There will be another and seperate percentage calculator that will include pre-treatment
	try:
		# Cammber_local = Chamber_df.sort_values(by = 'Frames_post_treatment', ascending= True) # Sort the values again, just in case they have gotten out of order
		Chamber_df_local = Chamber_df.copy()
		Chamber_df_local = Chamber_df_local[Chamber_df_local["Frames_post_treatment"] >= 0]

		Quant_unif_mKa = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKa']
		Quant_unif_mKO = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKO']


		frames_post_list= Quant_unif_mKa["Frames_post_treatment"].unique()
		# frame_index = Quant_unif_mKa["ImageID"].unique()
		percentage_reloc_p = pd.DataFrame([])

		try:
			mKa_prot = Quant_unif_mKa["Protein"].values[0]
		except:
			return("No mKa")
		try:
			mKO_prot = Quant_unif_mKO["Protein"].values[0]
		except:
			return("No mKO")
		for f in frames_post_list: # Loop through the list of post_treatment frame values
			t_mKa = Quant_unif_mKa[Quant_unif_mKa["Frames_post_treatment"] == f]
			t_mKO = Quant_unif_mKO[Quant_unif_mKO["Frames_post_treatment"] == f]

			reloc_t_mKa = t_mKa[t_mKa["Relocalized"] == 1]
			reloc_t_mKO = t_mKO[t_mKO["Relocalized"] == 1]

			percentage_reloc_t_mKa = (len(reloc_t_mKa)/len(t_mKa))*100
			percentage_reloc_t_mKO = (len(reloc_t_mKO)/len(t_mKO))*100

			t_less_reloc_mKa = t_mKa[t_mKa["Yet"] == 1]
			t_less_reloc_mKO = t_mKO[t_mKO["Yet"] == 1]

			percentage_moved_t_less_mKa = (len(t_less_reloc_mKa)/(len(t_mKa)))*100
			percentage_moved_t_less_mKO = ((len(t_less_reloc_mKO))/(len(t_mKO)))*100


			row = {
				"Frames_post_treatment" : [f],
				f"{mKa_prot}-perc_t" : [percentage_reloc_t_mKa],
				f"{mKO_prot}-perc_t" : [percentage_reloc_t_mKO],
				f"{mKa_prot}-perc_yet": [percentage_moved_t_less_mKa],
				f"{mKO_prot}-perc_yet": [percentage_moved_t_less_mKO]
			}
			f_row = pd.DataFrame(row)

			percentage_reloc_p = pd.concat([percentage_reloc_p, f_row])

		percentage_reloc_p.to_csv(f"percentage_reloc_pt_{chamber}.csv", index = False)
		return(None)
	except ZeroDivisionError:
		return("ZeroDivisionError")

def new_percentages_all_t(chamber, Chamber_df): # This function calculated the percentages with all the values. The only way that it differs is in that there is no subsetting for  frames after treatment
	# Function passed the chamber number and the corresponding df from the master function
	try:
		# Chamber_df_local = Chamber_df.sort_values(by = 'Frames_post_treatment', ascending= True) # Sort the values again just in case they have gotten out of order
		Chamber_df_local = Chamber_df.copy()
		Quant_unif_mKa = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKa']
		Quant_unif_mKO = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKO']


		frames_post_list= Quant_unif_mKa["Frames_post_treatment"].unique() #create a list of frames to test
		# frame_index = Quant_unif_mKa["ImageID"].unique()
		percentage_reloc_p = pd.DataFrame([])

		try:
			mKa_prot = Quant_unif_mKa["Protein"].values[0]
		except:
			return("No mKa")
		try:
			mKO_prot = Quant_unif_mKO["Protein"].values[0]
		except:
			return("No mKO")

		for f in frames_post_list: # Loop through the list of post_treatment frame values
			t_mKa = Quant_unif_mKa[Quant_unif_mKa["Frames_post_treatment"] == f]
			t_mKO = Quant_unif_mKO[Quant_unif_mKO["Frames_post_treatment"] == f]

			reloc_t_mKa = t_mKa[t_mKa["Relocalized"] == 1] #Create a subset of cells which are currently relocalized
			reloc_t_mKO = t_mKO[t_mKO["Relocalized"] == 1]

			percentage_reloc_t_mKa = (len(reloc_t_mKa)/len(t_mKa))*100 # Divide number of relocalized proteins in the current frame by the total number of cells in the current frame
			percentage_reloc_t_mKO = (len(reloc_t_mKO)/len(t_mKO))*100

			t_less_reloc_mKa = t_mKa[t_mKa["Yet"] == 1]
			t_less_reloc_mKO = t_mKO[t_mKO["Yet"] == 1]

			percentage_moved_t_less_mKa = (len(t_less_reloc_mKa)/(len(t_mKa)))*100
			percentage_moved_t_less_mKO = ((len(t_less_reloc_mKO))/(len(t_mKO)))*100


			row = {
				"Frames_post_treatment" : [f],
				f"{mKa_prot}-perc_t" : [percentage_reloc_t_mKa],
				f"{mKO_prot}-perc_t" : [percentage_reloc_t_mKO],
				f"{mKa_prot}-perc_yet": [percentage_moved_t_less_mKa],
				f"{mKO_prot}-perc_yet": [percentage_moved_t_less_mKO]
			}
			f_row = pd.DataFrame(row)

			percentage_reloc_p = pd.concat([percentage_reloc_p, f_row])

		percentage_reloc_p.to_csv(f"percentage_reloc_allt_{chamber}.csv", index = False)
		return(None)
	except ZeroDivisionError:
		return("ZeroDivisionError")

date_found = Chamber_index.iloc[0,0][-14:-4]

def percentages_bt_manager(chamber):
	try:
		Chamber_df = pd.read_csv(f"Chamber_{chamber}_{date_found}.csv")
	except FileNotFoundError:
		return(f"Col_merge file not found for {chamber}")

	Chamber_df.sort_values(by = 'Frames_post_treatment', ascending= True, inplace= True) # Sort the values by Frames_post_treatment to quantify the amount of relocalization so far
	#This is a new addtion to test whether the there has been relocalization yet

	# def complex_yet(x):
	# 	ind = x["Relocalized"].idxmax()
	# 	does = x.loc[ind]
	# 	x["yet"] = 0
	# 	x.loc[:ind, "yet"] = 0
	# 	x.loc[ind:, "yet"] = does
	# 	return

	def reloc_yet(x):
		ind = x.idxmax()
		does = x.loc[ind]
		x.loc[:ind] = 0
		x.loc[ind:] = does
		return(x)
	def workaround(ind):
		return(Chamber_df.loc[ind, "ImageID"])
	def does_workaround(ind):
		return(Chamber_df.loc[ind,"Relocalized"])

	Chamber_df["Yet"] = Chamber_df.groupby(by = "Cell_Barcode")["Relocalized"].transform(reloc_yet) # This repesents wether there has been relocaliztion yet
	Chamber_df["ind"] = Chamber_df.groupby(by = "Cell_Barcode")["Relocalized"].transform('idxmax')
	Chamber_df["Does"] = pd.Series(Chamber_df["ind"]).apply(does_workaround) #this will work for now but should make it come out of one of the other functions
	Chamber_df["When"] = pd.Series(Chamber_df["ind"]).apply(workaround)
	Chamber_df.drop(columns='ind', inplace = True)

	new_percentages_post_t(chamber, Chamber_df)
	new_percentages_all_t(chamber, Chamber_df)
	return(f"{chamber} Complete")

Parallel(n_jobs=pr, verbose = 100)(delayed(percentages_bt_manager)(p) for p in Chamber_index["Chamber"])


#%%
#############################################################################################

All_pos_allt_unified_prot_percentages = pd.DataFrame({
	"Frames_post_treatment": []
})
########NOTE## The output is All_pos_allt_unified_prot_percentages at ALL TIMEPOINTS
for c in Chamber_index["Chamber"]:
	try:
		percentage_reloc_p = pd.read_csv(f"percentage_reloc_allt_{c}.csv")
	except FileNotFoundError:
		print(f"Percentage fiile missing for {c}")

	All_pos_allt_unified_prot_percentages = pd.merge(All_pos_allt_unified_prot_percentages, percentage_reloc_p, how = 'outer', on='Frames_post_treatment')
	All_pos_allt_unified_prot_percentages["Time_post_treatment"] = All_pos_allt_unified_prot_percentages["Frames_post_treatment"] *time_per_frame

All_pos_allt_unified_prot_percentages.to_csv("All_pos_allt_percentage_sync.csv")
########NOTE## The output is All_pos_allt_unified_prot_percentages at ALL TIMEPOINTS


All_pos_pt_unified_prot_percentages = pd.DataFrame({
	"Frames_post_treatment": []
})


for c in Chamber_index["Chamber"]:
	try:
		percentage_reloc_p = pd.read_csv(f"percentage_reloc_pt_{c}.csv")
	except FileNotFoundError:
		print(f"Percentage fiile missing for {c}")

	All_pos_pt_unified_prot_percentages = pd.merge(All_pos_pt_unified_prot_percentages, percentage_reloc_p, how = 'outer', on='Frames_post_treatment')
	All_pos_pt_unified_prot_percentages["Time_post_treatment"] = All_pos_pt_unified_prot_percentages["Frames_post_treatment"] *time_per_frame
	# All_pos_pt_unified_prot_percentages.drop(columns='Frame', inplace = True)

# All_pos_pt_unified_prot_percentages = All_pos_pt_unified_prot_percentages.loc[:,~All_pos_pt_unified_prot_percentages.columns.duplicated()]

All_pos_pt_unified_prot_percentages.to_csv("All_pos_pt_percentage_sync.csv")

########################################################################################
#%%
list_prot_col_t = All_pos_pt_unified_prot_percentages.loc[:, (All_pos_pt_unified_prot_percentages.columns.str.endswith("-perc_t"))].columns

#Melt percentages for pt
All_pos_pt_t_percentages_melt = pd.melt(All_pos_pt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_t, var_name='Protein', value_name='Percentage_reloc')
All_pos_pt_t_percentages_melt["Time_post_treatment"] = All_pos_pt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
All_pos_pt_t_percentages_melt.to_csv("All_pos_pt_t_percentages_melt.csv")

#Melt percentages for allt
All_pos_allt_t_percentages_melt = pd.melt(All_pos_allt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_t, var_name='Protein', value_name='Percentage_reloc')
All_pos_allt_t_percentages_melt["Time_post_treatment"] = All_pos_allt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
All_pos_allt_t_percentages_melt.to_csv("All_pos_allt_t_percentages_melt")

#%%
list_prot_col_less = All_pos_pt_unified_prot_percentages.loc[:, (All_pos_pt_unified_prot_percentages.columns.str.endswith("-perc_yet"))].columns

#Melt percentages for pt_t less
All_pos_pt_t_less_percentages_melt = pd.melt(All_pos_pt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_less, var_name='Protein', value_name='Percentage_reloc_less')
All_pos_pt_t_less_percentages_melt["Time_post_treatment"] = All_pos_pt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
All_pos_pt_t_less_percentages_melt.to_csv("All_pos_pt_t_less_percentages_melt.csv")

#Melt percentages for all_t less
All_pos_allt_t_less_percentages_melt = pd.melt(All_pos_allt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_less, var_name='Protein', value_name='Percentage_reloc_less')
All_pos_allt_t_less_percentages_melt["Time_post_treatment"] = All_pos_allt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
All_pos_allt_t_less_percentages_melt.to_csv("All_pos_allt_t_less_percentages_melt.csv")

#%%
# Calculate the percent trajectory for all_t
def f_percent_trajectory(prot, perc): #This is set to run on the version which includes pre-treatment values. I determined that it is enlightening to know dynamics pre-treatment
	prot_subset = All_pos_allt_t_percentages_melt[(All_pos_allt_t_percentages_melt["Protein"]==prot)]
	prot_subset_pt = prot_subset[prot_subset["Frames_post_treatment"] >= 0]

	max_percent = np.max(prot_subset_pt["Percentage_reloc"])
	t_max = prot_subset[prot_subset["Percentage_reloc"] == max_percent]["Time_post_treatment"].values[0]

	percent_trajectory = (perc/max_percent) * 100
	return(max_percent, percent_trajectory, t_max)

applied_df= All_pos_allt_t_percentages_melt.apply(lambda row: f_percent_trajectory(row.Protein, row.Percentage_reloc), axis='columns', result_type='expand').rename(columns={0:"Max_percent", 1:"Percent_of_max_percent", 2:"t_max_percent"})
All_pos_t_percentages_melt= pd.concat([All_pos_allt_t_percentages_melt, applied_df], axis='columns')

All_pos_t_percentages_melt.to_csv("All_pos_allt_percentages_melt.csv", index = False)
#%%
#Calculate the trajectory for pt_t
def f_percent_trajectory(prot, perc): #This is set to run on the version which includes pre-treatment v alues. I determined that it is enlightening to know dynamics pre-treatment
	prot_subset = All_pos_pt_t_percentages_melt[(All_pos_pt_t_percentages_melt["Protein"]==prot)]
	prot_subset_pt = prot_subset[prot_subset["Frames_post_treatment"] >= 0]

	max_percent = np.max(prot_subset_pt["Percentage_reloc"])
	t_max = prot_subset[prot_subset["Percentage_reloc"] == max_percent]["Time_post_treatment"].values[0]

	percent_trajectory = (perc/max_percent) * 100 if max_percent != 0 else 0
	return(max_percent, percent_trajectory, t_max)

applied_df= All_pos_pt_t_percentages_melt.apply(lambda row: f_percent_trajectory(row.Protein, row.Percentage_reloc), axis='columns', result_type='expand').rename(columns={0:"Max_percent", 1:"Percent_of_max_percent", 2:"t_max_percent"})
All_pos_t_percentages_melt= pd.concat([All_pos_pt_t_percentages_melt, applied_df], axis='columns')


All_pos_t_percentages_melt.to_csv("All_pos_pt_percentages_melt.csv", index = False)
#%%
All_chambers_l = sorted(glob('Chamber_Col_*.csv'))
All_chambers_df = pd.concat((pd.read_csv(file) for file in All_chambers_l), ignore_index= True)

All_chambers_df.to_csv("ALL_CHAMBERS.csv")
#%%
# line = px.line(dataframe = )

#%%
import plotly.express as px

All_pos_t_percentages_melt["Protein"] = All_pos_t_percentages_melt["Protein"].astype('category')

Response_time_comparison = All_pos_t_percentages_melt[["Protein", "t_max_percent"]].drop_duplicates()
Response_time_comparison = Response_time_comparison.sort_values(by = 't_max_percent', ascending=True)

fig = px.bar(data_frame = Response_time_comparison,
		y= "Protein",
		x= 't_max_percent',
		orientation= 'h'
		).update_yaxes(categoryorder='total descending')
fig.update_xaxes(nticks=20)
fig.write_html("max_time.html")
fig
# Response_time_comparison
# sns.barplot(data = Response_time_comparison,
# 		x= "Protein",
# 		y= 't_max')

#%%
def remove(p):
	e = p.find("-")
	return(p[:e])

import plotly.express as px
Het_comparison = All_pos_t_percentages_melt[["Protein", "Max_percent"]].drop_duplicates()
Het_comparison = Het_comparison.sort_values(by = 'Max_percent')
Het_comparison["Protein"] = pd.Series(Het_comparison["Protein"]).apply(remove)
Het_comparison = Het_comparison[Het_comparison["Protein"] != 'FGV2']
Het_comparison["Sel"] = "reg"
Het_comparison.set_index("Protein", drop = False)
Het_comparison.loc[(Het_comparison.Protein == 'FLR1'),'Sel']='spec'
Het_comparison.loc[(Het_comparison.Protein == 'MRT4'),'Sel']='spec'

Percentage_ordered_graph = px.bar(data_frame= Het_comparison,
		y= "Protein",
		x= 'Max_percent',
		color = "Sel",
		orientation= 'h').update_yaxes(categoryorder='total descending')
Percentage_ordered_graph.update_xaxes(range=[0, 100])
Percentage_ordered_graph.update_layout(showlegend=False)
Percentage_ordered_graph.update_xaxes(nticks=20)
# sns.barplot(data = Het_comparison,
# 		x= "Protein",
# 		y= 'Max_percent')

Percentage_ordered_graph
Percentage_ordered_graph.write_html("percent_ordered{pt}.html")
#%%
t_less_penetrance = All_pos_pt_t_less_percentages_melt.groupby(by= "Protein").agg('max').reset_index()[["Protein", "Percentage_reloc_less"]]
t_less_penetrance["Protein"] = pd.Series(t_less_penetrance["Protein"]).apply(remove)
t_less_penetrance = t_less_penetrance[t_less_penetrance["Protein"] != 'FGV2'] ##### TEMP
t_less_penetrance["Sel"] = "reg"
t_less_penetrance.set_index("Protein", drop = False)

t_less_penetrance.loc[(t_less_penetrance.Protein == 'FLR1'),'Sel']='spec'
t_less_penetrance.loc[(t_less_penetrance.Protein == 'MRT4'),'Sel']='spec'

Pen_graph = px.bar(data_frame= t_less_penetrance,
		y= "Protein",
		x= 'Percentage_reloc_less',
		color = "Sel",
		orientation= 'h').update_yaxes(categoryorder='total descending')
Pen_graph.update_xaxes(range=[0, 100])
Pen_graph.update_layout(showlegend=False)
Pen_graph.write_html("Penetrance.html")
Pen_graph

#%%
# All_pos_t_percentages_melt
# t_traj_p = px.line(All_pos_t_percentages_melt, x = 'Frames_post_treatment', y = 'Percent_of_max_percent', color = "Protein")

# t_traj_p.write_html("t_traj_p.html")
# #%%

# All_pos_t_less_percentages_melt
# # t_less_traj = px.line(All_pos_t_less_percentages_melt, x = 'Frames_post_treatment', y = 'Percentage_reloc_less', color = "Protein")

# t_less_traj.write_html("t_less_traj.html")

# #%%

# #%%

# t_loc= px.line(All_chambers_df, x = 'Frames_post_treatment', y = 'Loc_score', color = "Protein")