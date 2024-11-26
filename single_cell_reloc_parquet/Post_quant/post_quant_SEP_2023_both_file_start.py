#%%
#, This document was modified from the AUG2 tree. This is a major modifiacation to make try and deal with the DivideByZero error which removes some cells
from datetime import datetime
from decimal import DivisionByZero
import math
import pandas as pd
import numpy as np
import os
import datetime
from scipy import stats
from joblib import Parallel, delayed
from glob import glob
from scipy.stats import variation
import single_cell_reloc_parquet.global_functions.global_variables as gv
from joblib import Parallel, delayed

#%% # Below was just copied from Quantificaiton.py becuase it was working at last check.
def Quantification_index_er():#? subset = False, subset_by = '', subset_collection = ''): Decided to just grab global variables instead of pasing local
	#! commented out because only Quant_ALL files are present
	Quantification_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith(".parquet") and name.startswith("Quantification_d"):
				Quantification_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			elif name.endswith(".csv") and name.startswith("Quantification_d"):
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

	def f_file_ext(x):
		extension = x[x.find(".")+1:]
		return(extension)

	# Quantification_index["Pad"] =
	Quantification_index["PositionID"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Position_ID)
	Quantification_index["Date"] = pd.Series(Quantification_index["PositionID"]).apply(f_expdate)
	Quantification_index["Frame"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Frame)
	Quantification_index['File_ext'] = pd.Series(Quantification_index.iloc[:,0]).apply(f_file_ext)

	#* NEWLY ADDED
	Quantification_index["Mod_epoch"] = pd.Series(Quantification_index["Path"]).apply(Mod_epoch)
	Quantification_index["Mod_date"] = pd.Series(Quantification_index["Mod_epoch"]).apply(epoch_convert)

	Quantification_index["Max_epoch_frame"] = Quantification_index.groupby("Frame")["Mod_epoch"].transform('max') #? This isn't sylish but works
	Quantification_index = Quantification_index.loc[Quantification_index["Mod_epoch"] == Quantification_index["Max_epoch_frame"]]
	Quantification_index.drop(columns=["Max_epoch_frame"])

	if Global_variables['subset'] == True:
		Quantification_index = Quantification_index[Quantification_index[Global_variables['subset_by']].isin(Global_variables['subset_collection'])]
		print(f"Post-quant will be performed for positions which correspond to the following {Global_variables['subset_by']}: {Global_variables['subset_collection']}")
	else:
		pass
	Quantification_index.sort_values(by = "PositionID", inplace = True)
	Quantification_index.to_parquet("Quantification_index.parquet")
	return(Quantification_index)

#%% Combine frames into posision
def combine_pos(pos):
	Quant_frame_comb = pd.DataFrame([])
	subset = Quantification_index.loc[Quantification_index["PositionID"] == pos]
	subset_s = subset.sort_values(by = "Frame")
	for qf in range(len(subset)):
		try:
			q = pd.read_parquet(subset_s.iloc[qf,0])
		except:
			try:
				q = pd.read_csv(subset_s.iloc[qf,0])
			except:
				pass
		Quant_frame_comb = pd.concat([Quant_frame_comb, q])
	Quant_frame_comb.to_parquet(f"Quant_{pos}_ALL.parquet")
	return(f"{pos} frame merging is complete")

# for pos in positions:
# 	Quant_frame_comb = pd.DataFrame([])
# 	subset = Quantification_index.loc[Quantification_index["PositionID"] == pos]
# 	subset.sort_values(by = "Frame", inplace = True)
# 	for qf in range(len(subset)):
# 		q = pd.read_parquet(subset.iloc[qf,0])
# 		Quant_frame_comb = pd.concat([Quant_frame_comb, q])
# 	Quant_frame_comb.to_parquet(f"Quant_{pos}_ALL.parquet")
# 	print(f"{pos} merging is complete")

#%% Generate an index of Quant_ALL
# # Create an index of Quantification dfs
############################################

def Quant_ALL_index_er():
	year = str(datetime.datetime.now().year) #this assumes the analysis is completed within the same decade of starting. A reasonable expectation, but something to check if running aroudn the new year
	decade = year[:3]
	del year

	def f_Position_ID_qALLi(z): #* This was a stupid thing to leave in the code
		start = z.find('ant_')+4
		end = z.find("_ALL")
		return(z[start:end])

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

	def f_file_ext(x):
		extension = x[x.find(".")+1:]
		return(extension)


	Quant_ALL_index  = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:

			if name.endswith("_ALL.parquet") and name.startswith("Quant"): # fix naming
				Quant_ALL_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			elif name.endswith("_ALL.csv") and name.startswith("Quant"): # fix naming
				Quant_ALL_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		break #This makes the program run non-recursively and not decend into daughter folders

	Quant_ALL_index = pd.DataFrame(Quant_ALL_index)
	Quant_ALL_index["PositionID"] = pd.Series(Quant_ALL_index.iloc[:,0]).apply(f_Position_ID_qALLi)
	Quant_ALL_index['File_ext'] = pd.Series(Quant_ALL_index.iloc[:,0]).apply(f_file_ext)

	#* NEWEST
	Quant_ALL_index["Mod_epoch"] = pd.Series(Quant_ALL_index["Path"]).apply(Mod_epoch)
	Quant_ALL_index["Mod_date"] = pd.Series(Quant_ALL_index["Mod_epoch"]).apply(epoch_convert)

	Quant_ALL_index["Max_epoch_pos"] = Quant_ALL_index.groupby("PositionID")["Mod_epoch"].transform('max') #? This isn't sylish but works
	Quant_ALL_index = Quant_ALL_index.loc[Quant_ALL_index["Mod_epoch"] == Quant_ALL_index["Max_epoch_pos"]]
	Quant_ALL_index.drop(columns=["Max_epoch_pos"])

	Quant_ALL_index.sort_values(by = "PositionID", inplace = True)
	Quant_ALL_index.to_parquet("Quant_ALL_index.parquet")
	return(Quant_ALL_index)

#%%
# Cell_index_timing = info_simple[["track_index", "track_start_frame"]]
# os.chdir(microfluidics_results)
# Quant_ALL_index = pd.read_parquet("Quant_ALL_index.parquet") #*#* Does the file need to read in to work?
# Quant_ALL_index.drop(Quant_ALL_index.columns[Quant_ALL_index.columns.str.contains('Unnamed',case = False)],axis = 1, inplace= True)
#%%

def f_convert_date(d):
	d = str(d)
	s = d.find(decade)+5
	return("d" + d[s:])

# Condition_information["Time"] = Condition_information["Time.1"] # This is required if there is a "Time" for aquisition, and timing for treatment
# Condition_information = Condition_information[["Date", "Run Number", "Time", "MapID (Col_Range)", "Myo1Marker", "Protein", "Time"]]
#%%
## Below is a far more complex implementation of Myo1 determination
#Instead of being just based on the max and min of each channel, this version also takes into account the start frame and also the dynamic range

def Strain_ID_multiplex(p, multiplex = True):
	os.chdir(post_path)
	try:
		open_path = Quant_ALL_index.iloc[p,0]
		#* I should have already filtered for the most recent version in Quant_ALL
		if open_path.endswith('.parquet'):
			Quant_ALL = pd.read_parquet(Quant_ALL_index.iloc[p,0])
		elif open_path.endswith('.csv'):
			Quant_ALL = pd.read_csv(Quant_ALL_index.iloc[p,0])

		Quant_ALL.drop(Quant_ALL.columns[Quant_ALL.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
		print(p)

		frame_max= np.max(Quant_ALL["Frame"])

		try: trackedCells = np.unique(Quant_ALL['Cell_Barcode'])
		except KeyError:
			return(f"Fail on index {p} at trackedCells Stage")

		if len(trackedCells) < 50:       ### Not sure if it would be a good idea to make sure that the position has enough cells. This does not gauruntee enough cells in future factor determination.
			return("There are less than 50 cells in time series")

		# def f_start_frame_test(cell):
		# 	start_frame = info_simple.loc[cell]["track_start_frame"].values[0]
		# 	return(start_frame)

		#* Modified what was originally written by Brandon

		Myo1_info= pd.DataFrame([])
		for i in trackedCells: #* This is a very easy line to miss! The function is bassically just a position specific manager for making cell measurements so that more calculations can be done in parrallel with the proper rate of turnover
			trackSubset = Quant_ALL.loc[Quant_ALL['Cell_Barcode'] == i]

			#only accept cells that have been tracked in the experiment for more than 15 frames (105 minutes)
			# if len(trackSubset) > 15:  ######## The was a check to make sure that cells were present for more than 15 frams. This is would throw out new cells born within the last 105 minutes
			validTrackedCells = i
			# else:
				# continue

			range_direct = 0 #* This line was not in use

			if len(trackSubset[(trackSubset["Progen_bud"] == 1) & (trackSubset["Frame"] != 1)])>0:
				mKAstart_bud_actv = trackSubset[trackSubset["Progen_bud"] ==1]['x99thPercentile_Diff_background_mKate'].values[0]
				mKOstart_bud_actv = trackSubset[trackSubset["Progen_bud"] == 1]['x99thPercentile_Diff_background_mKO'].values[0]
			else:
				try:
					mKAstart_bud_actv = trackSubset[trackSubset["Frame"] == frame_max]['x99thPercentile_Diff_background_mKate'].values[0]
					mKOstart_bud_actv = trackSubset[trackSubset["Frame"] == frame_max]['x99thPercentile_Diff_background_mKO'].values[0]
				except IndexError:
					mKAstart_bud_actv = 0
					mKOstart_bud_actv = 0

			mKOmax = trackSubset.iloc[trackSubset['x99thPercentile_Diff_background_mKO'].argmax(),:]
			mKOmin = trackSubset.iloc[trackSubset['x99thPercentile_Diff_background_mKO'].argmin(),:]

			if mKOmax['Frame'] > mKOmin['Frame']:
				KOdirection = 'Increase'
			else:
				KOdirection = 'Decrease'

			if mKOmin['x99thPercentile_Diff_background_mKO'] >0:
				mKOsignalChange = mKOmax['x99thPercentile_Diff_background_mKO']/mKOmin['x99thPercentile_Diff_background_mKO'] #* are some min values of the difference too close to zero?????
			else: #? Not sure which of the options below should be used
				# mKOsignalChange = mKOmax['x99thPercentile_Diff_background_mKO']/mKOmin['factor_mKO_background_Avg'] #* This may work but it is not as accurate as the one below
				mKOsignalChange = (mKOmax['x99thPercentile_Diff_background_mKO'] + mKOmin['factor_mKO_background_Avg'])/mKOmin['factor_mKO_background_Avg']

			mKOsignalChange = abs(mKOsignalChange)#/64.05)

			mKAmax = trackSubset.iloc[trackSubset['x99thPercentile_Diff_background_mKate'].argmax(),:]
			mKAmin = trackSubset.iloc[trackSubset['x99thPercentile_Diff_background_mKate'].argmin(),:]
			if mKAmax['Frame'] > mKAmin['Frame']: # for now this is just just a what is greater. Must modifiy to a threshold
				KAdirection = 'Increase'
			else:
				KAdirection = 'Decrease'

			if mKAmin['x99thPercentile_Diff_background_mKate'] > 0:
				mKAsignalChange = mKAmax['x99thPercentile_Diff_background_mKate']/mKAmin['x99thPercentile_Diff_background_mKate']
			else: #? Not sure which of the options below should be used. It seems as though the latter is better
				# mKAsignalChange = mKAmax['x99thPercentile_Diff_background_mKate']/mKAmin['factor_mKate_background_Avg'] #* This is not as accurate as below
				mKAsignalChange = (mKAmax['x99thPercentile_Diff_background_mKate'] + mKAmin['factor_mKate_background_Avg'])/mKAmin['factor_mKate_background_Avg']

			mKAsignalChange = abs(mKAsignalChange)#/25)


			# if range_direct == 0: #* This line was not in use
			#, With a 25% margin of confidence over alternative, determine the bud-neck type - ie. the strian for co-imaged strains. Branched decision
			if mKAstart_bud_actv > 1.25*mKOstart_bud_actv: #This was just changed back to comparison, but threshold of 2000 was being used before ## 14-12-21 Added a 25% confidence over the lower value
				Myo1ID = 'Myo1_mKa'
				foldChangeKO = mKOsignalChange
				foldChangeKA = mKAsignalChange
				boolRange = 0
				boolProgen = 1

			elif mKOstart_bud_actv> 1.25*mKAstart_bud_actv:
				Myo1ID = 'Myo1_mKO'
				foldChangeKO = mKOsignalChange
				foldChangeKA = mKAsignalChange
				boolRange = 0
				boolProgen = 1
			else:
				# Now determine whether the given tracked cells is mKO of mKa Myo1. This is only for now as this can be determined in TracX
				if mKAsignalChange > 1.25*mKOsignalChange: #This was just changed back to comparison, but threshold of 2000 was being used before
					Myo1ID = 'Myo1_mKa'
					foldChangeKO = mKOsignalChange
					foldChangeKA = mKAsignalChange
					boolRange = 1
					boolProgen = 0

				elif mKOsignalChange > 1.25*mKAsignalChange:
					Myo1ID = 'Myo1_mKO'
					foldChangeKO = mKOsignalChange
					foldChangeKA = mKAsignalChange
					boolRange = 1
					boolProgen = 0

				else:
					Myo1ID = 'Cannot Be Determined'
					foldChangeKO = mKOsignalChange
					foldChangeKA = mKAsignalChange
					boolRange = 0
					boolProgen = 0

			#* Genarate a row of measurements to be concatenated onto the growing dataframe
			measurements = {
				"TrackID_valid" : [validTrackedCells],
				"Myo1Identity" : [Myo1ID],
				"mKO_foldChange" : [foldChangeKO],
				"mKO_direction" : [KOdirection],
				"mKA_foldChange" : [foldChangeKA],
				"mKa_direction" : [KAdirection],
				"byProgen_bud": [boolProgen],
				"byRange": [boolRange]
			}

			measurements = pd.DataFrame(measurements)
			Myo1_info = pd.concat([Myo1_info, measurements])

		try: Quant_FIN_primary = pd.merge(Quant_ALL, Myo1_info, left_on="Cell_Barcode", right_on = "TrackID_valid") #* This is an inner join. The indices should have same values as Myo1_info is derived from Quant_ALL.
		except KeyError:
			return(f"index {p} of Quant_ALL_index may be too short. Cell was not found in >15 frames. Pos skipped")

		Pos = Quant_ALL["Cell_Barcode"].values[0] # This is a temporary and false value, based on the first cell barcode as all within df will have the same position
		Pos = Pos[0:Pos.find("c")]

		s = Pos.find("r")
		Run_n = int(Pos[s+1 : s+2])

		Run_info = Condition_information[(Condition_information["Date"] == Quant_ALL.iloc[0,:]["Date"]) & (Condition_information["Run Number"] == Run_n)]

		# DateCond = Run_info["Date"].values[0]
		Time_treat = Run_info["Time"].values[0]

		Col = int(Pos[Pos.find("p")+1:-4])  # This is the column of the postion which is being tested
		Pos_info = Run_info[Run_info["MapID (Col_Range)"]>=Col].reset_index (drop=True) # Col HERE represents the NAME and NOT the information. eg. info column 10,20,30 etc.  >= pos_col 10
		Pos_info = Pos_info.iloc[0:2] # THIS IS A VERY IMPORTANT LINE
		Prot_strain1 = Pos_info[Pos_info["Myo1Marker"] == "mKO"]["Protein"].values[0]
		Prot_strain2 = Pos_info[Pos_info["Myo1Marker"] == "mKa"]["Protein"].values[0]
		Col_info = Pos_info["MapID (Col_Range)"].values[0]

		def Protein_label_multi(myo1): #, The purpose of this funtion is to associate the determinined Myo1 major fluorescence to the apppropriate protien
			myo1 = myo1.values[0] # This is for the barcode grouped version. It should run faster
			if myo1 == "Myo1_mKO":
				return(Prot_strain1)
			elif myo1 == "Myo1_mKa":
				return(Prot_strain2)
			else:
				return(myo1)

		def Is_treated(frame): #, The purpose of this function is to label each frame by the relative time from treatment in frames. These values can later be converted to time values by multiplying by the appropriate scaling value.
			frame = frame.values[0] # Again, this is for the barcode grouped version
			f_treated = Time_treat/ time_per_frame
			if frame < f_treated:
				return (0)
			if frame >= f_treated:
				return(1)

		Quant_FIN_primary["Col_info"] = int(Col_info)
		# Quant_FIN_primary["Protein"] = pd.Series(Quant_FIN_primary["Myo1Identity"]).apply(Protein_label_multi)
		Quant_FIN_primary["Protein"] = Quant_FIN_primary.groupby(by = ["Cell_Barcode"])["Myo1Identity"].transform(Protein_label_multi) # This is the grouped run. This should be faster because it is only applyong once per cell

		# Quant_FIN_primary["Protein"] = pd.Series
		# (Quant_FIN_primary["Myo1Identity"]).apply(Protein_label_multi) #fix the lookup table

		# Quant_FIN_primary["Is_treated"] = pd.Series(Quant_FIN_primary["Frame"]).apply(Is_treated)
		Quant_FIN_primary["Is_treated"] = Quant_FIN_primary.groupby(by = ["ImageID"])["Frame"].transform(Is_treated)

		def f_Frames_post_treatment_shift(i, i_pt):
			v = i - i_pt
			return(v)

		psuedo_map = Quant_FIN_primary[["ImageID", "Is_treated", "Frame"]].drop_duplicates() # Grab just the imageID, "Is_treated", and "Frame"
		psuedo_map.sort_values(by = ["Frame"], inplace=True)   ####  Add sorting to make sure that the index comparison will be correct

		psuedo_map.reset_index(inplace=True, drop = True) #Reindex to make sure that the index numbers are correct and can be used for comparison
		i_pt = psuedo_map[psuedo_map["Is_treated"] == 1].index[0] - 1 # This isn't the most pretty way to go about but works for now
		psuedo_map["Frames_post_treatment"] = pd.Series(psuedo_map.index).apply(lambda x: f_Frames_post_treatment_shift(x, i_pt))
		psuedo_map.drop(columns=["Frame", "Is_treated"], inplace= True) #* This corrects the Frame_x issue issue the existed in prior versions
		Quant_FIN_primary = pd.merge(Quant_FIN_primary, psuedo_map, how = "left", on='ImageID')
		Quant_FIN_primary.to_parquet(f"Quant_{Pos}_primary.parquet")
	except IndexError:
		return(f"IndexError on {p}")
	except ValueError:
		return(f"ValueError on {p}")
# def Protein_label_multi(myo1):
# 	if myo1 == "Myo1_mKO":
# 		return(BHY131)
# 	if myo1 == "Myo1_mKa":
# 		return(BHY175)

# def Is_treated(frame):
# 	f_treated = Time_treat/ time_per_frame
# 	if frame < f_treated:
# 		return (0)
# 	if frame >= f_treated:
# 		return(1)




# for p in range(len(Quant_ALL_index)):
# 	Strain_ID_multiplex(p)


# Strain_ID_multiplex(0)


#%%
# def f_Position_ID_qALLi(z): #! Commmented out SEP10
# 	start = z.find('ant_')+4
# 	end = z.find("_prim")
# 	return(z[start:end])

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
def Quant_prim_index_er():
	Quant_prim_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith("mary.parquet") and name.startswith("Quant"): # fix naming
				Quant_prim_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		break #This makes the program run non recursively and not decend into daughter folders
	def f_Position_ID_qALLi(z):
		start = z.find('ant_')+4
		end = z.find("_prim")
		return(z[start:end])


	Quant_prim_index = pd.DataFrame(Quant_prim_index)
	Quant_prim_index["PositionID"] = pd.Series(Quant_prim_index.iloc[:,0]).apply(f_Position_ID_qALLi)
	# Quant_prim_index["Position"] = pd.Series(Quant_prim_index["Path"]).apply(f_Position_ID_qALLi)
	# Quant_prim_index["Run"] = pd.Series(Quant_prim_index["Position"]).apply(f_run)
	# Quant_prim_index["Col"] = pd.Series(Quant_prim_index["Position"]).apply(f_col)
	# Quant_prim_index["Date"] = pd.Series(Quant_prim_index["Position"]).apply(f_expdate)

	Quant_prim_index.sort_values(by = "PositionID", inplace = True)
	Quant_prim_index.to_parquet("Quant_prim_index.parquet")#? index = False)
	Quant_prim_index.to_csv("Quant_prim_index.csv")	#. Temp also save as csv to see if this is different
	return(Quant_prim_index)

##### MOVED THE FRAME JUSTIFICATION UPSTREAM
#%%
####### THIS IS ALL CELLS USED. THE ONLY VERSION THAT DOES NOT IS DEC16_SUBSET
####### This is the two sides Locp = 1_score test
# percentile = str(input("What is the the parameter to check Loc Change on?"))







#%%
def MI_move(p, percentile, perc_path): #* The Quant all INdex
	try: #TODO: Rename Frame_x to just Frame
		#* The line directly below reads in just the required information to avoid another slicing of a cumbersome df
		Quant_prim = pd.read_parquet(Quant_prim_index.iloc[p,0])#. , usecols = ['Cell_Barcode', 'ImageID', 'Date', 'Frame', 'Unique_Frame', 'factor_median_OBJ_GFP', 'factor_mean_OBJ_GFP', 'x90thPercentile_norm_OBJ_Median_GFP', 'x95thPercentile_norm_OBJ_Median_GFP', 'x99thPercentile_norm_OBJ_Median_GFP', 'Progen_bud', 'x90thPercentile_GFP_RAW', 'x95thPercentile_GFP_RAW', 'x99thPercentile_GFP_RAW', 'max_GFP_RAW', 'max_mKa_RAW', 'max_mKO_RAW', 'TrackID_valid', 'Myo1Identity', 'mKO_foldChange', 'mKO_direction', 'mKA_foldChange', 'mKa_direction', 'byProgen_bud', 'byRange', 'Col_info', 'Protein', 'Is_treated', 'Frames_post_treatment']) #* From here on are newly added columns to determine the cell stage. Could do a pd.merge at the end instead.
		# "x80thPercentile_Diff_background_mKate", "x90thPercentile_Diff_background_mKate", "x99thPercentile_Diff_background_mKate", "x80thPercentile_Diff_background_mKO", "x90thPercentile_Diff_background_mKO", "x99thPercentile_Diff_background_mKO", "averageIntensity_mKO_Frame", "averageIntesntiy_mKO_Background", "averageIntensity_mKO_Object", "mKO_spread", "averageIntensity_mKate_Frame", "averageIntesntiy_mKate_Background", "averageIntensity_mKate_Object", "mKate_spread", "factor_median _OBJ_KO", "factor_mean_OBJ_KO", "factor_total_OBJ_KO", "factor_mKO_background_Med", "factor_mKO_background_Avg", "factor_mKO_background_Tot", "factor_median_OBJ_mKate", "factor_mean_OBJ_mKate", "factor_total_OBJ_mKate", "factor_mKate_background_Med", "factor_mKate_background_Avg", "factor_mKate_background_Tot", "x60thPercentile_mKa_RAW", "x80thPercentile_mKa_RAW", "x90thPercentile_mKa_RAW", "x95thPercentile_mKa_RAW", "x99thPercentile_mKa_RAW", "max_mKa_RAW", "x60thPercentile_mKO_RAW", "x80thPercentile_mKO_RAW", "x90thPercentile_mKO_RAW", "x95thPercentile_mKO_RAW", "x99thPercentile_mKO_RAW", "max_mKO_RAW"])

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
		Quant_f__twe = Quant_prim.loc[(Quant_prim["Frame"] <= 20) & (Quant_prim["Is_treated"] == 0)].copy()
		Median_pop = Quant_f__twe.groupby(by = "Cell_Barcode").agg('min').groupby(by = "Myo1Identity").agg('median', numeric_only = True)
		Std_pop = Quant_f__twe.groupby(by = "Cell_Barcode").agg('min').groupby(by = "Myo1Identity").agg('std', numeric_only = True)

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
		Quant_no.to_csv(f"{perc_path}/Dropped_dead_{Pos}.csv")

		Quant_prim.set_index(["Cell_Barcode", "ImageID"], inplace = True) #! Changed index

		fluor_perc = (f"x{percentile}thPercentile_norm_OBJ_Median_GFP")
		fluor_perc_r = (f"x{percentile}thPercentile_GFP_RAW")

		UNT_mKa_all = Quant_prim[(Quant_prim["Myo1Identity"] == "Myo1_mKa") & (Quant_prim["Is_treated"] == 0)]
		UNT_mKa = UNT_mKa_all[UNT_mKa_all["Frames_post_treatment"] >= -10] # This is the subset of pre-treatment, 3 frames before treatment. May be able to expand but not clear that it is needed #! TESTING THE EFFECT OF USING FRAMES MORE THAN -10 INSTEAD OF -3
		#TODO: Decide finally if I am using 10 frames (75 minutes or 3 frames (29.5 minues))

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

		MMS_mKa = Quant_prim.loc[(Quant_prim["Myo1Identity"] == "Myo1_mKa") & (Quant_prim["Is_treated"] == 1)].copy()
		# MMS_mKa.loc[:,"mKa_factor_upper"] = UNT_GFP_mKa_factor_upper #. This is the old format, but might be causing slice error later on
		MMS_mKa["mKa_factor_upper"] = UNT_GFP_mKa_factor_upper
		MMS_mKa["mKa_factor_lower"] = UNT_GFP_mKa_factor_lower

		MMS_mKO = Quant_prim.loc[(Quant_prim["Myo1Identity"] == "Myo1_mKO") & (Quant_prim["Is_treated"] == 1)].copy()
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

		UNT_mKa.to_parquet(f"{perc_path}/UNT_{Pos}mKa.parquet") # Output the cells and information on which parameters are derived
		UNT_mKO.to_parquet(f"{perc_path}/UNT_{Pos}mKO.parquet")

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

		unif_mKa = Quant_prim.loc[Quant_prim["Myo1Identity"] == "Myo1_mKa"].copy()
		unif_mKO = Quant_prim.loc[Quant_prim["Myo1Identity"] == "Myo1_mKO"].copy()

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

		MMS_mKa[["Loc_score", "Upper", "Lower"]] = MMS_mKa.apply(reloc_score_mKa, axis = 1, result_type = "expand") #. This was also in .loc form, but changed back to this
		MMS_mKO[["Loc_score", "Upper", "Lower"]] = MMS_mKO.apply(reloc_score_mKO, axis = 1, result_type = 'expand')

		MMS_mKa["Relocalized"] = pd.Series(MMS_mKa["Loc_score"]).apply(loc_score_to_binary)
		MMS_mKO["Relocalized"] = pd.Series(MMS_mKO["Loc_score"]).apply(loc_score_to_binary)

		MMS_mKa.to_parquet(f"{perc_path}/MMS_{Pos}mKa_wReloc.parquet") #These are just the post_treatment observations
		MMS_mKO.to_parquet(f"{perc_path}/MMS_{Pos}mKO_wReloc.parquet")

		# unif_mKa = pd.concat([UNT_mKa_all, MMS_mKa])
		# unif_mKO = pd.concat([UNT_mKO_all, MMS_mKO])

		CV_table_mKa = unif_mKa.groupby('Frames_post_treatment').Loc_score.agg([lambda x: variation(x), 'count'])
		CV_table_mKa.rename(columns  = {"lambda_x": "CV_sp", "count" : "cell_count"}, inplace = True) # 02/24/21 made this inplace
		unif_mKa = pd.merge(unif_mKa, CV_table_mKa, left_on= "Frames_post_treatment", right_index=True)

		CV_table_mKO  = unif_mKO.groupby('Frames_post_treatment').Loc_score.agg([lambda x: variation(x), 'count'])
		CV_table_mKO.rename(columns  = {"lambda_x": "CV_sp", "count" : "cell_count"})
		unif_mKO = pd.merge(unif_mKO, CV_table_mKa, left_on= "Frames_post_treatment", right_index=True)

		unif_mKa.to_parquet(f"{perc_path}/unif_mKa_{Pos}.parquet")
		unif_mKO.to_parquet(f"{perc_path}/unif_mKO_{Pos}.parquet")

		Movement_treat_course = pd.concat([unif_mKa, unif_mKO]) # was previously names Movement_pseudo
		Movement_treat_course.to_parquet(f"{perc_path}/Movement_treat_course_{Pos}.parquet")

		return(f"The movement calculation and identification for position {Pos} is complete. Based on {percentile}th percentile. Saved in {perc_path}")
	except:
		return(f"The movement calculation and identification for position {Pos} HAS FAILED. NOT SAVED IN {perc_path}")

#%%
# for p in range(1): #len(Quant_prim_index)):
	# MI_move(p)

# #, Testing this new code to make sure that all cells presnet at treatment are
# Tracked_subset_pres_end = Tracked_subset[Tracked_subset["Cell_Barcode"].isin(Tracked_subset.iloc[-1,:]["Cell_Barcode"])]
# unif_mKa["Max_frame_pos"] = unif_mKa.groupby(["Unique_Pos"])["Frame"].transform('max') #* This is the max frame for a given position

# #!WhyTF is Unique_frame being ouput as an obbject and not a string?!
# df.groupby("Cell_Barcode")["Frame"].transform('max') #* This is the max frame for a given cell barcode
# df["Max_frame_pos"] = df.groupby(["Unique_Pos"])["Frame"].transform('max')


# def compare:
# 	if max(ser["Frame"]):
# 		return(0)
# 	if max(ser["Frame"]):
# 		return(1)


# df["Pres_end"] = df.apply(lambda x: x["Frame"] - x["Max_frame_pos"])

# df = df.loc[df["Pres_end"] == 1]
# dropped = df.loc[df["Pres_end"] != 1]

# dropped =
# df.drop(columns = "Max_frame_pos")

# df.groupby(["Cell_Barcode", "Unique_Frame"])




#%%

def movement_unif_index_er(perc_path):
	os.chdir(perc_path)
	movement_unif_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.startswith("Move") and name.endswith(".parquet") : # fix naming
				movement_unif_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		break #This makes the program run non recursively and not decend into daughter folders
	movement_unif_index = pd.DataFrame(movement_unif_index)

	def f_col(p):
		Col = int(p[p.find("p")+1:-4])
		r = int(Col)
		return(Col)

	def f_Position_move(z):
		start = z.find('course_')+7 #Note: The shift forward is dependant upon how you write out position
		end = z.find(".parquet")   #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
		return(z[start:end])

	def f_expdate(x):
		dend = x.find('r')
		expdate = x[0:dend] ### d****|r
		return(expdate)

	def f_run(x):
		s = x.find("r")+1
		e = s+1
		r = x[s:e] ### d****r*|
		r = int(r)
		return(r)


	movement_unif_index["Position"] = pd.Series(movement_unif_index["Path"]).apply(f_Position_move)
	movement_unif_index["Run"] = pd.Series(movement_unif_index["Position"]).apply(f_run)
	movement_unif_index["Col"] = pd.Series(movement_unif_index["Position"]).apply(f_col)
	movement_unif_index["Date"] = pd.Series(movement_unif_index["Position"]).apply(f_expdate)
	movement_unif_index.sort_values(by = ["Col"], inplace= True)  # the sorting is not required but makes it easier to track the progress of the jobs
	movement_unif_index.reset_index(inplace= True, drop = True) #even though it is dropped when saving, leaving this here in case I do not want to read in later. Would be faster but don't want to change anything right now
	movement_unif_index.to_parquet("index_Movement.parquet", index = False)
	return(movement_unif_index)

#%%
def f_Prot_combine_df(col_i, perc_path): #This function runs based on the information file and searched for corresponding data. Done this way becasue it allows for parrallel joining without interferance
	try:
		col_inf = Col_list.iloc[(col_i), :] # col_inf here represents the condtion information search value NOT the name derived
		col_inf_minus = Col_list.iloc[(col_i - 1), :]
		matches = movement_unif_index.loc[(movement_unif_index["Date"] == col_inf["Date"])  & (movement_unif_index["Run"] == col_inf["Run Number"]) & (col_inf_minus["Date"] == col_inf["Date"]) & (movement_unif_index["Col"] > col_inf_minus["MapID (Col_Range)"]) & (movement_unif_index["Col"] <= col_inf["MapID (Col_Range)"])] #* Changed to .loc while searching for mismatch issue. Should not make a difference but will be more stable
		if len(matches) == 0: #Test to see the position being tested corresponds to the first column of the run being tested
			smallmap = Col_list[["Date", "Run Number", "MapID (Col_Range)"]]
			firstof_this_run = smallmap[(smallmap["Run Number"] == col_inf["Run Number"]) & (smallmap["Date"] == col_inf["Date"]) ].iloc[0]# This is a bad way of doing it but was quick to write
			matches = movement_unif_index[(movement_unif_index["Date"] == firstof_this_run["Date"]) & (movement_unif_index["Run"].astype(int) == int(firstof_this_run["Run Number"])) & (movement_unif_index["Col"].astype(int) <= int(firstof_this_run["MapID (Col_Range)"]))] # Take the first row because it is a single ended comparison
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
				p = pd.read_parquet(pos)
				concated = pd.concat([concated, p])


			CV_table= concated.groupby(['Protein', 'Frames_post_treatment']).Loc_score.agg([lambda x: variation(x), 'count']) #Calculated the varation and the count of cells for each protein in each postion file
			CV_table.rename(columns  = {"lambda_0": "CV_apos", "count" : "cell_count"}) # Should rename, but does not seem to be working
			concated = pd.merge(concated, CV_table, left_on= ['Protein', 'Frames_post_treatment'], right_index=True)

			concated.to_parquet(f"{perc_path}/Chamber_Col_{Date}_r{Run}_Ch{Col}_{today}.parquet")
			return("Complete")
	except: #This is bad practice but should be fine for now
		return(f"Failed for {Col}")

#%%
def var_Prot_combine(col_i, perc_path): #This function runs based on the information file and searched for corresponding data. Done this way becasue it allows for parrallel joining without interferance
	try:
		col_inf = Col_list.iloc[(col_i), :] # col_inf here represents the condtion information search value NOT the name derived
		col_inf_minus = Col_list.iloc[(col_i - 1), :]
		matches = movement_unif_index[(movement_unif_index["Date"] == col_inf["Date"])  & (movement_unif_index["Run"] == col_inf["Run Number"]) & (col_inf_minus["Date"] == col_inf["Date"]) & (movement_unif_index["Col"] > col_inf_minus["MapID (Col_Range)"]) & (movement_unif_index["Col"] <= col_inf["MapID (Col_Range)"])]
		if len(matches) == 0: #Test to see the position being tested corresponds to the first column of the run being tested
			smallmap = Col_list[["Date", "Run Number", "MapID (Col_Range)"]]
			firstof_this_run = smallmap[(smallmap["Run Number"] == col_inf["Run Number"]) & (smallmap["Date"] == col_inf["Date"]) ].iloc[0]# This is a bad way of doing it but was quick to write
			matches = movement_unif_index[(movement_unif_index["Date"] == firstof_this_run["Date"]) & (movement_unif_index["Run"].astype(int) == int(firstof_this_run["Run Number"])) & (movement_unif_index["Col"].astype(int) <= int(firstof_this_run["MapID (Col_Range)"]))] # Take the first row because it is a single ended comparison
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
				p = pd.read_parquet(pos)
				concated = pd.concat([concated, p])


			CV_table= concated.groupby(['Protein', 'Frames_post_treatment']).Loc_score.agg([lambda x: variation(x), 'count']) #Calculated the varation and the count of cells for each protein in each postion file
			CV_table.rename(columns  = {"lambda_0": "CV_apos", "count" : "cell_count"}) # Should rename, but does not seem to be working
			concated = pd.merge(concated, CV_table, left_on= ['Protein', 'Frames_post_treatment'], right_index=True)

			concated.to_parquet(f"{perc_path}/Chamber_Col_{Date}_r{Run}_Ch{Col}_{today}.parquet")
			return("Complete")
	except: #This is bad practice but should be fine for now
		return(f"Failed for {Col}")


#%%

def Chamber_index_er(perc_path):
	os.chdir(perc_path)
	## Create a list of combined column files
	Chamber_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
				if name.endswith("index.parquet"):
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

	def f_date_find(path): #
			if '.parquet' in path:
				name_end = path.find('.parquet')
			elif '.csv' in path:
				name_end = path.find('.csv')
			date_start = name_end - 10
			return(path[date_start:name_end])

	Chamber_index["Chamber"] = pd.Series(Chamber_index.iloc[:,0]).apply(f_chamber)
	Chamber_index['Date_created'] = pd.Series(Chamber_index.iloc[:,0]).apply(f_date_find)
	Chamber_index.sort_values(by = "Date_created", inplace = True, ascending=False)
	latest_date = Chamber_index["Date_created"].iloc[0]
	Chamber_index = Chamber_index.loc[Chamber_index['Date_created'] == latest_date]
	Chamber_index. sort_values(by = "Chamber", inplace = True)
	# Chamber_index = Chamber_index.loc[]

	Chamber_index.to_parquet("Chamber_index.parquet")
	return(Chamber_index)

# %%

#%%
# This is the percentage localization. Below, the percentages are calculated both including and excluding the pre-treatment values

# try: Chamber_index
# except NameError:
# 	os.chdir(post_path)
# 	Chamber_index = pd.read_parquet("Chamber_index.parquet")

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
			try:
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
				percentage_reloc_p.to_parquet(f"percentage_reloc_pt_{chamber}.parquet", index = False)
			except ZeroDivisionError:
				continue

		return(None)
	except: #* Moved the ZeroDivisionError up and added this handler 13/09/2023
		return(f"Failure on {chamber}")

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
			try:
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
			except ZeroDivisionError:
				continue

		percentage_reloc_p.to_parquet(f"percentage_reloc_allt_{chamber}.parquet", index = False)
		return(None)
	except:
		return(f"Failure on {chamber}")


def percentages_bt_manager(chamber, date_found):
	try:
		Chamber_df = pd.read_parquet(f"Chamber_{chamber}_{date_found}.parquet")
	except FileNotFoundError:
		return(f"Col_merge file not found for {chamber}")

	Chamber_df.reset_index(inplace = True, drop = False)
	Chamber_df.sort_values(by = 'Frames_post_treatment', ascending= True, inplace= True) #! This line is absolutely critical.
	# Sort the values by Frames_post_treatment to quantify the amount of relocalization so far
	#This is a new addtion to test whether the there has been relocalization yet

	# def complex_yet(x):
	# 	ind = x["Relocalized"].idxmax()
	# 	does = x.loc[ind]
	# 	x["yet"] = 0
	# 	x.loc[:ind, "yet"] = 0
	# 	x.loc[ind:, "yet"] = does
	# 	return

	def reloc_yet(x):
		ind = x.idxmax() #* Get the first occurrence of max value. This is 1 when relocalised, so the first time in the 'Relocalized' series where true
		does = x.loc[ind] #* Get the value of max value. ie. If max value is 0, then it never relocalizes, whereas if 1 then it does at some point
		x.loc[:ind] = 0 #* Set all values less than the max value as 0
		x.loc[ind:] = does #* Set all values after the max value as 1
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

	post_t_res = new_percentages_post_t(chamber, Chamber_df)
	all_t_res = new_percentages_all_t(chamber, Chamber_df)
	return(f"{chamber} Complete", post_t_res, all_t_res)

#%%
# Calculate the percent trajectory for all_t
def f_percent_trajectory(prot, perc): #This is set to run on the version which includes pre-treatment values. I determined that it is enlightening to know dynamics pre-treatment
	prot_subset = All_pos_allt_t_percentages_melt[(All_pos_allt_t_percentages_melt["Protein"]==prot)]
	prot_subset_pt = prot_subset[prot_subset["Frames_post_treatment"] >= 0]

	max_percent = np.max(prot_subset_pt["Percentage_reloc"])
	t_max = prot_subset[prot_subset["Percentage_reloc"] == max_percent]["Time_post_treatment"].values[0]

	percent_trajectory = (perc/max_percent) * 100
	return(max_percent, percent_trajectory, t_max)

#%%
#Calculate the trajectory for pt_t
def f_percent_trajectory(prot, perc): #This is set to run on the version which includes pre-treatment v alues. I determined that it is enlightening to know dynamics pre-treatment
	prot_subset = All_pos_pt_t_percentages_melt[(All_pos_pt_t_percentages_melt["Protein"]==prot)]
	prot_subset_pt = prot_subset[prot_subset["Frames_post_treatment"] >= 0]

	max_percent = np.max(prot_subset_pt["Percentage_reloc"])
	t_max = prot_subset[prot_subset["Percentage_reloc"] == max_percent]["Time_post_treatment"].values[0]

	percent_trajectory = (perc/max_percent) * 100 if max_percent != 0 else 0
	return(max_percent, percent_trajectory, t_max)


#%%
# line = px.line(dataframe = )
########, Setup for __main__ run
if __name__ == '__main__':
	#, Load in the global_variables as local variables + Change path to post_path to start
	# Global_variables = gv.global_manager()
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

	log_prefix = datetime.datetime.today().strftime('%y_%m_%d-%H_%M_%S')

	#* Convert to local variables within the script
	microfluidic_results = Global_variables['microfluidics_results']
	post_path = Global_variables['post_path']
	list_percetiles = Global_variables['percentiles']
	time_per_frame = Global_variables['timepoint_gap']
	decade = str(datetime.datetime.now().year)[:3] #! I assume that this will not still be state of the art and running New Year's Eve 2029
	pn = Global_variables["cpu_se"] #* This can be changed
	today = str(datetime.date.today())
	year = str(datetime.date.today().year)
	decade = year[:3]

	os.chdir(post_path)

	#, Run the first set of indexers and do subsetting if required. Combine all the small files for position_time into position files
	Quantification_index = Quantification_index_er()
	####### Function run
	if Global_variables['subset'] == True:
		print(f"Analysis will be performed on {Global_variables['subset_by']}s matching {Global_variables['subset_collection']}")
	else:
		pass

	positions = Quantification_index["PositionID"].unique()
	l = len(positions)
	if  l > pn:
		pr = pn
	else:
		pr = l

	######### Function run combining time_pos into pos files
	z = Parallel(n_jobs=pr, verbose = 100)(delayed(combine_pos)(p) for p in positions)
	# ['Good' if v is None else v for v in z]
	z = pd.Series(z)
	PosFrame_combine_file_name = log_prefix + "PosFrame_combine_log.csv"
	z.to_csv(PosFrame_combine_file_name)

	# #* 2023-09-05: Save the log to file
	# if len(z) > 0:
	# 	with open(PosFrame_combine_file_name, 'w+') as w:
	# 		write = csv.writer(w)
	# 		write.writerow(z)

	# PosFrame_combine_file = open(PosFrame_combine_file_name, 'w+')
	# PosFrame_combine_file.write(z) #. This is the line that it failed on with Nikon. I will change this to csv for the mean time.
	# PosFrame_combine_file.close()


	#########

	os.chdir(post_path)
	###### Function run to index
	Quant_ALL_index = Quant_ALL_index_er()
	######

	# #, Get the condition information for associating positions to proteins
	os.chdir(Global_variables['microfluidics_results'])
	Condition_information = pd.read_excel("MicrofluidicsMap_wCol.xlsx", sheet_name='ProteinLocations', dtype={'Date' : str})[["Date", "Run Number", "Time_true", "MapID (Col_Range)", "Myo1Marker", "Protein", "Time"]]
	Condition_information = Condition_information.dropna() #* Get rid of the empty rows, the incomplete logs and the columns which are not used
	Condition_information.drop(Condition_information.columns[Condition_information.columns.str.contains('Unnamed',case = False)],axis = 1, inplace= True)
	Condition_information["Date"] = pd.Series(Condition_information["Date"]).apply(f_convert_date)
	# Condition_information["Time"] = Condition_information["Time.1"] # This is required if there is a "Time" for aquisition, and timing for treatment
	# Condition_information = Condition_information[["Date", "Run Number", "Time", "MapID (Col_Range)", "Myo1Marker", "Protein", "Time"]]

	Col_list = Condition_information[["Date", "MapID (Col_Range)", "Run Number"]].copy()
	Col_list = Col_list.drop_duplicates().reset_index(drop = True)
	Col_list_type_dict = {"Date": str,
					 	"MapID (Col_Range)": int,
						"Run Number": int}
	Col_list = Col_list.astype(Col_list_type_dict)
	# This is the set of combination for the above information components
	# movement_unif_index = pd.read_parquet("index_Movement.parquet")
	#

	#, Begin the next stage for getting Strain_IDs
	os.chdir(post_path)
	l = len(Quant_ALL_index)
	if l < pn:
		pr = l
	else:
		pr = pn # 8 #pn - 1

	os.chdir(post_path)
	###### Function run
	z = Parallel(n_jobs=pr, verbose = 100)(delayed(Strain_ID_multiplex)(p) for p in range(len(Quant_ALL_index)))
	z = pd.Series(z)
	strain_file_name = log_prefix + "Strain_ID_log.csv"
	z.to_csv(strain_file_name)
	#* 2023-09-05: Save the log to file

	Quant_prim_index = Quant_prim_index_er()

	#, Begin the actual movement calculation. This is percetile dependant, so it run for both and put into differnt folders

	for percentile in list_percetiles:
		path_save = os.path.join(post_path, f"{percentile}th_percentile")
		try:
			os.mkdir(path_save)
		except FileExistsError:
			print("Folder already exists")

		##### Run the Movement for percetile and save in individual folder
		z = Parallel(n_jobs=pr, verbose = 100)(delayed(MI_move)(p = p, percentile = percentile,perc_path = path_save) for p in range(len(Quant_prim_index)))
		z = pd.Series(z)
		move_file_name = log_prefix + "MI_move_log.csv"
		z.to_csv(move_file_name)

		#*Index
		movement_unif_index = movement_unif_index_er(perc_path = path_save)
		positions = movement_unif_index["Position"].unique()

		if len(movement_unif_index) > 20: #. Arbitrarily set but reasonable considering most machines will have 8+ cores and fewer than 20 cores. Unless you can make use of 1-2 cycles, it does not make sense to lauch parrallel flow
			Parallel(n_jobs=pr, verbose = 100)(delayed(f_Prot_combine_df)(col_i, perc_path = path_save) for col_i in range(len(Col_list))) #The parallel version should only be used when the number of positions being combined is large
		else:
			Col_len = len(Col_list)
			for col_i in range(Col_len):
				res = f_Prot_combine_df(col_i, perc_path = path_save)
				print(f"{col_i+1}/{Col_len}")
				print(res)
		Chamber_index = Chamber_index_er(perc_path = path_save)

		def f_date_find(path): #* This is just for the combining function.
			if '.parquet' in path:
				name_end = path.find('.parquet')
			elif '.csv' in path:
				name_end = path.find('.csv')
			date_start = name_end - 10
			return(path[date_start:name_end])

		date_found = f_date_find(Chamber_index.iloc[0,0])
		# date_found = Chamber_index.iloc[0,0][-14:-4] #! This is an old legacy version which was likely done for a RUSH version and not removed. Left here until final

		os.chdir(path_save)
		z = Parallel(n_jobs=pr, verbose = 100)(delayed(percentages_bt_manager)(p, date_found) for p in Chamber_index["Chamber"])
		z = pd.Series(z)

		perecentages_bt_log_file_name = log_prefix + "Strain_ID_log.csv"
		z.to_csv(perecentages_bt_log_file_name)

		#########################################
		All_pos_allt_unified_prot_percentages = pd.DataFrame({
			"Frames_post_treatment": []
		})
		########NOTE## The output is All_pos_allt_unified_prot_percentages at ALL TIMEPOINTS

		for c in Chamber_index["Chamber"]:
			try:
				percentage_reloc_p = pd.read_parquet(f"percentage_reloc_allt_{c}.parquet") #. This has no date portion
			except FileNotFoundError:
				print(f"Percentage fiile missing for {c}")
				continue

			All_pos_allt_unified_prot_percentages = pd.merge(All_pos_allt_unified_prot_percentages, percentage_reloc_p, how = 'outer', on='Frames_post_treatment')
			All_pos_allt_unified_prot_percentages["Time_post_treatment"] = All_pos_allt_unified_prot_percentages["Frames_post_treatment"] *time_per_frame

		All_pos_allt_unified_prot_percentages.to_parquet("All_pos_allt_percentage_sync.parquet")
		########NOTE## The output is All_pos_allt_unified_prot_percentages at ALL TIMEPOINTS


		All_pos_pt_unified_prot_percentages = pd.DataFrame({
			"Frames_post_treatment": []
		})

		for c in Chamber_index["Chamber"]:
			try:
				percentage_reloc_p = pd.read_parquet(f"percentage_reloc_pt_{c}.parquet")
			except FileNotFoundError:
				print(f"Percentage fiile missing for {c}")

			All_pos_pt_unified_prot_percentages = pd.merge(All_pos_pt_unified_prot_percentages, percentage_reloc_p, how = 'outer', on='Frames_post_treatment')
			All_pos_pt_unified_prot_percentages["Time_post_treatment"] = All_pos_pt_unified_prot_percentages["Frames_post_treatment"] *time_per_frame
			# All_pos_pt_unified_prot_percentages.drop(columns='Frame', inplace = True)

		# All_pos_pt_unified_prot_percentages = All_pos_pt_unified_prot_percentages.loc[:,~All_pos_pt_unified_prot_percentages.columns.duplicated()]

		All_pos_pt_unified_prot_percentages.to_parquet("All_pos_pt_percentage_sync.parquet")
		#######################

		list_prot_col_t = All_pos_pt_unified_prot_percentages.loc[:, (All_pos_pt_unified_prot_percentages.columns.str.endswith("-perc_t"))].columns

		#Melt percentages for pt
		All_pos_pt_t_percentages_melt = pd.melt(All_pos_pt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_t, var_name='Protein', value_name='Percentage_reloc')
		All_pos_pt_t_percentages_melt["Time_post_treatment"] = All_pos_pt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
		All_pos_pt_t_percentages_melt.to_parquet("All_pos_pt_t_percentages_melt.parquet")

		#Melt percentages for allt
		All_pos_allt_t_percentages_melt = pd.melt(All_pos_allt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_t, var_name='Protein', value_name='Percentage_reloc')
		All_pos_allt_t_percentages_melt["Time_post_treatment"] = All_pos_allt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
		All_pos_allt_t_percentages_melt.to_parquet("All_pos_allt_t_percentages_melt")

		list_prot_col_less = All_pos_pt_unified_prot_percentages.loc[:, (All_pos_pt_unified_prot_percentages.columns.str.endswith("-perc_yet"))].columns

		#Melt percentages for pt_t less
		All_pos_pt_t_less_percentages_melt = pd.melt(All_pos_pt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_less, var_name='Protein', value_name='Percentage_reloc_less')
		All_pos_pt_t_less_percentages_melt["Time_post_treatment"] = All_pos_pt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
		All_pos_pt_t_less_percentages_melt.to_parquet("All_pos_pt_t_less_percentages_melt.parquet")

		#Melt percentages for all_t less
		All_pos_allt_t_less_percentages_melt = pd.melt(All_pos_allt_unified_prot_percentages, id_vars=['Frames_post_treatment'], value_vars=list_prot_col_less, var_name='Protein', value_name='Percentage_reloc_less')
		All_pos_allt_t_less_percentages_melt["Time_post_treatment"] = All_pos_allt_t_percentages_melt["Frames_post_treatment"] *time_per_frame
		All_pos_allt_t_less_percentages_melt.to_parquet("All_pos_allt_t_less_percentages_melt.parquet")

		applied_df= All_pos_allt_t_percentages_melt.apply(lambda row: f_percent_trajectory(row.Protein, row.Percentage_reloc), axis='columns', result_type='expand').rename(columns={0:"Max_percent", 1:"Percent_of_max_percent", 2:"t_max_percent"})
		All_pos_t_percentages_melt= pd.concat([All_pos_allt_t_percentages_melt, applied_df], axis='columns')

		All_pos_t_percentages_melt.to_parquet("All_pos_allt_percentages_melt.parquet", index = False)
		del All_pos_t_percentages_melt
		del applied_df

		applied_df= All_pos_pt_t_percentages_melt.apply(lambda row: f_percent_trajectory(row.Protein, row.Percentage_reloc), axis='columns', result_type='expand').rename(columns={0:"Max_percent", 1:"Percent_of_max_percent", 2:"t_max_percent"})
		All_pos_t_percentages_melt= pd.concat([All_pos_pt_t_percentages_melt, applied_df], axis='columns')

		All_pos_t_percentages_melt.to_parquet("All_pos_pt_percentages_melt.parquet", index = False)

		All_chambers_l = sorted(glob('Chamber_Col_*.parquet'))
		All_chambers_df = pd.concat((pd.read_parquet(file) for file in All_chambers_l), ignore_index= True)

		All_chambers_df.to_parquet("ALL_CHAMBERS.parquet")

		#* Just a simple garbage collector to make sure that no variables have been carried over between percentage loops
		variables_remove = ['All_chambers_df', 'All_chambers_l', 'applied_df', 'movement_unif_index', 'positions', 'Chamber_index', 'date_found', 'All_pos_pt_unified_prot_percentages', 'All_pos_pt_unified_prot_percentages', 'list_prot_col_t', 'All_pos_pt_t_percentages_melt', 'All_pos_allt_t_percentages_melt', 'list_prot_col_less', 'All_pos_pt_t_less_percentages_melt', 'All_pos_allt_t_less_percentages_melt', 'All_pos_allt_unified_prot_percentages']
		for v in variables_remove:
			try:
				exec(f'del {v}')
			except NameError:
				pass
##########################?

#! Commented out for testing and because it produces the incompete informaiton
# 	#, Create special figures
# 	import plotly.express as px
# 	All_pos_t_percentages_melt["Protein"] = All_pos_t_percentages_melt["Protein"].astype('category')

# 	Response_time_comparison = All_pos_t_percentages_melt[["Protein", "t_max_percent"]].drop_duplicates()
# 	Response_time_comparison = Response_time_comparison.sort_values(by = 't_max_percent', ascending=True)

# 	fig = px.bar(data_frame = Response_time_comparison,
# 			y= "Protein",
# 			x= 't_max_percent',
# 			orientation= 'h'
# 			).update_yaxes(categoryorder='total descending')
# 	fig.update_xaxes(nticks=20)
# 	fig.write_html("max_time.html")
# 	fig
# 	# Response_time_comparison
# 	# sns.barplot(data = Response_time_comparison,
# 	# 		x= "Protein",
# 	# 		y= 't_max')

# 	special_proteins = input("Are there any special protein to generate figures for?")


# #%%

# #%%
# def remove(p):
# 	e = p.find("-")
# 	return(p[:e])

# import plotly.express as px
# Het_comparison = All_pos_t_percentages_melt[["Protein", "Max_percent"]].drop_duplicates()
# Het_comparison = Het_comparison.sort_values(by = 'Max_percent')
# Het_comparison["Protein"] = pd.Series(Het_comparison["Protein"]).apply(remove)
# Het_comparison = Het_comparison[Het_comparison["Protein"] != 'FGV2']
# Het_comparison["Sel"] = "reg"
# Het_comparison.set_index("Protein", drop = False)
# Het_comparison.loc[(Het_comparison.Protein == 'FLR1'),'Sel']='spec'
# Het_comparison.loc[(Het_comparison.Protein == 'MRT4'),'Sel']='spec'

# Percentage_ordered_graph = px.bar(data_frame= Het_comparison,
# 		y= "Protein",
# 		x= 'Max_percent',
# 		color = "Sel",
# 		orientation= 'h').update_yaxes(categoryorder='total descending')
# Percentage_ordered_graph.update_xaxes(range=[0, 100])
# Percentage_ordered_graph.update_layout(showlegend=False)
# Percentage_ordered_graph.update_xaxes(nticks=20)
# # sns.barplot(data = Het_comparison,
# # 		x= "Protein",
# # 		y= 'Max_percent')

# Percentage_ordered_graph
# Percentage_ordered_graph.write_html("percent_ordered{pt}.html")
# #%%
# t_less_penetrance = All_pos_pt_t_less_percentages_melt.groupby(by= "Protein").agg('max').reset_index()[["Protein", "Percentage_reloc_less"]]
# t_less_penetrance["Protein"] = pd.Series(t_less_penetrance["Protein"]).apply(remove)
# t_less_penetrance = t_less_penetrance[t_less_penetrance["Protein"] != 'FGV2'] ##### TEMP
# t_less_penetrance["Sel"] = "reg"
# t_less_penetrance.set_index("Protein", drop = False)

# t_less_penetrance.loc[(t_less_penetrance.Protein == 'FLR1'),'Sel']='spec'
# t_less_penetrance.loc[(t_less_penetrance.Protein == 'MRT4'),'Sel']='spec'

# Pen_graph = px.bar(data_frame= t_less_penetrance,
# 		y= "Protein",
# 		x= 'Percentage_reloc_less',
# 		color = "Sel",
# 		orientation= 'h').update_yaxes(categoryorder='total descending')
# Pen_graph.update_xaxes(range=[0, 100])
# Pen_graph.update_layout(showlegend=False)
# Pen_graph.write_html("Penetrance.html")
# Pen_graph

# #%%
# # All_pos_t_percentages_melt
# # t_traj_p = px.line(All_pos_t_percentages_melt, x = 'Frames_post_treatment', y = 'Percent_of_max_percent', color = "Protein")

# # t_traj_p.write_html("t_traj_p.html")
# # #%%

# # All_pos_t_less_percentages_melt
# # # t_less_traj = px.line(All_pos_t_less_percentages_melt, x = 'Frames_post_treatment', y = 'Percentage_reloc_less', color = "Protein")

# # t_less_traj.write_html("t_less_traj.html")

# # #%%

# # #%%

# # t_loc= px.line(All_chambers_df, x = 'Frames_post_treatment', y = 'Loc_score', color = "Protein")