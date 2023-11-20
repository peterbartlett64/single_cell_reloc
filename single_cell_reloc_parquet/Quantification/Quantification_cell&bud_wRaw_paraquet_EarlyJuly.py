#%%
# IMPORT ALL PACKAGES REQUIRED
from datetime import datetime
from numpy.lib.function_base import median, quantile
import pandas as pd
import numpy as np
from pandas.core.algorithms import isin, unique
from pandas.io.pytables import to_hdf
from pytest import Instance
import scipy as scipy
import matplotlib.pyplot as plt
#import seaborn as sns
import scipy.io as scipyio
from PIL import Image
import glob
import os
import datetime
#import progressbar
import csv
import statistics
from scipy import stats
from sklearn.cluster import KMeans
import matplotlib
from scipy.io import loadmat
import single_cell_reloc_parquet.global_functions.global_variables as gv

def switch_slash(path):
	new = path.replace(os.sep, '/')
	return (new)

if __name__ == "__main__":
	temp_store = "C:/Users/Peter/Desktop/DATA_TEST"
	# Global_variables = gv.global_manager()
	Global_variables = {'microfluidics_results': gv.slash_switch(input("microfluidics_results?")),
						'subset': True,
						'subset_by':'Unique_frame',
						'subset_collection': ['d0214r1p010200f0023'], #* This is for testing
						# 'subset_by':'Unique_pos',
						# 'subset_collection': ['d0214r1p010200', 'd0214r1p010300', 'd0214r1p030200', 'd0214r1p030300', 'd0214r1p040200', 'd0214r1p040300', 'd0214r1p050200', 'd0214r1p050300', 'd0214r1p070200', 'd0214r1p070300', 'd0214r1p080200', 'd0214r1p080300', 'd0214r1p090200', 'd0214r1p090300', 'd0214r1p110200', 'd0214r1p110300', 'd0214r1p120200', 'd0214r1p120300', 'd0214r1p130200', 'd0214r1p130300', 'd0214r1p150200', 'd0214r1p150300', 'd0214r1p160200', 'd0214r1p160300', 'd0214r2p010200', 'd0214r2p010300', 'd0214r2p030200', 'd0214r2p030300', 'd0214r2p040200', 'd0214r2p040300', 'd0214r2p050200', 'd0214r2p050300', 'd0214r2p070200', 'd0214r2p070300', 'd0214r2p080200', 'd0214r2p080300', 'd0214r2p090200', 'd0214r2p090300', 'd0214r2p110200', 'd0214r2p110300', 'd0214r2p120200', 'd0214r2p120300', 'd0214r2p130200', 'd0214r2p130300', 'd0214r2p150200', 'd0214r2p150300', 'd0214r2p160200', 'd0214r2p160300', 'd0215r1p010200', 'd0215r1p010300', 'd0215r1p030200', 'd0215r1p030300', 'd0215r1p040200', 'd0215r1p040300', 'd0215r1p050200', 'd0215r1p050300', 'd0215r1p070200', 'd0215r1p070300', 'd0215r1p080200', 'd0215r1p080300', 'd0215r1p090200', 'd0215r1p090300', 'd0215r1p110200', 'd0215r1p110300', 'd0215r1p120200', 'd0215r1p120300', 'd0215r1p130200', 'd0215r1p130300', 'd0215r1p150200', 'd0215r1p150300', 'd0215r1p160200', 'd0215r1p160300'],
						'cpu_se': 16
	}
	microfluidics_results = Global_variables['microfluidics_results']

def mask_load(m):
	mask = loadmat(m)
	data = mask['data']
	return (data)

def bud_mask_load(m):
	mask = loadmat(m)
	data = mask['px_data']
	return (data)

def simp_mask_load(m):
	mask = loadmat(m)
	simp = mask['segmentationMask']
	simp[simp >=1 ] = 1
	small = simp
	return (small)

os.chdir(microfluidics_results)


#, Files which are stored under the "microfluidics_results"
#<Indices derived from the analysis folder?
# imgIndex = pd.read_parquet('imgIndex.parquet') #Index of all images to be loaded in
imgIndex = pd.read_csv('imgIndex.csv') #Index of all images to be loaded in
#<Indices derived from segmentaion folder>
# Allmasks = pd.read_parquet('Allmasks.parquet') # Index of masks which show the cell space for internal fluorescence
Allmasks = pd.read_csv('Allmasks.csv') # Index of masks which show the cell space for internal fluorescence

# Allmasks_exp = pd.read_parquet('Allmasks_exp.parquet') #Index of masks which show the expanded cell space for budneck fluorescence
Allmasks_exp = pd.read_csv('Allmasks_exp.csv') #Index of masks which show the expanded cell space for budneck fluorescence
# info_index = pd.read_parquet("info_index.parquet")
info_index = pd.read_csv("info_index.csv", on_bad_lines = 'warn').set_index("Pos")
# Info_all = pd.read_parquet('info_simple.parquet', error_bad_lines=False, warn_bad_lines= True)

# pos_list = imgIndex["Unique_pos"].unique()
# trackDataIndex= pd.read_hdf("MASTER.h5", "trackIndex")
#%% This is to double check that the test versions have been removed
def remove_crap(str_in):
    if "test" in str_in:
        return None
    elif '.sync' in str_in:
        return None
    else:
        return (str_in)

imgIndex['Path'] = pd.Series(imgIndex['Path']).apply(remove_crap)
imgIndex = imgIndex.dropna()
Allmasks['Path'] = pd.Series(Allmasks['Path']).apply(remove_crap)
Allmasks = Allmasks.dropna()
Allmasks_exp['Path'] = pd.Series(Allmasks_exp['Path']).apply(remove_crap)
Allmasks_exp = Allmasks_exp.dropna()


#%%

imgIndex_TP = imgIndex.reset_index(drop = False)
imgIndex_TP.set_index(["Unique_frame", "Channel"], inplace = True, drop = True)
imgIndex_TP.sort_index(inplace=True)
# imgIndex_TP.reset_index(inplace = True)

try:
	segIndex_TP = Allmasks.set_index(["Unique_frame"], drop = False)
except:
	segIndex_TP = Allmasks
segIndex_TP.sort_index(inplace=True)

try:
	bud_segIndex_TP = Allmasks_exp.set_index(["Unique_frame"], drop = False )
except:
	bud_segIndex_TP = Allmasks_exp
bud_segIndex_TP.sort_index(inplace=True)

#. Also read in missing files to make sure not wasting time on incomplete tracks
# Missing_files = pd.read_csv("Files_missing_df.csv", index_col= 0)

today = datetime.date.today()
#%%
# This function is for quantification of a single frame and will output a frame specific csv file.
def Quant_frame(i):
	try:
		frame_n = int(i[i.find("f")+1:])
		unique_frame = i
		unique_pos =  i[:i.find("f")]
		objMeas_DF = pd.DataFrame([])
		date = i[i.find("d"):(i.find("r"))]
		try:
			GFP_path = imgIndex_TP.loc[(i, "GFP"), ["Path"]].values[0]
		except KeyError:
			GFP_path = imgIndex_TP.loc[(i, "gre"), ["Path"]].values[0]
		mKO_path = imgIndex_TP.loc[(i, "mKO"), ["Path"]].values[0]
		mKate_path = imgIndex_TP.loc[(i, "mKa"), ["Path"]].values[0]

		raw_GFP = np.array(Image.open(GFP_path))
		raw_mKO = np.array(Image.open(mKO_path))
		raw_mKate = np.array(Image.open(mKate_path))

		object_pxMAT = mask_load(segIndex_TP.loc[i,["Path"]].values[0])
		bud_pxMAT = bud_mask_load(bud_segIndex_TP.loc[i,["Path"]].values[0])

		GFP_background_Raw = raw_GFP[object_pxMAT == 0]
		factor_GFP_background_Med = np.median(GFP_background_Raw)
		factor_GFP_background_Avg = np.mean(GFP_background_Raw)
		factor_GFP_background_Tot = np.sum(GFP_background_Raw)

		mKO_background_Raw = raw_mKO[bud_pxMAT == 0] # laod in the mKO area based on the expanded mask which should better capture the budneck
		backgroundKO_val = np.percentile(mKO_background_Raw, 75) # determine background KO at 75th percentile
		factor_mKO_background_Med = np.median(mKO_background_Raw)
		factor_mKO_background_Avg = np.mean(mKO_background_Raw)
		factor_mKO_background_Tot = np.sum(mKO_background_Raw)

		mKate_background_Raw = raw_mKate[bud_pxMAT == 0] # laod in the mKa non-cell space based on the expanded mask which should better capture the budneck
		backgroundKA_val = np.percentile(mKate_background_Raw, 75) #determine the backgound mKa at 75th percentile
		factor_mKate_background_Med = np.median(mKate_background_Raw)
		factor_mKate_background_Avg = np.mean(mKate_background_Raw)
		factor_mKate_background_Tot = np.sum(mKate_background_Raw)

		# object_pxMAT_copy = simp_mask_load(Allmask.iloc[k,0])  ### This is a legeacy line in case using the sparse function
		## directly in the mask load function

		def pos_info_read(unique_pos):
			path =  info_index.loc[unique_pos, "Path"]
			#! The below must be left as read_csv because it is taking in tabular data generated by tracking
			Info_all = pd.read_csv(path, sep = "\t", usecols = ['cell_frame', 'cell_index', 'cell_majoraxis', 'cell_minoraxis', 'cell_area', 'cell_volume', 'cell_perimeter', 'cell_eccentricity', 'cell_fractionOfGoodMembranePixels', 'cell_mem_area', 'cell_mem_volume', 'cell_nuc_radius', 'cell_nuc_area', 'cell_pole1_age', 'cell_pole2_age', 'cell_timepoint', 'track_index', 'track_fingerprint_real_distance', 'track_age', 'track_parent', 'track_parent_frame', 'track_parent_prob', 'track_parent_score', 'track_generation', 'track_lineage_tree', 'track_cell_cycle_phase', 'track_assignment_fraction', 'track_start_frame', 'track_end_frame', 'track_index_corrected', 'track_has_bud', 'track_budneck_total'])

			Info_all['Unique_pos'] = unique_pos

			Info_all["Unique_frame"] = Info_all["Unique_pos"] + "f" + Info_all["cell_frame"].astype(str).str.zfill(4)
			Info_all["Circ(major/minor)"] = Info_all["cell_majoraxis"]/Info_all["cell_minoraxis"]
			Info_all["Unique_cell"] = Info_all['Unique_pos'] + "c" + Info_all["track_index"].astype(str).str.zfill(4) ########### SUPER IMPORTANT. Correction to track_index from cell_index
			Info_all.set_index(['Unique_cell'], inplace = True)
			Info_all["Track_length"] = Info_all["track_end_frame"] - Info_all["track_start_frame"]

			Info_all["Circ(major/minor)"] = Info_all["cell_majoraxis"]/Info_all["cell_minoraxis"]
			Info_all["Track_length"] = Info_all["track_end_frame"] - Info_all["track_start_frame"]
			Info_all = Info_all.loc[:,["track_index", "Unique_frame", "track_start_frame", "Track_length", "cell_area", "cell_volume", "cell_perimeter", "Circ(major/minor)"]] # Further subset the file
			# too_short_drop = info_simple[info_simple["Track_length"] < 6].reset_index()["Unique_cell"] ## This was subsetting based on the
			# curdate = str(datetime.date.today())
			# too_short_drop.to_parquet(f"Dropped_bc_short_{curdate}")
			# info_simple = info_simple[info_simple["Track_length"] >= 6] #Filter info to those which are present for more than 6 frames. For now, this has been set as the arbritrary lower limit of frames required to determine current localization.
			Info_all.reset_index(inplace=True)
			return(Info_all)

		Info_all = pos_info_read(unique_pos)
		Cell_index_TP = Info_all[Info_all["Unique_frame"] == i].reset_index(drop =True)

		for k in Cell_index_TP["track_index"]:
			try:
				# row = Cell_index_TP.iloc[k ,:] #! This is WRONG
				row = Cell_index_TP.loc[Cell_index_TP['track_index'] == k].iloc[0] #? This returns as a series as before
				c_barcode = row["Unique_cell"]
				tc = row["track_index"]

				bud_actv = Info_all[Info_all["Unique_cell"] == c_barcode]["track_start_frame"].values[0]
				# if bud_actv > 2:
				# 	bud_actv = bud_actv - 2

				if frame_n == bud_actv:
					progen_bud = 1
				else:
					progen_bud = 0

				# start_frame = row["track_start_frame"]
				# if start_frame == f"{unique_frame[:-4]}{frame_n-2}":#Substract 2 frames and test if this is the start frame for the cell
				# 	start_min2_myo_mKa = raw_mKate
				# 	start_min2_myo_mKO = raw_mKate
				# GFP_aura_background_Raw = raw_GFP[bud_pxMAT == tc]

				# load in the cell specific fluorescence area
				GFP_object = raw_GFP[object_pxMAT == tc] #load in the GFP values based on the regular cell space
				mKO_object = raw_mKO[bud_pxMAT == tc] #load in the mKO values based on the expanded cell space
				mKate_object = raw_mKate[bud_pxMAT == tc] #load in the mKa values based on the expanded cell space

				#####
				factor_median_OBJ_GFP = np.median(GFP_object)
				factor_mean_OBJ_GFP = np.mean(GFP_object)
				factor_total_OBJ_GFP = np.sum(GFP_object)

				factor_median_OBJ_KO = np.median(mKO_object)
				factor_mean_OBJ_KO = np.mean(mKO_object)
				factor_total_OBJ_KO = np.sum(mKO_object)

				factor_median_OBJ_mKate = np.median(mKate_object)
				factor_mean_OBJ_mKate = np.mean(mKate_object)
				factor_total_OBJ_mKate = np.sum(mKate_object)

				#######################################
				averageIntensity_GFP_Frame = np.mean(raw_GFP) #average intesnity of frame including the CELLS and NON-CELL SPACES
				averageIntesntiy_GFP_Background = np.mean(GFP_background_Raw) #average intensity of NON-CELL space
				averageIntensity_GFP_Object= np.mean(GFP_object) #average intensity of GFP within cell

				## Normalizing by the median, mean and sum of intensity of the cell object
				# normGFPbyMedian_object_less_aura = (GFP_object - GFP_aura_background_Raw)/factor_median_OBJ_GFP

				# normGFPbyMedian_object_less_background = (GFP_object - GFP_background_Raw)/factor_median_OBJ_GFP


				normGFPbyMedian_object = GFP_object/factor_median_OBJ_GFP
				normGFPbyMean_object = GFP_object/factor_mean_OBJ_GFP
				normGFPbyTotal_intensity_object = GFP_object/factor_total_OBJ_GFP

				## Normalizing by the median mean and sum of intensity of the non-cell space
				normGFPbyBackground_median = GFP_object/factor_GFP_background_Med
				normGFPbyBackground_mean = GFP_object/factor_GFP_background_Avg
				normGFPbyBackground_total = GFP_object/factor_GFP_background_Tot


				## Decile bins for normalized against OBJECT (cell) MEDIAN
				# x10thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 10, axis = 0)
				# x20thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 20, axis = 0)
				# x30thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 30, axis = 0)
				# x40thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 40, axis = 0)
				# x50thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 50, axis = 0)
				# x60thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 60, axis = 0)
				# x70thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 70, axis = 0)
				x80thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 80, axis = 0)
				x90thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 90, axis = 0)
				x95thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 95, axis = 0)
				x99thPercentile_norm_OBJ_Median_GFP = np.percentile(normGFPbyMedian_object, 99, axis = 0) # This is for finding foci

				## Decile bins for normalized against OBJECT (cell) MEAN
				# x10thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 10, axis = 0)
				# x20thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 20, axis = 0)
				# x30thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 30, axis = 0)
				# x40thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 40, axis = 0)
				# x50thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 50, axis = 0)
				# x60thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 60, axis = 0)
				# x70thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 70, axis = 0)
				# x80thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 80, axis = 0)
				x90thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 90, axis = 0)
				x95thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 95, axis = 0)
				x99thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(normGFPbyMean_object, 99, axis = 0)  # This is for finding foci

				## Decile bins for normalized against OBJECT (cell) TOTAL
				# x20thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 20, axis = 0)
				# x10thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 10, axis = 0)
				# x30thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 30, axis = 0)
				# x40thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 40, axis = 0)
				# x50thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 50, axis = 0)
				# x60thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 60, axis = 0)
				# x70thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 70, axis = 0)
				# x80thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 80, axis = 0)
				x90thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 90, axis = 0)
				x95thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 95, axis = 0)
				x99thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(normGFPbyTotal_intensity_object, 99, axis = 0)  # This is for finding foci

				## Decile bins normalized against BACKGROUND intesntity MEDIAN
				# x10thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 10, axis = 0)
				# x20thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 20, axis = 0)
				# x30thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 30, axis = 0)
				# x40thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 40, axis = 0)
				# x50thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 50, axis = 0)
				# x60thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 60, axis = 0)
				# x70thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 70, axis = 0)
				# x80thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 80, axis = 0)
				x90thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 90, axis = 0)
				x95thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 95, axis = 0)
				x99thPercentile_norm_BKGRND_Median_GFP = np.percentile(normGFPbyBackground_median, 99, axis = 0)  # This is for finding foci

				##  Decile bins normalized against BACKGROUND intenstiy MEAN
				# x10thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 10, axis = 0)
				# x20thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 20, axis = 0)
				# x30thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 30, axis = 0)
				# x40thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 40, axis = 0)
				# x50thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 50, axis = 0)
				# x60thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 60, axis = 0)
				# x70thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 70, axis = 0)
				# x80thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 80, axis = 0)
				x90thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 90, axis = 0)
				x95thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 95, axis = 0)
				x99thPercentile_norm_BKGRND_Mean_GFP = np.percentile(normGFPbyBackground_mean, 99, axis = 0)  # This is for finding foci

				##  Decile bins normalized against BACKGROUND intensity TOTAL
				# x10thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 10, axis = 0)
				# x20thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 20, axis = 0)
				# x30thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 30, axis = 0)
				# x40thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 40, axis = 0)
				# x50thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 50, axis = 0)
				# x60thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 60, axis = 0)
				# x70thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 70, axis = 0)
				# x80thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 80, axis = 0)
				x90thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 90, axis = 0)
				x95thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 95, axis = 0)
				x99thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(normGFPbyBackground_total, 99, axis = 0)  # This is for finding foci

				GFP_spread = len(GFP_object)

				##### Metrics from KO intensity
				averageIntensity_mKO_Frame = np.mean(raw_mKO) #average intesnity of frame including the CELLS and NON-CEll SPACES
				averageIntesntiy_mKO_Background = np.mean(mKO_background_Raw) #average intensity of NON-CELL space
				averageIntensity_mKO_Object= np.mean(mKO_object) #average intensity of GFP within cell

				segCellmKOvalues = mKO_object - backgroundKO_val # subtract background
				finalKOvalues = segCellmKOvalues.copy()#*segCellmKOvalues*segCellmKOvalues # improve signal-noise ratio

				normKObyMedian_object = mKO_object/factor_median_OBJ_mKate
				normKObyMean_object = mKO_object/factor_mean_OBJ_mKate
				normKObyTotal_intensity_object = mKO_object/factor_total_OBJ_mKate

				# normKObyBackground_median = mKO_object/factor_mKO_background_Med
				# normKObyBackground_mean = mKO_object/factor_mKO_background_Avg
				# normKObyBackground_total = mKO_object/factor_mKO_background_Tot

				x80thPercentile_Diff_background_mKO = np.percentile(finalKOvalues, 80)
				x90thPercentile_Diff_background_mKO = np.percentile(finalKOvalues, 90)
				x99thPercentile_Diff_background_mKO = np.percentile(finalKOvalues, 99)


				x60thPercentile_GFP_RAW = np.percentile(GFP_object,60)
				x80thPercentile_GFP_RAW = np.percentile(GFP_object,80)
				x90thPercentile_GFP_RAW = np.percentile(GFP_object,90)
				x95thPercentile_GFP_RAW = np.percentile(GFP_object,95)
				x99thPercentile_GFP_RAW = np.percentile(GFP_object,99)
				max_GFP_RAW = max(GFP_object)

				x60thPercentile_mKa_RAW = np.percentile(mKate_object,60)
				x80thPercentile_mKa_RAW = np.percentile(mKate_object,80)
				x90thPercentile_mKa_RAW = np.percentile(mKate_object,90)
				x95thPercentile_mKa_RAW = np.percentile(mKate_object,95)
				x99thPercentile_mKa_RAW = np.percentile(mKate_object,99)
				max_mKa_RAW = max(mKate_object)

				x60thPercentile_mKO_RAW = np.percentile(mKO_object,60)
				x80thPercentile_mKO_RAW = np.percentile(mKO_object,80)
				x90thPercentile_mKO_RAW = np.percentile(mKO_object,90)
				x95thPercentile_mKO_RAW = np.percentile(mKO_object,95)
				x99thPercentile_mKO_RAW = np.percentile(mKO_object,99)
				max_mKO_RAW = max(mKO_object)

				# x90thPercentile_norm_OBJ_Median_mKO = np.percentile(normKObyMedian_object, 90)
				# x99thPercentile_norm_OBJ_Median_mKO = np.percentile(normKObyMedian_object, 99)
				# x90thPercentile_norm_OBJ_Mean_mKO = np.percentile(normKObyMean_object, 90)
				# x99thPercentile_norm_OBJ_Mean_mKO = np.percentile(normKObyMean_object, 99)
				# x90thPercentile_norm_OBJ_Total_intensity_mKO = np.percentile(normKObyTotal_intensity_object, 90)
				# x99thPercentile_norm_OBJ_Total_intensity_mKO = np.percentile(normKObyTotal_intensity_object, 99)

				# x90thPercentile_norm_BKGRND_Median_mKO = np.percentile(normKObyBackground_median, 90)
				# x99thPercentile_norm_BKGRND_Median_mKO = np.percentile(normKObyBackground_median, 99)
				# x90thPercentile_norm_BKGRND_Mean_mKO = np.percentile(normKObyBackground_mean, 90)
				# x99thPercentile_norm_BKGRND_Mean_mKO = np.percentile(normKObyBackground_mean, 99)
				# x90thPercentile_norm_BKGRND_Total_intensity_mKO = np.percentile(normKObyBackground_total, 90)
				# x99thPercentile_norm_BKGRND_Total_intensity_mKO = np.percentile(normKObyBackground_total, 99)


				mKO_spread = len(mKO_object)

				#### Metrics from mKate intensity
				averageIntensity_mKate_Frame = np.mean(raw_mKate) #average intesnity of frame including the CELLS and NON-CEll SPACES
				averageIntesntiy_mKate_Background = np.mean(mKate_background_Raw) #average intensity of NON-CELL space
				averageIntensity_mKate_Object= np.mean(mKate_object) #average intensity of GFP within cell

				# normKatebyMedian_object = mKO_object/factor_median_OBJ_mKate
				# normKatebyMean_object = mKO_object/factor_mean_OBJ_mKate
				# normKatebyTotal_intensity_object = mKO_object/factor_total_OBJ_mKate

				segCellmKAvalues = mKate_object - backgroundKA_val
				finalKAvalues = segCellmKAvalues.copy()#*segCellmKAvalues*segCellmKAvalues

				normKatebyBackground_median = mKate_object/factor_mKate_background_Med
				normKatebyBackground_mean = mKate_object/factor_mKate_background_Avg
				normKatebyBackground_total = np.sum(mKate_object)/factor_mKate_background_Tot

				#Percentiles of Cubic
				x80thPercentile_Diff_background_mKate = np.percentile(finalKAvalues, 80)
				x90thPercentile_Diff_background_mKate = np.percentile(finalKAvalues, 90)
				x99thPercentile_Diff_background_mKate = np.percentile(finalKAvalues, 99)

				# x90thPercentile_norm_OBJ_Median_mKate = np.percentile(normKatebyMedian_object, 90)
				# x99thPercentile_norm_OBJ_Median_mKate = np.percentile(normKatebyMedian_object, 99)
				# x90thPercentile_norm_OBJ_Mean_mKate = np.percentile(normKatebyMean_object, 90)
				# x99thPercentile_norm_OBJ_Mean_mKate = np.percentile(normKatebyMean_object, 99)
				# x90thPercentile_norm_OBJ_Total_intensity_mKate = np.percentile(normKatebyTotal_intensity_object, 90)
				# x99thPercentile_norm_OBJ_Total_intensity_mKate = np.percentile(normKatebyTotal_intensity_object, 99)

				# x90thPercentile_norm_BKGRND_Median_mKate = np.percentile(normKatebyBackground_median, 90)
				# x99thPercentile_norm_BKGRND_Median_mKate = np.percentile(normKatebyBackground_median, 99)
				# x90thPercentile_norm_BKGRND_Mean_mKate = np.percentile(normKatebyBackground_mean, 90)
				# x99thPercentile_norm_BKGRND_Mean_mKate = np.percentile(normKatebyBackground_mean, 99)
				# x90thPercentile_norm_BKGRND_Total_intesnity_mKate = np.percentile(normKatebyBackground_total, 90)
				# x99thPercentile_norm_BKGRND_Total_intensity_mKate = np.percentile(normKatebyBackground_total, 99)

				mKate_spread = len(mKate_object)

				# Create the metric row
				object_measurements = {
									'Cell_Barcode' : c_barcode,
									'ImageID' : i,
									'Date' : date,
									'Frame': frame_n,
									'Unique_Frame': unique_frame, # This was added on Novemebr 30th because it is required for merging with info
									# 'x10thPercentile_norm_OBJ_Median_GFP' : [x10thPercentile_norm_OBJ_Median_GFP],
									# 'x20thPercentile_norm_OBJ_Median_GFP' : [x20thPercentile_norm_OBJ_Median_GFP],
									# 'x30thPercentile_norm_OBJ_Median_GFP' : [x30thPercentile_norm_OBJ_Median_GFP],
									# 'x40thPercentile_norm_OBJ_Median_GFP' : [x40thPercentile_norm_OBJ_Median_GFP],
									# 'x50thPercentile_norm_OBJ_Median_GFP' : [x50thPercentile_norm_OBJ_Median_GFP],
									# 'x60thPercentile_norm_OBJ_Median_GFP' : [x60thPercentile_norm_OBJ_Median_GFP],
									# 'x70thPercentile_norm_OBJ_Median_GFP' : [x70thPercentile_norm_OBJ_Median_GFP],
									'x80thPercentile_norm_OBJ_Median_GFP' : [x80thPercentile_norm_OBJ_Median_GFP],
									'x90thPercentile_norm_OBJ_Median_GFP' : [x90thPercentile_norm_OBJ_Median_GFP],
									'x95thPercentile_norm_OBJ_Median_GFP' : [x95thPercentile_norm_OBJ_Median_GFP],
									'x99thPercentile_norm_OBJ_Median_GFP' : [x99thPercentile_norm_OBJ_Median_GFP],

									# 'x10thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x10thPercentile_norm_OBJ_Mean_Intensity_GFP],
									# 'x20thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x20thPercentile_norm_OBJ_Mean_Intensity_GFP],
									# 'x30thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x30thPercentile_norm_OBJ_Mean_Intensity_GFP],
									# 'x40thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x40thPercentile_norm_OBJ_Mean_Intensity_GFP],
									# 'x50thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x50thPercentile_norm_OBJ_Mean_Intensity_GFP],
									# 'x60thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x60thPercentile_norm_OBJ_Mean_Intensity_GFP],
									# 'x70thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x70thPercentile_norm_OBJ_Mean_Intensity_GFP],
									# 'x80thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x80thPercentile_norm_OBJ_Mean_Intensity_GFP],
									'x90thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x90thPercentile_norm_OBJ_Mean_Intensity_GFP],
									'x95thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x95thPercentile_norm_OBJ_Mean_Intensity_GFP],
									'x99thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x99thPercentile_norm_OBJ_Mean_Intensity_GFP],

									# 'x10thPercentile_norm_OBJ_Total_Intensity_GFP' : [x10thPercentile_norm_OBJ_Total_Intensity_GFP],
									# 'x20thPercentile_norm_OBJ_Total_Intensity_GFP' : [x20thPercentile_norm_OBJ_Total_Intensity_GFP],
									# 'x30thPercentile_norm_OBJ_Total_Intensity_GFP' : [x30thPercentile_norm_OBJ_Total_Intensity_GFP],
									# 'x40thPercentile_norm_OBJ_Total_Intensity_GFP' : [x40thPercentile_norm_OBJ_Total_Intensity_GFP],
									# 'x50thPercentile_norm_OBJ_Total_Intensity_GFP' : [x50thPercentile_norm_OBJ_Total_Intensity_GFP],
									# 'x60thPercentile_norm_OBJ_Total_Intensity_GFP' : [x60thPercentile_norm_OBJ_Total_Intensity_GFP],
									# 'x70thPercentile_norm_OBJ_Total_Intensity_GFP' : [x70thPercentile_norm_OBJ_Total_Intensity_GFP],
									# 'x80thPercentile_norm_OBJ_Total_Intensity_GFP' : [x80thPercentile_norm_OBJ_Total_Intensity_GFP],
									'x90thPercentile_norm_OBJ_Total_Intensity_GFP' : [x90thPercentile_norm_OBJ_Total_Intensity_GFP],
									'x95thPercentile_norm_OBJ_Total_Intensity_GFP' : [x95thPercentile_norm_OBJ_Total_Intensity_GFP],
									'x99thPercentile_norm_OBJ_Total_Intensity_GFP' : [x99thPercentile_norm_OBJ_Total_Intensity_GFP],

									# 'x10thPercentile_norm_BKGRND_Median_GFP' : [x10thPercentile_norm_BKGRND_Median_GFP],
									# 'x20thPercentile_norm_BKGRND_Median_GFP' : [x20thPercentile_norm_BKGRND_Median_GFP],
									# 'x30thPercentile_norm_BKGRND_Median_GFP' : [x30thPercentile_norm_BKGRND_Median_GFP],
									# 'x40thPercentile_norm_BKGRND_Median_GFP' : [x40thPercentile_norm_BKGRND_Median_GFP],
									# 'x50thPercentile_norm_BKGRND_Median_GFP' : [x50thPercentile_norm_BKGRND_Median_GFP],
									# 'x60thPercentile_norm_BKGRND_Median_GFP' : [x60thPercentile_norm_BKGRND_Median_GFP],
									# 'x70thPercentile_norm_BKGRND_Median_GFP' : [x70thPercentile_norm_BKGRND_Median_GFP],
									# 'x80thPercentile_norm_BKGRND_Median_GFP' : [x80thPercentile_norm_BKGRND_Median_GFP],
									'x90thPercentile_norm_BKGRND_Median_GFP' : [x90thPercentile_norm_BKGRND_Median_GFP],
									'x95thPercentile_norm_BKGRND_Median_GFP' : [x95thPercentile_norm_BKGRND_Median_GFP],
									'x99thPercentile_norm_BKGRND_Median_GFP' : [x99thPercentile_norm_BKGRND_Median_GFP],

									# 'x10thPercentile_norm_BKGRND_Mean_GFP' : [x10thPercentile_norm_BKGRND_Mean_GFP],
									# 'x20thPercentile_norm_BKGRND_Mean_GFP' : [x20thPercentile_norm_BKGRND_Mean_GFP],
									# 'x30thPercentile_norm_BKGRND_Mean_GFP' : [x30thPercentile_norm_BKGRND_Mean_GFP],
									# 'x40thPercentile_norm_BKGRND_Mean_GFP' : [x40thPercentile_norm_BKGRND_Mean_GFP],
									# 'x50thPercentile_norm_BKGRND_Mean_GFP' : [x50thPercentile_norm_BKGRND_Mean_GFP],
									# 'x60thPercentile_norm_BKGRND_Mean_GFP' : [x60thPercentile_norm_BKGRND_Mean_GFP],
									# 'x70thPercentile_norm_BKGRND_Mean_GFP' : [x70thPercentile_norm_BKGRND_Mean_GFP],
									# 'x80thPercentile_norm_BKGRND_Mean_GFP' : [x80thPercentile_norm_BKGRND_Mean_GFP],
									'x90thPercentile_norm_BKGRND_Mean_GFP' : [x90thPercentile_norm_BKGRND_Mean_GFP],
									'x95thPercentile_norm_BKGRND_Mean_GFP' : [x95thPercentile_norm_BKGRND_Mean_GFP],
									'x99thPercentile_norm_BKGRND_Mean_GFP' : [x99thPercentile_norm_BKGRND_Mean_GFP],

									# 'x10thPercentile_norm_BKGRND_Total_intensity_GFP' : [x10thPercentile_norm_BKGRND_Total_intensity_GFP],
									# 'x20thPercentile_norm_BKGRND_Total_intensity_GFP' : [x20thPercentile_norm_BKGRND_Total_intensity_GFP],
									# 'x30thPercentile_norm_BKGRND_Total_intensity_GFP' : [x30thPercentile_norm_BKGRND_Total_intensity_GFP],
									# 'x40thPercentile_norm_BKGRND_Total_intensity_GFP' : [x40thPercentile_norm_BKGRND_Total_intensity_GFP],
									# 'x50thPercentile_norm_BKGRND_Total_intensity_GFP' : [x50thPercentile_norm_BKGRND_Total_intensity_GFP],
									# 'x60thPercentile_norm_BKGRND_Total_intensity_GFP' : [x60thPercentile_norm_BKGRND_Total_intensity_GFP],
									# 'x70thPercentile_norm_BKGRND_Total_intensity_GFP' : [x70thPercentile_norm_BKGRND_Total_intensity_GFP],
									# 'x80thPercentile_norm_BKGRND_Total_intensity_GFP' : [x80thPercentile_norm_BKGRND_Total_intensity_GFP],
									'x90thPercentile_norm_BKGRND_Total_intensity_GFP' : [x90thPercentile_norm_BKGRND_Total_intensity_GFP],
									'x95thPercentile_norm_BKGRND_Total_intensity_GFP' : [x95thPercentile_norm_BKGRND_Total_intensity_GFP],
									'x99thPercentile_norm_BKGRND_Total_intensity_GFP' : [x99thPercentile_norm_BKGRND_Total_intensity_GFP],

									'averageIntensity_GFP_Frame' : [averageIntensity_GFP_Frame],
									'averageIntesntiy_GFP_Background' : [averageIntesntiy_GFP_Background],
									'averageIntensity_GFP_Object' : [averageIntensity_GFP_Object],
									'GFP_spread': [GFP_spread],
									# 'GFP_background': [GFP_background_Raw],

									'Progen_bud': [progen_bud],

									'x80thPercentile_Diff_background_mKate' : [x80thPercentile_Diff_background_mKate],
									'x90thPercentile_Diff_background_mKate' : [x90thPercentile_Diff_background_mKate],
									'x99thPercentile_Diff_background_mKate' : [x99thPercentile_Diff_background_mKate],

									# 'x90thPercentile_norm_OBJ_Median_mKO' : [x90thPercentile_norm_OBJ_Median_mKO],
									# 'x99thPercentile_norm_OBJ_Median_mKO' : [x99thPercentile_norm_OBJ_Median_mKO],
									# 'x90thPercentile_norm_OBJ_Mean_mKO' : [x90thPercentile_norm_OBJ_Mean_mKO],
									# 'x99thPercentile_norm_OBJ_Mean_mKO' : [x99thPercentile_norm_OBJ_Mean_mKO],
									# 'x90thPercentile_norm_OBJ_Total_intensity_mKO' : [x90thPercentile_norm_OBJ_Total_intensity_mKO],
									# 'x99thPercentile_norm_OBJ_Total_intensity_mKO' : [x99thPercentile_norm_OBJ_Total_intensity_mKO],
									# 'x90thPercentile_norm_BKGRND_Median_mKO' : [x90thPercentile_norm_BKGRND_Median_mKO],
									# 'x99thPercentile_norm_BKGRND_Median_mKO' : [x99thPercentile_norm_BKGRND_Median_mKO],
									# 'x90thPercentile_norm_BKGRND_Mean_mKO' : [x90thPercentile_norm_BKGRND_Mean_mKO],
									# 'x99thPercentile_norm_BKGRND_Mean_mKO' : [x99thPercentile_norm_BKGRND_Mean_mKO],
									# 'x90thPercentile_norm_BKGRND_Total_intensity_mKO' : [x90thPercentile_norm_BKGRND_Total_intensity_mKO],
									# 'x99thPercentile_norm_BKGRND_Total_intensity_mKO' : [x99thPercentile_norm_BKGRND_Total_intensity_mKO],

									'x80thPercentile_Diff_background_mKO' : [x80thPercentile_Diff_background_mKO],
									'x90thPercentile_Diff_background_mKO' : [x90thPercentile_Diff_background_mKO],
									'x99thPercentile_Diff_background_mKO' : [x99thPercentile_Diff_background_mKO],

									'averageIntensity_mKO_Frame' : [averageIntensity_mKO_Frame],
									'averageIntesntiy_mKO_Background' : [averageIntesntiy_mKO_Background],
									'averageIntensity_mKO_Object' : [averageIntensity_mKO_Object],
									'mKO_spread': [mKO_spread],
									# 'mKO_background': [mKO_background_Raw],

									# 'x90thPercentile_norm_OBJ_Median_mKate' : [x90thPercentile_norm_OBJ_Median_mKate],
									# 'x99thPercentile_norm_OBJ_Median_mKate' : [x99thPercentile_norm_OBJ_Median_mKate],
									# 'x90thPercentile_norm_OBJ_Mean_mKate' : [x90thPercentile_norm_OBJ_Mean_mKate],
									# 'x99thPercentile_norm_OBJ_Mean_mKate' : [x99thPercentile_norm_OBJ_Mean_mKate],
									# 'x90thPercentile_norm_OBJ_Total_intensity_mKate' : [x90thPercentile_norm_OBJ_Total_intensity_mKate],
									# 'x99thPercentile_norm_OBJ_Total_intensity_mKate' : [x99thPercentile_norm_OBJ_Total_intensity_mKate],
									# 'x90thPercentile_norm_BKGRND_Median_mKate' : [x90thPercentile_norm_BKGRND_Median_mKate],
									# 'x99thPercentile_norm_BKGRND_Median_mKate' : [x99thPercentile_norm_BKGRND_Median_mKate],
									# 'x90thPercentile_norm_BKGRND_Mean_mKate' : [x90thPercentile_norm_BKGRND_Mean_mKate],
									# 'x99thPercentile_norm_BKGRND_Mean_mKate' : [x99thPercentile_norm_BKGRND_Mean_mKate],
									# 'x90thPercentile_norm_BKGRND_Total_Intesnity_mKate' : [x90thPercentile_norm_BKGRND_Total_intesnity_mKate],
									# 'x99thPercentile_norm_BKGRND_Total_intensity_mKate' : [x99thPercentile_norm_BKGRND_Total_intensity_mKate],

									# 'x80thPercentile_Diff_background_mKa' : [x80thPercentile_Diff_background_mKate],
									# 'x90thPercentile_Diff_background_mKa' : [x90thPercentile_Diff_background_mKate],
									# 'x99thPercentile_Diff_background_mKa' : [x99thPercentile_Diff_background_mKate],

									'averageIntensity_mKate_Frame' : [averageIntensity_mKate_Frame],
									'averageIntesntiy_mKate_Background' : [averageIntesntiy_mKate_Background],
									'averageIntensity_mKate_Object' : [averageIntensity_mKate_Object],
									'mKate_spread': [mKate_spread],
									# 'mKate_background': [mKate_background_Raw],

									# Show calcultion factors specific to the cell and frame number (timepoint)
									'factor_median_OBJ_GFP' : [factor_median_OBJ_GFP],
									'factor_mean_OBJ_GFP' : [factor_mean_OBJ_GFP],
									'factor_total_OBJ_GFP' : [factor_total_OBJ_GFP],
									'factor_GFP_background_Med' : [factor_GFP_background_Med],
									'factor_GFP_background_Avg' : [factor_GFP_background_Avg],
									'factor_GFP_background_Tot' : [factor_GFP_background_Tot],

									'factor_median _OBJ_KO' : [factor_median_OBJ_KO],
									'factor_mean_OBJ_KO' : [factor_mean_OBJ_KO],
									'factor_total_OBJ_KO' : [factor_total_OBJ_KO],
									'factor_mKO_background_Med' : [factor_mKO_background_Med],
									'factor_mKO_background_Avg' : [factor_mKO_background_Avg],
									'factor_mKO_background_Tot' : [factor_mKO_background_Tot],

									'factor_median_OBJ_mKate' : [factor_median_OBJ_mKate],
									'factor_mean_OBJ_mKate' : [factor_mean_OBJ_mKate],
									'factor_total_OBJ_mKate' : [factor_total_OBJ_mKate],
									'factor_mKate_background_Med' : [factor_mKate_background_Med],
									'factor_mKate_background_Avg' : [factor_mKate_background_Avg],
									'factor_mKate_background_Tot' : [factor_mKate_background_Tot],

									'x80thPercentile_GFP_RAW':[x80thPercentile_GFP_RAW],
									'x60thPercentile_GFP_RAW':[x60thPercentile_GFP_RAW],
									'x90thPercentile_GFP_RAW':[x90thPercentile_GFP_RAW],
									'x95thPercentile_GFP_RAW':[x95thPercentile_GFP_RAW],
									'x99thPercentile_GFP_RAW':[x99thPercentile_GFP_RAW],
									'max_GFP_RAW' : [max_GFP_RAW],
									'x60thPercentile_mKa_RAW' : [x60thPercentile_mKa_RAW],
									'x80thPercentile_mKa_RAW' : [x80thPercentile_mKa_RAW],
									'x90thPercentile_mKa_RAW' : [x90thPercentile_mKa_RAW],
									'x95thPercentile_mKa_RAW' : [x95thPercentile_mKa_RAW],
									'x99thPercentile_mKa_RAW' : [x99thPercentile_mKa_RAW],
									'max_mKa_RAW' : [max_mKa_RAW],
									'x60thPercentile_mKO_RAW': [x60thPercentile_mKO_RAW],
									'x80thPercentile_mKO_RAW': [x80thPercentile_mKO_RAW],
									'x90thPercentile_mKO_RAW': [x90thPercentile_mKO_RAW],
									'x95thPercentile_mKO_RAW': [x95thPercentile_mKO_RAW],
									'x99thPercentile_mKO_RAW': [x99thPercentile_mKO_RAW],
									'max_mKO_RAW' : [max_mKO_RAW]
									}

				object_measurements = pd.DataFrame(object_measurements)
				objMeas_DF = pd.concat([objMeas_DF, object_measurements])
			except:
				print(f"Fail on cell {k} for {i}")
				pass
		# os.chdir(microfluidics_results)
		# # objMeas_DF.to_hdf('MASTER.h5', key = i, mode='r+')
		objMeas_DF.to_parquet(f"Quantification_{i}_{today}.parquet")
		# objMeas_DF.to_csv(f"Quantification_{i}_{today}.csv")
		# objMeas_DF.to_parquet(f"Quantification_{i}_{today}.parquet")


		done = f"{i} completed"
		return (done)
	except:
		return(f"Failure on {i}")

from joblib import Parallel, delayed

# p = "d0221r1p500200"


# imgIndex_TPt = imgIndex_TP[imgIndex_TP["Unique_pos"] == p]['Unique_frame'].reset_index(drop = True)

# Run the quantification only on the frames which have been properly tracked. Below an index of frames which have been properly tracked is created
# Allmasks = Allmasks[Allmasks["Unique_pos"] != "d0218r2p180200"]
# Allmasks = Allmasks.iloc[763:]
#%%

try:
	path = os.path.join(microfluidics_results, str(today))
	os.mkdir(path)
	print("Creating folder")
except FileExistsError:
	path = os.path.join(microfluidics_results, str(today))
	print("Folder already exists")
	pass

# path = "D:/Microfluidics/Missing_RESULTS2/2022-06-15"

### THis is to dertermine the files that have already been completed. Could make recursive to test for all dates, but at this stage do not want to because there have been several tests within the "RESULTS" root
# os.chdir("D:/Microfluidics/RESULTS/2022-03-09")5360

os.chdir('E:/Microfluidics/RESULTS/2023-07-19')
Quantification_index = []

count = 0
for root, dirs, files, in os.walk(os.getcwd()):
	for name in files:
		if name.endswith(".parquet") and name.startswith("Quantification_d"): # fix naming
			Quantification_index.append({'Path': os.path.join(root, name)})
			count = count + 1
			print(count, end="\r")
		elif name.endswith(".csv") and name.startswith("Quantification_d"): # fix naming
			Quantification_index.append({'Path': os.path.join(root, name)})
			count = count + 1
			print(count, end="\r")
		else:
			pass
	# break #.This makes the program run non-recursively and not decend into daughter folders
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

try:
	Quantification_index["PositionID"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Position_ID)
	Quantification_index["Date"] = pd.Series(Quantification_index["PositionID"]).apply(f_expdate)
	Quantification_index["Frame"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Frame)
	# Quantification_index["Mix" ] =
	# Quantification_index = Quantifica tion_index[Quantification_index["PositionID"] != "ind"] # This works for now
	Quantification_index.to_parquet("Quantification_index.parquet")
	Completed_quant_index = Quantification_index["Frame"]
except:
	print("Failure feature one of extraction functionss")


#< 	# Completed_quant_index = pd.Series([]) #* This was to accept if there are not completed series.
#< 	#.This is a MAJOR issue if one of the above functions jsut fails

Allmasks_exp_sub = Allmasks_exp
# Allmasks_exp_sub = Allmasks_exp_sub.loc[~(Allmasks_exp_sub["Unique_frame"].isin(Completed_quant_index))] #. This is temp commented out


os.chdir(path)
# del Quantification_index
# del Completed_quant_index
# del count

# Quantification_index = []
# count = 0
# for root, dirs, files, in os.walk(os.getcwd()):
# 	for name in files:
# 		if name.endswith(".parquet") and name.startswith("Quantification_d"): # fix naming
# 			Quantification_index.append({'Path': os.path.join(root, name)})
# 			count = count + 1
# 			print(count, end="\r")
# 		elif name.endswith(".csv") and name.startswith("Quantification_d"): # fix naming
# 			Quantification_index.append({'Path': os.path.join(root, name)})
# 			count = count + 1
# 			print(count, end="\r")
# 		else:
# 			pass
# 	# break #.This makes the program run non-recursively and not decend into daughter folders

# if len(Quantification_index) == 0: #* This is a the newer version which should be a lot safer.
# 	Completed_quant_index = pd.Series([])

# Quantification_index = pd.DataFrame(Quantification_index)

# year = str(datetime.datetime.now().year)
# decade = year[:3] # This assumes that the analysis is done within the same decade as starting
# del year

# def f_Position_ID(z):
# 	start = z.find('tion_')+5 #Note: The shift forward is dependant upon how you write out position
# 	end = z.find("_"+ decade)-5   #Assume that this pipeline will only be used this century! Make sure that the 'n' is present in the array to confirm no place 0 has been lost
# 	return(z[start:end])

# def f_Frame(z):
# 	start = z.find('tion_')+5
# 	end = z.find('_' + "20")
# 	return(z[start:end])

# def f_expdate(x):
# 	dend = x.find('r')
# 	expdate = x[0:dend]
# 	return(expdate)

# try:
# 	Quantification_index["PositionID"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Position_ID)
# 	Quantification_index["Date"] = pd.Series(Quantification_index["PositionID"]).apply(f_expdate)
# 	Quantification_index["Frame"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Frame)
# 	# Quantification_index["Mix" ] =
# 	# Quantification_index = Quantifica tion_index[Quantification_index["PositionID"] != "ind"] # This works for now
# 	Quantification_index.to_parquet("Quantification_index.parquet")
# 	Completed_quant_index = Quantification_index["Frame"]
# except:
# 	print("Failure feature one of extraction functionss")

# # Allmasks_exp_sub.reset_index(inplace=True, drop = False)

# if len(Completed_quant_index) > 0:
# 	Allmasks_exp_sub = Allmasks_exp_sub.loc[~(Allmasks_exp_sub["Unique_frame"].isin(Completed_quant_index))]


#* Check that also not included in missing_files

#< Allmasks_exp_sub = Allmasks_exp_sub.loc[~(Allmasks_exp_sub["Unique_pos"].isin(Missing_files["Unique_pos"]))]
#< Allmasks_exp_sub = Allmasks_exp_sub.loc[~(Allmasks_exp_sub["Date"].isin(skip_date))]


if Global_variables['subset'] == True:
	Allmasks_exp_sub = Allmasks_exp_sub.loc[Allmasks_exp_sub[Global_variables['subset_by']].isin(Global_variables['subset_collection'])]
else:
	pass

# Allmasks_exp_sub = Allmasks_exp_sub.loc[Allmasks_exp_sub["Date"] == 'd0212'].sort_values(by ="Unique_frame") #TEMP

# Allmasks = Allmasks[(Allmasks["Unique_pos"] == 'd0222r1p1200200')]
MaskIndexTPt = Allmasks_exp_sub["Unique_frame"]
MaskIndexTPt =  MaskIndexTPt.unique()
# MaskIndexTPt = MaskIndexTPt[~MaskIndexTPt

# del Quantification_index # Clean up the RAM. This could get large for a big microfluidics experiement.
#%%
pn = os.cpu_count()
pr = int(input(f"Machine has {pn} cores. How many do you wish to run Quant with?"))

print(f"function_starting with {pr} cores")
Parallel(n_jobs=pr, verbose = 100)(delayed(Quant_frame)(i) for i in MaskIndexTPt)


# #Pad specific mapping.
# Quantification_index  = []
# count = 0
# for root, dirs, files, in os.walk(os.getcwd()):
# 	for name in files:
# 		print(name)
# 		if name.endswith(".parquet") and name.startswith("Quantification"): # fix naming
# 			Quantification_index.append({'Path': os.path.join(root, name)})
# 			count = count + 1
# 			print(count, end="\r")
# 		else:
# 			pass
# Quantification_index = pd.DataFrame(Quantification_index)


# def f_Position_ID(z):
# 	start = z.find('tion_')+5 #Note: The shift forward is dependant upon how you write out position
# 	end = z.find("_"+decade)-5   #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
# 	return(z[start:end])

# def f_Frame(z):
# 	start = z.find('_time')+5
# 	end = z.find('_time')+9
# 	return(z[start:end])

# def f_expdate(x):
# 	dend = x.find('r')
# 	expdate = x[0:dend]
# 	return(expdate)


# # Quantification_index["Pad"] =
# Quantification_index["PositionID"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Position_ID)
# Quantification_index["Date"] = pd.Series(Quantification_index["PositionID"]).apply(f_expdate)
# # Quantification_index["Mix" ] =
# Quantification_index = Quantification_index[Quantification_index["PositionID"] != "ind"] # This works for now
# Quantification_index.to_parquet("Quantification_index.parquet")


# from joblib import Parallel, delayed

# def combine_pos(pos):
# 	Quant_frame_comb = pd.DataFrame([])
# 	subset = Quantification_index[Quantification_index["PositionID"] == pos]
# 	for qf in range(len(subset)):
# 		q = pd.read_parquet(subset.iloc[qf,0])
# 		Quant_frame_comb = pd.concat([Quant_frame_comb, q])
# 	Quant_frame_comb.to_parquet(f"Quant_{pos}_ALL.parquet")
# 	return(f"{pos} merging is complete")

# pn = os.cpu_count()
# pr = pn - 1

# positions = Quantification_index["PositionID"].unique()
# Parallel(n_jobs=pr)(delayed(combine_pos)(p) for p in positions)


# %%
