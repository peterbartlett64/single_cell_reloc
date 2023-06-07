import pandas as pd
from scipy import stats
import numpy as np
from decimal import DivisionByZero

def I_move(p, percentile):
	try: #TODO: Rename Frame_x to just Frame
		Quant_prim = pd.read_csv(Quant_prim_index.iloc[p,0], usecols = ['Cell_Barcode', 'ImageID', 'Date', 'Frame', 'Unique_Frame', 'factor_median_OBJ_GFP', 'factor_mean_OBJ_GFP', 'x90thPercentile_norm_OBJ_Median_GFP', 'x95thPercentile_norm_OBJ_Median_GFP', 'x99thPercentile_norm_OBJ_Median_GFP', 'Progen_bud', 'x90thPercentile_GFP_RAW', 'x95thPercentile_GFP_RAW', 'x99thPercentile_GFP_RAW', 'max_GFP_RAW', 'max_mKa_RAW', 'max_mKO_RAW', 'TrackID_valid', 'Myo1Identity', 'mKO_foldChange', 'mKO_direction', 'mKA_foldChange', 'mKa_direction', 'byProgen_bud', 'byRange', 'Col_info', 'Protein', 'Is_treated', 'Frames_post_treatment', #* From here on are newly added columns to determine the cell stage. Could do a pd.merge at the end instead.
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
		Quant_f__twe = Quant_prim.loc[(Quant_prim["Frame"] <= 20) & (Quant_prim["Is_treated"] == 0)].copy()
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
