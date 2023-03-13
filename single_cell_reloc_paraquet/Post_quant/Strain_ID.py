import os
import pandas as pd

def Strain_ID_multiplex(p, multiplex = True):
	try:
		Quant_ALL = pd.read_paraquet(Quant_ALL_index.iloc[p,0])
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
			else: #? Not sure which of the options below should be used
				# mKAsignalChange = mKAmax['x99thPercentile_Diff_background_mKate']/mKAmin['factor_mKate_background_Avg'] #* This is not as accurate as below
				mKAsignalChange = (mKAmax['x99thPercentile_Diff_background_mKate'] + mKAmin['factor_mKate_background_Avg'])/mKAmin['factor_mKate_background_Avg']

			mKAsignalChange = abs(mKAsignalChange)#/25)


			# if range_direct == 0: #* This line was not in use
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

		try: Quant_FIN_primary = pd.merge(Quant_ALL, Myo1_info, left_on="Cell_Barcode", right_on = "TrackID_valid")
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
		Quant_FIN_primary.to_paraquet(f"Quant_{Pos}_primary.paraquet")
	except IndexError:
		return(f"IndexError on {p}")
	except ValueError:
		return(f"ValueError on {p}")

def Strain_ID_single(p, recover = False, retreat = False): #! This is not functional yet
	try:
		Quant_ALL = pd.read_paraquet(Quant_ALL_index.iloc[p,0])
		Quant_ALL.drop(Quant_ALL.columns[Quant_ALL.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)

		# frame_max = np.max(Quant_ALL["Frame"])  This is not needed for the Strain_single as the final frame intenstity is not taken into account. In the multiplex version, they Myo1 intensity is used to get stain identity

		try: trackedCells = np.unique(Quant_ALL['Cell_Barcode'])
		except KeyError:
			return(f"Fail on index {p} at trackedCells Stage")

		if len(trackedCells) < 50:       ### Not sure if it would be a good idea to make sure that the position has enough cells. This does not gauruntee enough cells in future factor determination.
			return("There are less than 50 cells in time series")

		Pos = Quant_ALL["Cell_Barcode"].values[0] # This is a temporary and false value, based on the first cell barcode as all within df will have the same position
		Pos = Pos[0:Pos.find("c")]

		s = Pos.find("r")
		Run_n = int(Pos[s+1 : s+2])

		Run_info = Condition_information[(Condition_information["Date"] == Quant_ALL.iloc[0,:]["Date"]) & (Condition_information["Run Number"] == Run_n)]

		# DateCond = Run_info["Date"].values[0]
		Time_treat = Run_info["Time_treat"].values[0]
		if recover == True:
			Time_recover = Run_info["Time_recover"].values[0]
		if retreat == True:
			Time_retreat = Run_info["Time_retreat"].values[0]

		Col = int(Pos[Pos.find("p")+1:-4])  # This is the column of the postion which is being tested
		Pos_info = Run_info[Run_info["MapID (Col_Range)"]>=Col].reset_index (drop=True) # Col HERE represents the NAME and NOT the information. eg. info column 10,20,30 etc.  >= pos_col 10
		Pos_info = Pos_info.iloc[0] # THIS IS A VERY IMPORTANT LINE to make sure that only the
		Prot_strain = Pos_info["Protein"].values[0]
		Col_info = Pos_info["MapID (Col_Range)"].values[0]

		def Section_stage(frame, recovery = recover):
			frame = frame.values[0] # Again, this is for the barcode grouped version

			f_treated = Time_treat/time_per_frame
			if recovery == True and retreat == False:
				f_recovery = Time_recover/time_per_frame
				if frame < f_treated:
					return (0)
				elif frame >= f_treated and frame < f_recovery:
					return(1)
				elif frame >= f_recovery:
					return(2)
			elif recovery == False and retreat == False:
				if frame < f_treated:
					return (0)
				elif frame >= f_treated and < f_recovery:
					return(1)
			elif recover == False and retreat == True: #* This is impossible, so I could either handle with
				return(float(NaN)) #* The function could be written to return an error message, but that would take a lot of downstream handling
			elif recovery == True and retreat == True:
				f_recovery = Time_recover/time_per_frame
				f_retreated = Time_retreat/time_per_frame
				if frame < f_treated:
					print(0)
				elif frame >= f_treated and frame < f_recovery:
					return(1)
				elif frame > f_recovery and frame < f_retreated:
					return(2)
				elif frame >= f_retreated:
					return(3)
			else:
				pass

		Quant_FIN_primary["Col_info"] = int(Col_info)
		# Quant_FIN_primary["Protein"] = pd.Series(Quant_FIN_primary["Myo1Identity"]).apply(Protein_label_multi)
		Quant_FIN_primary["Protein"] = Prot_strain
		# Quant_FIN_primary["Protein"] = pd.Series
		# (Quant_FIN_primary["Myo1Identity"]).apply(Protein_label_multi) #fix the lookup table

		# Quant_FIN_primary["Section_stage"] = pd.Series(Quant_FIN_primary["Frame"]).apply(Section_stage)
		Quant_FIN_primary["Section_stage"] = Quant_FIN_primary.groupby(by = ["ImageID"])["Frame"].transform(Section_stage)

		def f_Frames_post_treatment_shift(i, i_pt):
			v = i - i_pt
			return(v)

		psuedo_map = Quant_FIN_primary[["ImageID", "Is_treated",'In_recovery', "Frame"]].drop_duplicates() # Grab just the imageID, "Is_treated", and "Frame"
		psuedo_map.sort_values(by = ["Frame"], inplace=True)   ####  Add sorting to make sure that the index comparison will be correct

		psuedo_map.reset_index(inplace=True, drop = True) #Reindex to make sure that the index numbers are correct and can be used for comparison
		i_pt = psuedo_map[psuedo_map["Is_treated"] == 1].index[0] - 1 # This isn't the most pretty way to go about but works for now
		psuedo_map["Frames_post_treatment"] = pd.Series(psuedo_map.index).apply(lambda x: f_Frames_post_treatment_shift(x, i_pt))

		psuedo_map.drop(columns = ["Is_treated"], inplace = True)
		Quant_FIN_primary = pd.merge(Quant_FIN_primary, psuedo_map, how = "left", on='ImageID')

		Quant_FIN_primary.to_paraquet(f"Quant_{Pos}_primary.paraquet")
	except IndexError:
		return(f"IndexError on {p}")
	except ValueError:
		return(f"ValueError on {p}")

