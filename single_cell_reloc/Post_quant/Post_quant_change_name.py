import os
import pandas as pd

global_variables.global_variables as global_variables

def Strain_ID_single(p, recover = False, retreat = False): #! This is not functional yet
	try:
		Quant_ALL = pd.read_csv(Quant_ALL_index.iloc[p,0])
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

		Quant_FIN_primary.to_csv(f"Quant_{Pos}_primary.csv")
	except IndexError:
		return(f"IndexError on {p}")
	except ValueError:
		return(f"ValueError on {p}")

if __name__ == "__main__":
	global_variables()
	Quant_all_index()



