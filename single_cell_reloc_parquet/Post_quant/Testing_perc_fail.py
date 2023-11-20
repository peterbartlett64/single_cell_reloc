#%%
#! This is for testing the perc fail. Last part of post quant
#Todo: If a more stable version is found, then move to post_quant

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
#%%
def f_percentage_index():
	percentage_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith(".parquet") and name.startswith("percentage"):
				percentage_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		# break #This makes the program run non-recursively and not decend into daughter folders

	percentage_index = pd.DataFrame(percentage_index)

	def f_Col_ID(z):
		start = z.find('_Col_')+1 #Note: The shift forward is dependant upon how you write out position
		end = z.find(".parquet")   #Assume that this pipeline will only be used this century! Make sure that the 'n' is present in the array to confirm no place 0 has been lost
		return(z[start:end])

	# percentage_index["Pad"] =
	percentage_index["Col"] = pd.Series(percentage_index.iloc[:,0]).apply(f_Col_ID)
	return(percentage_index)


#%%
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
########, Setup for __main__ run
if __name__ == '__main__':
	#, Load in the global_variables as local variables + Change path to post_path to start
	#. Global_variables = gv.global_manager()
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': 'E:/Microfluidics/MOST FINAL', #. gv.slash_switch(input("Post quant path?")) , #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': int(math.floor(os.cpu_count()*0.7)),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True}

	post_path = Global_variables['post_path']
	list_percentiles = Global_variables['percentiles']
	for percentile in list_percentiles:
		path_save = os.path.join(post_path, f"{percentile}th_percentile")
		try:
			os.mkdir(path_save)
		except FileExistsError:
			print("Folder already exists")

		os.chdir(path_save) #* I don't know how this line was missing????


		# Chamber_index = Chamber_index_er(perc_path = path_save)

		# def f_date_find(path): #* This is just for the combining function.
		# 	if '.parquet' in path:
		# 		name_end = path.find('.parquet')
		# 	elif '.csv' in path:
		# 		name_end = path.find('.csv')
		# 	date_start = name_end - 10
		# 	return(path[date_start:name_end])

		# date_found = f_date_find(path = Chamber_index.iloc[0,0])
		# # date_found = Chamber_index.iloc[0,0][-14:-4] #! This is an old legacy version which was likely done for a RUSH version and not removed. Left here until final

		# os.chdir(path_save)
		# z = Parallel(n_jobs=pr, verbose = 100)(delayed(percentages_bt_manager)(p, date_found) for p in Chamber_index["Chamber"])
		# z = pd.Series(z)

		# perecentages_bt_log_file_name = log_prefix + "Strain_ID_log.csv"
		# z.to_csv(perecentages_bt_log_file_name)

		#########################################
		All_pos_allt_unified_prot_percentages = pd.DataFrame({
			"Frames_post_treatment": []
		})
		########NOTE## The output is All_pos_allt_unified_prot_percentages at ALL TIMEPOINTS

		time_per_frame = Global_variables['timepoint_gap']
		percentage_index = f_percentage_index()

		for c in percentage_index['Col']:
		# for c in Chamber_index["Chamber"]:
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

		for c in percentage_index['Col']:
			try:
				percentage_reloc_p = pd.read_parquet(f"percentage_reloc_pt_{c}.parquet")
			except FileNotFoundError:
				print(f"Percentage fiile missing for {c}")

			All_pos_pt_unified_prot_percentages = pd.merge(All_pos_pt_unified_prot_percentages, percentage_reloc_p, how = 'outer', on='Frames_post_treatment', suffixes=('', '_y'))
			All_pos_pt_unified_prot_percentages.drop(All_pos_pt_unified_prot_percentages.filter(regex='_y$').columns, axis=1, inplace=True)
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