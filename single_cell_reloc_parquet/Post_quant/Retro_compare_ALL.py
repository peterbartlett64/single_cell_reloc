
#, This is to check the effect of taking postitions which are present at the end of the timeseries
#%%
import pandas as pd
import os
import datetime
from single_cell_reloc_parquet.global_functions import global_variables as gv
import math
#%%
#! THIS FILE IS NEWER THAN 'RETRO_COMPARE.PY'

#* For now, only going to read in the one file manually. Should be rather easy to convert up into the index call version
# os.chdir("D:/Sandbox")
# Global_variables['post_path '] = "D:/Sandbox"

year = str(datetime.date.today().year)

#%%
if __name__ == '__main__':
	# Global_Variables = gv.global_manager()
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': 'D:/ALL_FINAL', #. gv.slash_switch(input("Post quant path?")) , #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': int(math.floor(os.cpu_count()*0.7)),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True,
	'figures_root': 'D:/Figures_root'}


	post_combined_path = os.path.join(Global_variables['post_path'], 'Combined_by_perc')
	os.chdir(post_combined_path)

	try:
		info_file = pd.read_parquet('info_merged.parquet')
	except FileNotFoundError:
		os.chdir(Global_variables['microfluidics_results'])
		info_index = pd.read_parquet('info_index.parquet')
		info_index.drop(info_index.columns[info_index.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True) #* This is just in case it somehow made it through

		def drop_test(path):
			if path.__contains__('test'): #* This is a regular string contains
				return(None)
			else:
				return(path)

		info_index['Path'] = pd.Series(info_index['Path']).apply(drop_test)
		info_index.dropna(inplace=True)

		info_list = list(info_index['Path'])

		def f_combine(prefix, suffix):
			return(str(prefix) + "c" + str(suffix).zfill(4))

		c =0
		for file_path in info_list:
			if file_path.endswith('.txt'):
				if c == 0:
					# Load the Parquet file into a DataFrame
					pos = file_path[file_path.find("TrackingResults_d02")+16:file_path.find(".txt")-11]
					info_file = pd.read_csv(file_path, sep = "\t", usecols=["cell_frame", "cell_index", "cell_center_x", "cell_center_y", "cell_majoraxis", "cell_minoraxis", "cell_orientation", "cell_area", "cell_volume", "cell_perimeter", "cell_eccentricity", "cell_fractionOfGoodMembranePixels", "cell_mem_area", "cell_mem_volume", "cell_nuc_radius", "cell_nuc_area", "track_index", "track_start_frame", "track_end_frame"]).reset_index(drop = False)
					info_file = info_file.astype({'track_index': 'str'})
					# info_file['track_index'] = info_file + df['track_index']
					info_file['track_index'] = pd.Series(info_file['track_index']).apply(lambda x: f_combine(prefix = pos, suffix = x))
				else:
					# Load the Parquet file into a DataFrame
					pos = file_path[file_path.find("_d02")+1:file_path.find(".txt")-11]
					df = pd.read_csv(file_path, sep = "\t", usecols=["cell_frame", "cell_index", "cell_center_x", "cell_center_y", "cell_majoraxis", "cell_minoraxis", "cell_orientation", "cell_area", "cell_volume", "cell_perimeter", "cell_eccentricity", "cell_fractionOfGoodMembranePixels", "cell_mem_area", "cell_mem_volume", "cell_nuc_radius", "cell_nuc_area", "track_index", "track_start_frame", "track_end_frame"]).reset_index(drop = False)
					df['track_index'] = pd.Series(df['track_index']).apply(lambda x: f_combine(prefix = pos, suffix = x))
					# Append the data to the info_file DataFrame
					info_file = pd.concat([info_file, df], ignore_index=True)
				print(c)
				c += 1

		info_file['Pos'] = pd.Series(info_file["track_index"]).apply(lambda x: x[: x.find("c")])
		info_file.to_parquet('info_merged.parquet')

	# input('pause here')

	loc_stuff = pd.read_parquet("Final_wAbund.parquet") #*This could be changed to the version which has the localization

	#, Merge files now
	merged = pd.merge(loc_stuff, info_file, how='left', left_on=['Cell_Barcode','Frame'], right_on = ['track_index','cell_frame'])
	merged['track_length'] =  merged['track_end_frame'] - merged['track_start_frame']
	# merged.to_parquet("test_merged.parquet")


	frame_max_ser = merged.groupby('Unique_pos').Frame.agg('max')
	cell_end =pd.DataFrame(merged.groupby('Cell_Barcode').track_end_frame.agg('max')).reset_index(drop = False)
	cell_end['Unique_pos'] = pd.Series(cell_end['Cell_Barcode']).apply(lambda x: x[:x.find("c")])
	check_df = pd.merge(cell_end, pd.DataFrame(frame_max_ser), left_on = 'Unique_pos', right_index=True)
	check_df['diff'] = check_df['Frame'] - check_df['track_end_frame']

	def binary_diff(diff):
		if diff == 0:
			return(1)
		else:
			return(0)

	check_df['pres_end'] = pd.Series(check_df['diff']).apply(binary_diff)
	r_df = check_df[['Cell_Barcode', 'pres_end']].set_index('Cell_Barcode', drop = True)
	del check_df
	merged = merged.merge(r_df, left_on = 'Cell_Barcode', right_index=True)
	merged['c_count_end_nl'] = merged.groupby('Unique_pos').Cell_Barcode.transform(lambda x: len(x.unique()))
	merged.to_parquet("test_merged.parquet")

	merged = merged.loc[merged['track_length'] > 5]
	merged['c_count_end_p5'] = merged.groupby('Unique_pos').Cell_Barcode.transform(lambda x: len(x.unique()))
	merged.to_parquet("merged_data_final.parquet")


	#* Instead of subsetting, it might be a better idea to create a categorical variable about whethere it is present at the end.
	#< same_list = check_df.loc[check_df['diff'] == 0]['Cell_Barcode']
	#< subset = merged.loc[merged['Cell_Barcode'].isin(same_list)]

	# def susbet_func(merged_grp):
	# 	subset = merged_grp.loc[merged_grp['track_end_frame'] == merged_grp['Frame'].unique()[-1]] #* Check that the cells are present in the final frame.
	# 	subset = subset.loc[subset['track_length'] >= 5]

	# subset_df = merged.loc[merged.groupby("Unique_pos").track_end_frame]


	# subset_df = merged.groupby('Unique_pos').apply(susbet_func)
	# subset_df = merged.groupby('Unqiqu_pos').track_end_frame.
	# subset_df.to_parquet('subset_end.parquet', index=True)


#%%
input('somehow worked!')
subset = merged.loc[merged['track_end_frame'] == merged['Frame'].unique()[-1]] #* Check that the cells are present in the final frame.
subset = subset.loc[subset['track_length'] >= 5]
#%%
subset.to_parquet('subset_end.parquet', index=True)
#%%
subset = pPosd.read_parquet("D:\ALL_FINAL\Combined_by_perc\merged_data_final.parquet")
#%%
test = subset.loc[subset['Protein'] == 'EXO70'].reset_index(drop = True)
# print (100 * test.groupby(['Protein','V']).size() / len(df.index))
#%%
new_percs = pd.DataFrame(subset.groupby(['Protein', 'Frames_post_treatment']).Yet.value_counts(normalize=True) *100).rename(columns={'Yet': 'Yet_perc'}).reset_index(drop = False)
new_percs = new_percs.loc[new_percs['Yet'] == 1]

new_percs.to_parquet("D:/ALL_FINAL/Combined_by_perc/new_percs.parquet")

#%%
# This is a modification of Percent Trajectory concept which was calculated in Brandon't pipleine. Percent response is by nature a population metric but this is calculated and stored at the single cell level.
#The following series of functions produces the variabbles
#	1. Gobal percent repsonse for each protein (The number of cells which pass the threshold in total)
#	2. The time specific percentage of cells which are past the threshold
#	3. The percent trajectory of the max percentage response

#< def new_percentages_post_t(chamber, Chamber_df, date_found): #< This function calculated the percentages with just the post-treatment values
#< 	#< Function passed the chamber number and the corresponding df from the master function
#< 	#<Subset the data to just the post treatment
#< 	#<There will be another and seperate percentage calculator that will include pre-treatment
#< 	try:
#< 		#< Cammber_local = Chamber_df.sort_values(by = 'Frames_post_treatment', ascending= True) #< Sort the values again, just in case they have gotten out of order
#< 		Chamber_df_local = Chamber_df.copy()

#< 		#<<NOV23' ADD
#< 		Chamber_df_local_rem = Chamber_df_local[Chamber_df_local["Frames_post_treatment"] < 0]
#< 		#<>

#< 		Chamber_df_local = Chamber_df_local[Chamber_df_local["Frames_post_treatment"] >= 0]

#< 		#< def complex_yet(x):
#< 		#< 	ind = x["Relocalized"].idxmax()
#< 		#< 	does = x.loc[ind]
#< 		#< 	x["yet"] = 0
#< 		#< 	x.loc[:ind, "yet"] = 0
#< 		#< 	x.loc[ind:, "yet"] = does
#< 		#< 	return

#< 		def reloc_yet(x):
#< 			ind = x.idxmax() #<* Get the first occurrence of max value. This is 1 when relocalised, so the first time in the 'Relocalized' series where true
#< 			does = x.loc[ind] #<* Get the value of max value. ie. If max value is 0, then it never relocalizes, whereas if 1 then it does at some point
#< 			x.loc[:ind] = 0 #<* Set all values less than the max value as 0
#< 			x.loc[ind:] = does #<* Set all values after the max value as 1
#< 			return(x)

#< 		def workaround(ind):
#< 			return(Chamber_df_local.loc[ind, "ImageID"])
#< 		def does_workaround(ind):
#< 			return(Chamber_df_local.loc[ind,"Relocalized"])

#< 		Chamber_df_local["Yet"] = Chamber_df_local.groupby(by = "Cell_Barcode")["Relocalized"].transform(reloc_yet) #< This repesents wether there has been relocaliztion yet
#< 		Chamber_df_local["ind"] = Chamber_df_local.groupby(by = "Cell_Barcode")["Relocalized"].transform('idxmax')
#< 		Chamber_df_local["Does"] = pd.Series(Chamber_df_local["ind"]).apply(does_workaround) #<this will work for now but should make it come out of one of the other functions
#< 		Chamber_df_local["When"] = pd.Series(Chamber_df_local["ind"]).apply(workaround)
#< 		Chamber_df_local.drop(columns='ind', inplace = True)

#< 		#<< Added NOV23
#< 		#< Chamber_df_local_rem["Yet"] = 0
#< 		#< Chamber_df_local_rem["Does"] = 0
#< 		#< Chamber_df_local_rem["When"] = 0

#< 		#< Chamber_df_local_rem["Does"] = Chamber_df_local_rem.groupby(by = "Cell_Barcode")["Does"].transform('max')
#< 		#< Chamber_df_local_rem["When"] = Chamber_df_local_rem.groupby(by = "Cell_Barcode")["When"].transform('max')

#< 		#< Chamber_df_local = Chamber_df_local_rem.append(Chamber_df_local)#<. , ignore_index = True)
#< 		#< Chamber_df_local.to_parquet(f"Chamber_{chamber}_{date_found}.parquet") #<. This saves over the older version whcih was read in.
#< 		#<>



#< 		Quant_unif_mKa = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKa']
#< 		Quant_unif_mKO = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKO']


#< 		frames_post_list= Quant_unif_mKa["Frames_post_treatment"].unique()
#< 		#< frame_index = Quant_unif_mKa["ImageID"].unique()
#< 		percentage_reloc_p = pd.DataFrame([])

#< 		try:
#< 			mKa_prot = Quant_unif_mKa["Protein"].values[0]
#< 		except:
#< 			return("No mKa")
#< 		try:
#< 			mKO_prot = Quant_unif_mKO["Protein"].values[0]
#< 		except:
#< 			return("No mKO")
#< 		for f in frames_post_list: #< Loop through the list of post_treatment frame values
#< 			try:
#< 				t_mKa = Quant_unif_mKa[Quant_unif_mKa["Frames_post_treatment"] == f]
#< 				t_mKO = Quant_unif_mKO[Quant_unif_mKO["Frames_post_treatment"] == f]

#< 				reloc_t_mKa = t_mKa[t_mKa["Relocalized"] == 1]
#< 				reloc_t_mKO = t_mKO[t_mKO["Relocalized"] == 1]

#< 				percentage_reloc_t_mKa = (len(reloc_t_mKa)/len(t_mKa))*100
#< 				percentage_reloc_t_mKO = (len(reloc_t_mKO)/len(t_mKO))*100

#< 				t_less_reloc_mKa = t_mKa[t_mKa["Yet"] == 1]
#< 				t_less_reloc_mKO = t_mKO[t_mKO["Yet"] == 1]

#< 				percentage_moved_t_less_mKa = (len(t_less_reloc_mKa)/(len(t_mKa)))*100
#< 				percentage_moved_t_less_mKO = ((len(t_less_reloc_mKO))/(len(t_mKO)))*100


#< 				row = {
#< 					"Frames_post_treatment" : [f],
#< 					f"{mKa_prot}-perc_t" : [percentage_reloc_t_mKa],
#< 					f"{mKO_prot}-perc_t" : [percentage_reloc_t_mKO],
#< 					f"{mKa_prot}-perc_yet": [percentage_moved_t_less_mKa],
#< 					f"{mKO_prot}-perc_yet": [percentage_moved_t_less_mKO]
#< 				}
#< 				f_row = pd.DataFrame(row)

#< 				percentage_reloc_p = pd.concat([percentage_reloc_p, f_row])
#< 				percentage_reloc_p.to_parquet(f"percentage_reloc_pt_{chamber}.parquet", index = False)
#< 			except ZeroDivisionError:
#< 				continue

#< 		return(None)
#< 	except: #<* Moved the ZeroDivisionError up and added this handler 13/09/2023
#< 		return(f"Failure on {chamber}")

#< def new_percentages_all_t(chamber, Chamber_df): #< This function calculated the percentages with all the values. The only way that it differs is in that there is no subsetting for  frames after treatment
#< 	#< Function passed the chamber number and the corresponding df from the master function
#< 	try:
#< 		#< Chamber_df_local = Chamber_df.sort_values(by = 'Frames_post_treatment', ascending= True) #< Sort the values again just in case they have gotten out of order
#< 		Chamber_df_local = Chamber_df.copy()

#< 		#< def complex_yet(x):
#< 		#< 	ind = x["Relocalized"].idxmax()
#< 		#< 	does = x.loc[ind]
#< 		#< 	x["yet"] = 0
#< 		#< 	x.loc[:ind, "yet"] = 0
#< 		#< 	x.loc[ind:, "yet"] = does
#< 		#< 	return

#< 		def reloc_yet(x):
#< 			ind = x.idxmax() #<* Get the first occurrence of max value. This is 1 when relocalised, so the first time in the 'Relocalized' series where true
#< 			does = x.loc[ind] #<* Get the value of max value. ie. If max value is 0, then it never relocalizes, whereas if 1 then it does at some point
#< 			x.loc[:ind] = 0 #<* Set all values less than the max value as 0
#< 			x.loc[ind:] = does #<* Set all values after the max value as 1
#< 			return(x)

#< 		def workaround(ind):
#< 			return(Chamber_df_local.loc[ind, "ImageID"])
#< 		def does_workaround(ind):
#< 			return(Chamber_df_local.loc[ind,"Relocalized"])

#< 		Chamber_df_local["Yet"] = Chamber_df_local.groupby(by = "Cell_Barcode")["Relocalized"].transform(reloc_yet) #< This repesents wether there has been relocaliztion yet
#< 		Chamber_df_local["ind"] = Chamber_df_local.groupby(by = "Cell_Barcode")["Relocalized"].transform('idxmax')
#< 		Chamber_df_local["Does"] = pd.Series(Chamber_df_local["ind"]).apply(does_workaround) #<this will work for now but should make it come out of one of the other functions
#< 		Chamber_df_local["When"] = pd.Series(Chamber_df_local["ind"]).apply(workaround)
#< 		Chamber_df_local.drop(columns='ind', inplace = True)

#< 		Quant_unif_mKa = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKa']
#< 		Quant_unif_mKO = Chamber_df_local[Chamber_df_local["Myo1Identity"] == 'Myo1_mKO']


#< 		frames_post_list= Quant_unif_mKa["Frames_post_treatment"].unique() #<create a list of frames to test
#< 		#< frame_index = Quant_unif_mKa["ImageID"].unique()
#< 		percentage_reloc_p = pd.DataFrame([])

#< 		try:
#< 			mKa_prot = Quant_unif_mKa["Protein"].values[0]
#< 		except:
#< 			return("No mKa")
#< 		try:
#< 			mKO_prot = Quant_unif_mKO["Protein"].values[0]
#< 		except:
#< 			return("No mKO")

#< 		for f in frames_post_list: #< Loop through the list of post_treatment frame values
#< 			try:
#< 				t_mKa = Quant_unif_mKa[Quant_unif_mKa["Frames_post_treatment"] == f]
#< 				t_mKO = Quant_unif_mKO[Quant_unif_mKO["Frames_post_treatment"] == f]

#< 				reloc_t_mKa = t_mKa[t_mKa["Relocalized"] == 1] #<Create a subset of cells which are currently relocalized
#< 				reloc_t_mKO = t_mKO[t_mKO["Relocalized"] == 1]

#< 				percentage_reloc_t_mKa = (len(reloc_t_mKa)/len(t_mKa))*100 #< Divide number of relocalized proteins in the current frame by the total number of cells in the current frame
#< 				percentage_reloc_t_mKO = (len(reloc_t_mKO)/len(t_mKO))*100

#< 				t_less_reloc_mKa = t_mKa[t_mKa["Yet"] == 1]
#< 				t_less_reloc_mKO = t_mKO[t_mKO["Yet"] == 1]

#< 				percentage_moved_t_less_mKa = (len(t_less_reloc_mKa)/(len(t_mKa)))*100
#< 				percentage_moved_t_less_mKO = ((len(t_less_reloc_mKO))/(len(t_mKO)))*100


#< 				row = {
#< 					"Frames_post_treatment" : [f],
#< 					f"{mKa_prot}-perc_t" : [percentage_reloc_t_mKa],
#< 					f"{mKO_prot}-perc_t" : [percentage_reloc_t_mKO],
#< 					f"{mKa_prot}-perc_yet": [percentage_moved_t_less_mKa],
#< 					f"{mKO_prot}-perc_yet": [percentage_moved_t_less_mKO]
#< 				}
#< 				f_row = pd.DataFrame(row)

#< 				percentage_reloc_p = pd.concat([percentage_reloc_p, f_row])
#< 			except ZeroDivisionError:
#< 				continue

#< 		percentage_reloc_p.to_parquet(f"percentage_reloc_allt_{chamber}.parquet", index = False)
#< 		return(None)
#< 	except:
#< 		return(f"Failure on {chamber}")


#< def percentages_bt_manager(chamber, date_found):
#< 	try:
#< 		Chamber_df = subset
#< 		#< Chamber_df = pd.read_parquet(f"Chamber_{chamber}_{date_found}.parquet")
#< 	except FileNotFoundError:
#< 		return(f"Col_merge file not found for {chamber}")


#< 	Chamber_df.reset_index(inplace = True, drop = False)
#< 	Chamber_df.sort_values(by = 'Frames_post_treatment', ascending= True, inplace= True) #<! This line is absolutely critical.
#< 	#< Sort the values by Frames_post_treatment to quantify the amount of relocalization so far
#< 	#<This is a new addtion to test whether the there has been relocalization yet

#< 	#<..... Moved the reloc_yet into the respective functions.

#< 	post_t_res = new_percentages_post_t(chamber, Chamber_df, date_found= date_found) #<. Added datefound for save over because "YET" missing from output
#< 	#< all_t_res = new_percentages_all_t(chamber, Chamber_df)
#< 	all_t_res = "skipped"
#< 	return(f"{chamber} Complete", post_t_res, all_t_res)


#< percentages_bt_manager(chamber= "Col_d0218_r2_Ch60", date_found="2023-10-11")
