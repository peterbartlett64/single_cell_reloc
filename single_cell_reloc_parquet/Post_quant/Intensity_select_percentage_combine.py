#%%
import pandas as pd
import math
from scipy import stats
from joblib import Parallel, delayed
from datetime import date
import os
import psutil as p
import math
from datetime import datetime
from joblib import Parallel, delayed
from glob import glob
from single_cell_reloc_parquet.Post_quant.abundance_gen import Abundance_log_manager
# from tqdm import tqdm

#. This version was derived from the quick versions found in 'current' which is no longer current
#! Quite a bit of stuff was removed from this version. The most recent version before this was the quick version

#%% Create a merged Chamber index table

def var_Chamber_index(post_quant_root, percentile):
	post_qunat_folder_percentile = f"{percentile}th_percentile" #* There is no need to account for dates, as this will be within a dated post_quant folder
	full_path = os.path.join(os.path.join(post_quant_root, post_qunat_folder_percentile))

	year = str(datetime.today().year)
	def f_chamber(col):
		s = col.find("_Col")+1
		e = col.find("_" + year)
		Chamber = col[s:e]
		return(Chamber)

	def run_perc(path): #* This will give the percentage that that it was run for. Folder names can be variable but shoudl have 'th' immediately following the percentage number
		start = path.find("th") - 2
		end = path.find("th")
		#* This could be changed to just path[0:2]
		return(path[start:end])

	if os.path.exists(full_path) == True:
		pass
	else:
		t = 0
		while t == 0:
			post_qunat_folder_percentile = input(f'No folder for {percentile}th percentile found at {full_path}. Please input folder name relative to {post_quant_root}')
			full_path = os.path.join(os.path.join(post_quant_root, post_qunat_folder_percentile))
			if os.path.exists == True:
				t = 1
			else:
				pass
	Chamber_index_temp = []
	count = 0
	for root, dirs, files in os.walk(full_path):
		for name in files:
			if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
				if name.endswith("index.parquet"):
					pass
				else:
					Path_n = f"Path_{percentile}"
					Chamber_index_temp.append({Path_n: os.path.join(root, name)})
					count = count + 1
					print(count, end="\r")
			else:
				pass
		break #This makes the program run non recursively and not decend into daughter folders
	Chamber_index_temp = pd.DataFrame(Chamber_index_temp)
	Chamber_index_temp["Chamber"] = pd.Series(Chamber_index_temp.iloc[:,0]).apply(f_chamber)
	Chamber_index_temp["Percent_run"] = pd.Series(Chamber_index_temp.iloc[:,0]).apply(run_perc)
	Chamber_index_temp = pd.DataFrame(Chamber_index_temp)
	return(Chamber_index_temp)

#%% This was moved up and then retired. Appended index vs. merged
# #, Non-recursively import all the chamber values and join into 95th and 99th percentile tables
# Chamber_index = []
# count = 0
# for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_99))):
# 	for name in files:
# 		if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
# 			if name.endswith("index.parquet"):
# 				pass
# 			else:
# 				Chamber_index.append({'Path': os.path.join(root, name)})
# 				count = count + 1
# 				print(count, end="\r")
# 		else:
# 			pass
# 	break #This makes the program run non recursively and not decend into daughter folders
# for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_95))):
# 	for name in files:
# 		if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
# 			if name.endswith("index.parquet"):
# 				pass
# 			else:
# 				Chamber_index.append({'Path': os.path.join(root, name)})
# 				count = count + 1
# 				print(count, end="\r")
# 		else:
# 			pass
# 	break #This makes the program run non recursively and not decend into daughter folders

# Chamber_index = pd.DataFrame(Chamber_index)


# def f_chamber(col):
# 	s = col.find("_Col")+1
# 	e = col.find("_" + year)
# 	Chamber = col[s:e]
# 	return(Chamber)

# def run_perc(path): #* This will give the percentage that that it was run for. Folder names can be variable but shoudl have 'th' immediately following the percentage number
# 	start = path.find("th") - 2
# 	end = path.find("th")
# 	return(path[start:end])

# Chamber_index["Chamber"] = pd.Series(Chamber_index.iloc[:,0]).apply(f_chamber)
# Chamber_index["Percent_run"] = pd.Series(Chamber_index.iloc[:,0]).apply(run_perc) #*Extract the run perctage from the path. The folder will contain "xxth"
# Chamber_index.to_parquet(f"Chamber_index_combined.parquet")


#%%
#, NEW version of Abundance genration.

# def Abundance_calc_manager(local_prot_df:pd.DataFrame): #* The abundacne could be incorporated in the post-quant function, but it could slow down. It also makes the process less modular and therefore makes timing of computer location a restraint
# 	#, Calculate the abundance for every cell at every timepoint. This is done with such density so that associations can be made to the current abundance and the Loc_score
# 	local_prot_df["Abundance"] = pd.Series(local_prot_df["factor_median_OBJ_GFP"]).apply(lambda x: math.log(x)) #* Log transform absorbance values in a rough metric of molecules per cell. Not true measure because cannot related to baselenine dataset
# 	#* purposely overwrit the old version. Because the dataframes are quite large, it is excessive to keep the unmodified df in memory for a couple non-mutative lines
# 	#? Perfroming with the lambda function might be slightly faster because the function does not need to be declared for each core/run

# 	group_prot_frame = local_prot_df.groupby(["Protein", "Frame"]) #. This may need to be changed to jsut regular "Frame"

# 	# #* This seems like a good plac to grab the z_scores for the Loc_score and Abundance
# 	local_prot_df["z_score_Loc"] = group_prot_frame["Loc_score"].transform(stats.zscore)
# 	local_prot_df["z_score_Abund"] = group_prot_frame["Abundance"].transform(stats.zscore)
# 	return(local_prot_df)

#%%
#, This should be called after var_Chamber_index()
def Tag_percentage_run(r: int, percentile:int): #* This version should run faster than previous and be more stable
	row = Chamber_index.iloc[r] #* interate into the colum list

	#. Instead of consulting the information table again, which would have to open a file, open the the lower percentile and if that is correct continue on
	percentile_x = row[f"Path_{percentile}"]
	if os.path.exists(percentile_x) == True:
		pass
	else: #* A bit of error handling because the shift from pd.DataFrame to series or element can be unpredictable
		percentile_x = percentile_x[0]
		if os.path.exists == False:
			return(f"Could not find a path for {r}")
		else:
			pass

	percentile_x_file = pd.read_parquet(percentile_x)
	col_proteins = percentile_x_file.Protein.unique()

	fail = []
	for prot in col_proteins:
		try:
			selected = percentage_to_use_lib.loc[prot, "Selected_series"]
			if selected == percentile:
				prot_df = percentile_x_file.loc[percentile_x_file['Protein'] == prot].copy()
			else:
				prot_df_path = row[f"Path_{selected}"]
				if os.path.exists(prot_df_path) == True:
					pass
				else: #* A bit of error handling because the shift from pd.DataFrame to series or element can be unpredictable
					prot_df_path = prot_df_path[0]
					if os.path.exists == False:
						return(f"Could not find a path for {r}")
					else:
						pass
				prot_df = pd.read_parquet(prot_df_path)
				prot_df = prot_df.loc[prot_df['Protein']== prot]

			#Todo: change the below to proper name instead of the aliase given in the if '__main__'
			prot_df = Abundance_calc_manager(local_prot_df=prot_df)
			prot_df.to_parquet(f'{prot}_selected.parquet') #* This is just to keep track that the abundance has not been attached yet. Can be easily renamed with batch_rename
		except KeyError:
			fail.append(f"KeyError on {prot}")
			continue
		except OSError:
			fail.append(f"KeyError on {prot}")
	return(fail)
	# 	path = row["Path"]
	# 	percent = int(row["Percent_run"])

	# 	temp_df = pd.read_parquet(path)

	# 	proteins = temp_df["Protein"].unique()

	# 	try:
	# 		protein_one = proteins[0]
	# 		if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_one) & (percentage_to_use_lib["Selected_series"] == percent)]) >0:
	# 			prot_one_df = temp_df.loc[temp_df["Protein"] == protein_one]

	# 			prot_one_df.to_parquet(f"{protein_one}.parquet")
	# 		else:
	# 			pass
	# 	except Exception as e:
	# 		return(e, f"{path}")

	# 	try:
	# 		protein_two = proteins[1]
	# 		if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_two) & (percentage_to_use_lib["Selected_series"] == percent)])>0:
	# 			prot_two_df.to_parquet(f"{protein_two}.parquet")

	# 			prot_two_df = temp_df.loc[temp_df["Protein"] == protein_two]
	# 		else:
	# 			pass
	# 	except Exception as e:
	# 		return(e, f"{path}")

	# 	return(f"{path} complete")
	# except Exception as e:
	# 	return(e, f"{path}")

# def Tag_percentage_run(r):
# 	try:
# 		row = Chamber_index.loc[r] #* For list of chamber values
# 		path = row["Path"]
# 		percent = int(row["Percent_run"])

# 		temp_df = pd.read_parquet(path)

# 		proteins = temp_df["Protein"].unique()

# 		try:
# 			protein_one = proteins[0]
# 			if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_one) & (percentage_to_use_lib["Selected_series"] == percent)]) >0:
# 				prot_one_df = temp_df.loc[temp_df["Protein"] == protein_one]

# 				prot_one_df.to_parquet(f"{protein_one}.parquet")
# 			else:
# 				pass
# 		except Exception as e:
# 			return(e, f"{path}")

# 		try:
# 			protein_two = proteins[1]
# 			if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_two) & (percentage_to_use_lib["Selected_series"] == percent)])>0:
# 				prot_two_df.to_parquet(f"{protein_two}.parquet")

# 				prot_two_df = temp_df.loc[temp_df["Protein"] == protein_two]
# 			else:
# 				pass
# 		except Exception as e:
# 			return(e, f"{path}")

# 		return(f"{path} complete")
# 	except Exception as e:
# 		return(e, f"{path}")
#%%
# Chamber_index_95 = []
# count = 0
# for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_95))):
# 	for name in files:
# 		if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
# 			if name.endswith("index.parquet"):
# 				pass
# 			else:
# 				Chamber_index_95.append({'Path_95': os.path.join(root, name)})
# 				count = count + 1
# 				print(count, end="\r")
# 		else:
# 			pass
# 	break #This makes the program run non recursively and not decend into daughter folders
# Chamber_index_95 = pd.DataFrame(Chamber_index_95)
# Chamber_index_95["Chamber"] = pd.Series(Chamber_index_95.iloc[:,0]).apply(f_chamber)
# Chamber_index_95["Percent_run"] = pd.Series(Chamber_index_95.iloc[:,0]).apply(run_perc)

# Chamber_index_99 = []
# count = 0
# for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_99))):
# 	for name in files:
# 		if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
# 			if name.endswith("index.parquet"):
# 				pass
# 			else:
# 				Chamber_index_99.append({'Path_99': os.path.join(root, name)})
# 				count = count + 1
# 				print(count, end="\r")
# 		else:
# 			pass
# 	break #This makes the program run non recursively and not decend into daughter folders
# Chamber_index_99 = pd.DataFrame(Chamber_index_99)
# Chamber_index_99["Chamber"] = pd.Series(Chamber_index_99.iloc[:,0]).apply(f_chamber)
# Chamber_index_99["Percent_run"] = pd.Series(Chamber_index_99.iloc[:,0]).apply(run_perc)

# Chamber_index = pd.merge(Chamber_index_95,Chamber_index_99, left_on = 'Chamber', right_on = 'Chamber', how = 'outside')
# Chamber_index.to_parquet(f"Chamber_index_combined.parquet")
# #############################################>
#%% Main Run stuff
if __name__ == '__main__':
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': "D:/ALL_FINAL", #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': os.cpu_count(),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True}

	today = str(date.today())
	Abundance_calc_manager = Abundance_log_manager #.Set an alternative name for the imported function. This is temporary but just to speed up before the version in this script is removed

	post_quant_root = Global_variables['post_path']
	os.chdir(post_quant_root)

	percentage_to_use = pd.read_parquet("Final_combined_comparison.parquet") #* This is fine to be hard coded as it will be produced just prior with fixed name from the percentile selection script
	percentage_to_use_lib = percentage_to_use.loc[:, ["Protein", "Selected_series"]] #* Get a list of which series to be used

	def conv_str_number(str_number):#* this is needed right now but will be fixed in next version when switch to xxth_percentile
		# if str_number == 'ninetynine':
		# 	return(99)
		# if str_number == 'ninetyfive':
		# 	return(95)
		#. Corrected below
		return(int(str_number[0:2]))

	percentage_to_use_lib["Selected_series"] = pd.Series(percentage_to_use_lib["Selected_series"]).apply(conv_str_number)
	percentage_to_use_lib.set_index('Protein', inplace= True)
	percentage_to_use_lib.to_parquet('percentage_to_use_lib.parquet')

	# #* Here if the path is given explicitly, run
	# post_quant_sub_99 = slash_switch(input("Where are the 99th percentile files? Folder name"))
	# post_quant_sub_95 = slash_switch(input("Where are the 95th percentile files? Folder name"))
	# year = str(date.today().year)

	p_n = 0
	for p in Global_variables['percentiles']:
		if p_n == 0: #* This first call is just to deal with the first row, so that Chamber_index does not have to be defined outside
			Chamber_index= var_Chamber_index(post_quant_root = post_quant_root, percentile = p)
		else:
			Chamber_index_temp = var_Chamber_index(post_quant_root = post_quant_root, percentile = p)
			Chamber_index = pd.merge(Chamber_index, Chamber_index_temp, left_on = 'Chamber', right_on = 'Chamber', how = 'outer')
		p_n += 1

	Chamber_index = Chamber_index.dropna()
	Chamber_index.set_index('Chamber', inplace= True)
	# Chamber_index_na = Chamber_index.loc[~(Chamber_index.dropna())]
	# if len(Chamber_index_na) > 0:
	# 	print(Chamber_index_na)
	# state = input('Continue?')

	# if state == 'y' or state == 'yes':
	# 	pass
	# else:
	# 	print('Chamber index aysmetry and input cuased termination')
	# 	exit

	#, Create folder to copy the chosen files into
	folder_path = os.path.join(post_quant_root, "Combined_by_perc")
	try:
		os.mkdir(folder_path)
	except FileExistsError: #* If the folder already exists, it will not be overwritten. The exception will be caught
		pass
	os.chdir(folder_path) #* Set the directory to the folder path where the files will output to

	y = Parallel(n_jobs = Global_variables['cpu_se'], verbose = 100)(delayed(Tag_percentage_run)(r=ci, percentile = Global_variables['percentiles'][0]) for ci in range(len(Chamber_index)))

	print(y)
	pd.Series(y).to_csv("output_selection.csv")

	# try: #* This is the same except it has a different output folder name. This is kept here just to know what the previous standard was
	# 	new_folder = os.join(microfluidics_results, "Combined_results")
	# 	os.mkdir(new_folder)
	# except FileExistsError:
	# 	new_folder = os.join(microfluidics_results, "Combined_results")
	# 	pass
	# os.chdir(new_folder)
	#* Keep a version of the df just to know what files were used. This will be helpful for merging new versions
	# Chamber_index.to_parquet("Chamber_index.parquet", index = False)

	# #, Actually copy the chosen files
	# y = Parallel(n_jobs=Global_variables['cpu_se'], verbose= 100)(delayed(Tag_percentage_run)(r = c) for c in Chamber_index['Chamber'].unique())



