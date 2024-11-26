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

#. This version was derived from the quick versions found in 'current' which is no longer current
#! Quite a bit of stuff was removed from this version. The most recent version before this was the quick version

#%% Create a merged Chamber index table

def f_chamber(col):
	s = col.find("_Col")+1
	e = col.find("_" + year)
	Chamber = col[s:e]
	return(Chamber)

def run_perc(path): #* This will give the percentage that that it was run for. Folder names can be variable but shoudl have 'th' immediately following the percentage number
	start = path.find("th") - 2
	end = path.find("th")
	return(path[start:end])

def var_Chamber_index(post_quant_root, percentile):
	post_qunat_folder_percentile = f"{percentile}th_percentile" #* There is no need to account for dates, as this will be within a dated post_quant folder
	full_path = os.path.join(os.path.join(post_quant_root, post_qunat_folder_percentile))

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

	for root, dirs, files, in os.walk():
		for name in files:
			if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
				if name.endswith("index.parquet"):
					pass
				else:
					Path_n = f"Path_{percentile}"
					Chamber_index_95.append({Path_n: os.path.join(root, name)})
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

# #%% This was moved up and then retired. Appended index vs. merged
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


#%%This is for if the comaprison has already been completed
percentage_to_use = pd.read_parquet("Final_combined_comparison.parquet") #* This is fine to be hard coded as it will be produced just prior with fixed name from the percentile selection script
percentage_to_use_lib = percentage_to_use.loc[:, ["Protein", "Selected_series"]] #* Get a list of which series to be used

def conv_str_number(str_number): #* this is needed right now but will be fixed in next version
	if str_number == 'ninetynine':
		return(99)
	if str_number == 'ninetyfive':
		return(95)

percentage_to_use_lib["Selected_series"] = pd.Series(percentage_to_use_lib["Selected_series"]).apply(conv_str_number)
#%%
#* Why was this retired. This seems to be the merged type flow
# os.chdir(os.path.join(post_quant_root, post_quant_sub_99))
# All_chambers_l_99 = sorted(glob('Chamber_Col_*.parquet')) #* Get a unique list of the files, with thier roots

# os.chdir(os.path.join(post_quant_root, post_quant_sub_95))
# All_chambers_l_95 = sorted(glob('Chamber_Col_*.parquet'))

#%%

#? I don't know what this is?
# for p in percentage_to_use_lib["Protein"]:
# 	percentage = percentage_to_use_lib.loc[p, "Selected_series"]
# 	path = Chamber_index.loc[(Chamber_index["Protein"] == p) & (Chamber_index["Percent_run"]  == percentage)]


# for c in Chamber_index["Chamber"]:
# 	percentage = percentage_to_use_lib.loc[p, "Selected_series"]
# 	path = Chamber_index.loc[(Chamber_index["Protein"] == p) & (Chamber_index["Percent_run"]  == percentage)]
# 	chamber =

#%%
def Tag_percentage_run(r):
	try:
		row = Chamber_index.loc[r] #* For list of chamber values
		path = row["Path"]
		percent = int(row["Percent_run"])

		temp_df = pd.read_parquet(path)

		proteins = temp_df["Protein"].unique()

		try:
			protein_one = proteins[0]
			if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_one) & (percentage_to_use_lib["Selected_series"] == percent)]) >0:
				prot_one_df = temp_df.loc[temp_df["Protein"] == protein_one]

				prot_one_df.to_parquet(f"{protein_one}.parquet")
			else:
				pass
		except Exception as e:
			return(e, f"{path}")

		try:
			protein_two = proteins[1]
			if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_two) & (percentage_to_use_lib["Selected_series"] == percent)])>0:
				prot_two_df.to_parquet(f"{protein_two}.parquet")

				prot_two_df = temp_df.loc[temp_df["Protein"] == protein_two]
			else:
				pass
		except Exception as e:
			return(e, f"{path}")

		return(f"{path} complete")
	except Exception as e:
		return(e, f"{path}")

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

#%% Main Run stuff
if __name__ == '__main__':
	Global_variables = {'percentiles': [95,99],
					 	'cpu_se': 14} #* this is temporary and parameters for current experiment.
	#* Made the version versitile so that it does not assume any locations. Also more stable this way
	def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
		new = path.replace(os.sep, '/')
		return (new)
	today = str(date.today())


	post_quant_root = slash_switch(input("where are the files stored?"))

	#* Here if the path is given explicitly, run
	post_quant_sub_99 = slash_switch(input("Where are the 99th percentile files? Folder name"))
	post_quant_sub_95 = slash_switch(input("Where are the 95th percentile files? Folder name"))
	year = str(date.today().year)

	#. It is not critical to set the path here as it will be refered with full path in later calls
	os.chdir(post_quant_root) #* Temp set to the root path so that the information on which percentage to use can be read in

	#########################################< This the manual version
	Chamber_index_95 = []
	count = 0
	for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_95))):
		for name in files:
			if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
				if name.endswith("index.parquet"):
					pass
				else:
					Chamber_index_95.append({'Path_95': os.path.join(root, name)})
					count = count + 1
					print(count, end="\r")
			else:
				pass
		break #This makes the program run non recursively and not decend into daughter folders
	Chamber_index_95 = pd.DataFrame(Chamber_index_95)
	Chamber_index_95["Chamber"] = pd.Series(Chamber_index_95.iloc[:,0]).apply(f_chamber)
	Chamber_index_95["Percent_run"] = pd.Series(Chamber_index_95.iloc[:,0]).apply(run_perc)

	Chamber_index_99 = []
	count = 0
	for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_99))):
		for name in files:
			if name.startswith("Chamber") and name.endswith(f".parquet"): # fix naming
				if name.endswith("index.parquet"):
					pass
				else:
					Chamber_index_99.append({'Path_99': os.path.join(root, name)})
					count = count + 1
					print(count, end="\r")
			else:
				pass
		break #This makes the program run non recursively and not decend into daughter folders
	Chamber_index_99 = pd.DataFrame(Chamber_index_99)
	Chamber_index_99["Chamber"] = pd.Series(Chamber_index_99.iloc[:,0]).apply(f_chamber)
	Chamber_index_99["Percent_run"] = pd.Series(Chamber_index_99.iloc[:,0]).apply(run_perc)

	Chamber_index = pd.merge(Chamber_index_95,Chamber_index_99, left_on = 'Chamber', right_on = 'Chamber', how = 'outside')
	Chamber_index.to_parquet(f"Chamber_index_combined.parquet")
	#############################################>

	for p in Global_variables['percentiles']:
		p_n = 0
		if p_n == 0: #* This first call is just to deal with the first row, so that Chamber_index does not have to be defined outside
			Chamber_index= var_Chamber_index(post_quant_root = post_quant_root, percentile = p) #Todo: determine if post_quant root will be predefined
		else:
			Chamber_index_temp = var_Chamber_index(post_quant_root = post_quant_root, percentile = p) #Todo: determine if post_quant root will be predefined
			Chamber_index = pd.merge(Chamber_index, Chamber_index_temp, left_on = 'Chamber', right_on = 'Chamber', how = 'outside')
		p_n += 1

	Chamber_index_na = Chamber_index.isna()
	if len(Chamber_index_na) > 0:
		print(Chamber_index_na)
	state = input('Continue?')

	if state == 'y' or state == 'yes':
		pass
	else:
		print('Chamber index aysmetry and input cuased termination')
		exit

	#, Create folder to copy the chosen files into
	folder_path = os.path.join(post_quant_root, "Combined_by_perc")
	try:
		os.mkdir(folder_path)
	except FileExistsError: #* If the folder already exists, it will not be overwritten. The exception will be caught
		pass
	os.chdir(folder_path) #* Set the directory to the folder path where the files will output to

	# try: #* This is the same except it has a different output folder name. This is kept here just to know what the previous standard was
	# 	new_folder = os.join(microfluidics_results, "Combined_results")
	# 	os.mkdir(new_folder)
	# except FileExistsError:
	# 	new_folder = os.join(microfluidics_results, "Combined_results")
	# 	pass
	# os.chdir(new_folder)
	#* Keep a version of the df just to know what files were used. This will be helpful for merging new versions
	# Chamber_index.to_parquet("Chamber_index.parquet", index = False)

	#, Actually copy the chosen files
	y = Parallel(n_jobs=Global_variables['cpu_se'], verbose= 100)(delayed(Tag_percentage_run)(r = c) for c in Chamber_index['Chamber'].unique())
	print(y)
	pd.Series(y).to_parquet("output_selection.parquet") #Todo: This could be changed to text file


