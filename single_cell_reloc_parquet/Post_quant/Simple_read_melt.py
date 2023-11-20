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


#! Quite a bit of stuff was removed from this version. The most recent version before this was the quick version

#%%
today = str(date.today())

def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
    new = path.replace(os.sep, '/')
    return (new)

#%%
post_quant_root = slash_switch(input("where are the files stored?"))
post_quant_sub_99 = slash_switch(input("Where are the 99th percentile files? Folder name"))
post_quant_sub_95 = slash_switch(input("Where are the 95th percentile files? Folder name"))

os.chdir(post_quant_root) #* Temp set to the root path so that the information on which percentage to use can be read in
percentage_to_use = pd.read_csv("Final_combined_comparison.csv") #* This is fine to be hard coded as it will be produced just prior with fixed name
percentage_to_use_lib = percentage_to_use.loc[:, ["Protein", "Selected_series"]] #* Get a list of which series to be used

def conv_str_number(str_number): #* this is needed right now but will be fixed in next version
	if str_number == 'ninetynine':
		return(99)
	if str_number == 'ninetyfive':
		return(95)

percentage_to_use_lib["Selected_series"] = pd.Series(percentage_to_use_lib["Selected_series"]).apply(conv_str_number)

#%%

# os.chdir(os.path.join(post_quant_root, post_quant_sub_99))
# All_chambers_l_99 = sorted(glob('Chamber_Col_*.csv')) #* Get a unique list of the files, with thier roots

# os.chdir(os.path.join(post_quant_root, post_quant_sub_95))
# All_chambers_l_95 = sorted(glob('Chamber_Col_*.csv'))


#%%

#, Get the chambers and the percentage
Chamber_index = []
count = 0
for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_99))):
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
for root, dirs, files, in os.walk(os.path.join(os.path.join(post_quant_root, post_quant_sub_95))):
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
year = str(date.today().year)

def f_chamber(col):
	s = col.find("_Col")+1
	e = col.find("_" + year)
	Chamber = col[s:e]
	return(Chamber)

def run_perc(path): #* This will give the percentage that that it was run for. Folder names can be variable but shoudl have 'th' immediately following the percentage number
	start = path.find("th") - 2
	end = path.find("th")
	return(path[start:end])
Chamber_index["Chamber"] = pd.Series(Chamber_index.iloc[:,0]).apply(f_chamber)
Chamber_index["Percent_run"] = pd.Series(Chamber_index.iloc[:,0]).apply(run_perc) #*Extract the run perctage from the path. The folder will contain "xxth"
Chamber_index.to_csv(f"Chamber_index_combined.csv")

#%%
#* The below is to copy the files to a differnt folder.
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
# Chamber_index.to_csv("Chamber_index.csv", index = False)

#%%

# for p in percentage_to_use_lib["Protein"]:
# 	percentage = percentage_to_use_lib.loc[p, "Selected_series"]
# 	path = Chamber_index.loc[(Chamber_index["Protein"] == p) & (Chamber_index["Percent_run"]  == percentage)]


# for c in Chamber_index["Chamber"]:
# 	percentage = percentage_to_use_lib.loc[p, "Selected_series"]
# 	path = Chamber_index.loc[(Chamber_index["Protein"] == p) & (Chamber_index["Percent_run"]  == percentage)]
# 	chamber =

def Tag_percentage_run(r):
	try:
		row = Chamber_index.iloc[r]
		path = row["Path"]
		percent = int(row["Percent_run"])

		temp_df = pd.read_csv(path)
		proteins = temp_df["Protein"].unique()

		try:
			protein_one = proteins[0]
			if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_one) & (percentage_to_use_lib["Selected_series"] == percent)]) >0:
				prot_one_df = temp_df.loc[temp_df["Protein"] == protein_one]

				prot_one_df.to_csv(f"{protein_one}.csv")
			else:
				pass
		except Exception as e:
			return(e, f"{path}")

		try:
			protein_two = proteins[1]
			if len(percentage_to_use_lib.loc[(percentage_to_use_lib["Protein"] == protein_two) & (percentage_to_use_lib["Selected_series"] == percent)])>0:
				prot_two_df.to_csv(f"{protein_two}.csv")

				prot_two_df = temp_df.loc[temp_df["Protein"] == protein_two]
			else:
				pass
		except Exception as e:
			return(e, f"{path}")

		return(f"{path} complete")
	except Exception as e:
		return(e, f"{path}")

pr = 14

y = Parallel(n_jobs=pr, verbose= 100)(delayed(Tag_percentage_run)(r = c) for c in range(len(Chamber_index)))

print(y)
pd.Series(y).to_csv("output_selection.csv")
