#%% Load in the libraries
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
import numpy as np
import seaborn as sns
sns.set(rc={'figure.figsize':(11.7,8.27)})
#* Derived from "best_perc_from_AUG.py"

#....... This is a work in progress and is not variable. Not priority right now.
#Todo: complete the varaible function

def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
	new = path.replace(os.sep, '/')
	return (new)
#%% Define the global variables#%%
def join_percent_folder(percentiles, measure):
	count = 0
	if measure == 'Percentage_reloc':
		perc_name = 'All_pos_pt_percentages_melt.parquet'
	elif measure == 'Percentage_reloc_less':
		perc_name = 'All_pos_pt_t_less_percentages_melt.parquet'
	else:
		return('Invalid Measure')

	for p in percentiles:
		if count == 0:
			post_quant_perc_folder = f"{p}th_percentile"
			#Todo: make this so it can run for any percentages
			#, Read in the first pt_percentages file for 95th percentile
			temp_file = pd.read_parquet(slash_switch(
				os.path.join(post_quant_perc_folder,perc_name)))
			temp_file = temp_file.drop(columns = "Time_post_treatment")
			temp_file["Protein"]= temp_file["Protein"].str.replace('-perc_t','')
			temp_file = temp_file.set_index(["Protein", "Frames_post_treatment"])
			temp_file = temp_file.add_suffix('_ninetyfive')
			temp_file.drop(temp_file.columns[temp_file.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
			variable_cols = [f"{measure}_{p}"]
			master_temp = temp_file.copy()
		else:
			post_quant_perc_folder = f"{p}th_percentile"
			#Todo: make this so it can run for any percentages
			#, Read in the first pt_percentages file for 95th percentile
			temp_file = pd.read_parquet(slash_switch(
				os.path.join(post_quant_perc_folder,perc_name)))
			temp_file = temp_file.drop(columns = "Time_post_treatment")
			temp_file["Protein"]= temp_file["Protein"].str.replace('-perc_t','')
			temp_file = temp_file.set_index(["Protein", "Frames_post_treatment"])
			temp_file = temp_file.add_suffix('_ninetyfive')
			temp_file.drop(temp_file.columns[temp_file.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
			master_temp = pd.merge(master_temp, temp_file, left_index= True, right_index= True, how = 'outer') #* This is a new change 2023-09-04]]
			master_temp.reset_index(level= 'Protein', inplace= True)
			master_temp = master_temp.groupby("Protein").agg('max')
			variable_cols.append(f"{measure}_{p}")

		return(master_temp, variable_cols)


if __name__ == '__main__':
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': 'D:/MAD_1_5_comparision', #. gv.slash_switch(input("Post quant path?")) , #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': int(math.floor(os.cpu_count()*0.7)),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True}

	# log_prefix = datetime.datetime.today().strftime('%y_%m_%d-%H_%M_%S')
	log_prefix = datetime.today().strftime('%y_%m_%d-%H_%M_%S')

	#* Convert to local variables within the script
	microfluidic_results = Global_variables['microfluidics_results']
	post_path = Global_variables['post_path']
	list_percentiles = Global_variables['percentiles']
	time_per_frame = Global_variables['timepoint_gap']
	decade = str(datetime.now().year)[:3] #! I assume that this will not still be state of the art and running New Year's Eve 2029
	pn = Global_variables["cpu_se"] #* This can be changed
	today = str(date.today())
	year = str(date.today().year)
	decade = year[:3]

	post_path = Global_variables['post_path']
	os.chdir(post_path)

	perc_less = join_percent_folder(Global_variables['percentiles'], 'Percentage_reloc_less')

	perc_frame = join_percent_folder(Global_variables['percentiles'], 'Percentage_reloc')



#%%

variable_col_95_one = measure_one + "_ninetyfive"
variable_col_99_one = measure_one + "_ninetynine"

# #Define the possible conditoins
# conditions = [master_temp[variable_col_95] > master_temp[variable_col_99],
#               master_temp[variable_col_95] < master_temp[variable_col_99]]

# #Define the output choices
# choices = ['ninetyfive', 'ninetynine']

#* Make a copy of the dataframe to do manipulations on
temp_one = master_temp.copy()
#create new column with which has a larger percentage
# master['larger(selected)'] = np.select(conditions, choices, default='ninetynine')
#%%
# def result_give(val1, val2, selector):
# 	if selector == "ninetyfive":
# 		return(val1, selector)
# 	if selector == "ninetynine":
# 		return(val2, selector)

# final_one = master.apply(lambda x: result_give(x[variable_col_95], x[variable_col_99], x["larger(selected)"]), axis = 1, result_type= 'expand').rename(columns= {0: measure_one, 1: "Selected_series"})

#%%
# measure_two = input("What are the values being compared?(This is the merged so it should be Percentage_reloc_less)")
# file_Three = input("File One (95)") + ".parquet"
# file_Four = input("File Two (99)") + ".parquet"
# file_Three = slash_switch(input("File One (95)"))
# file_Four = slash_switch(input("File Two (99)"))
#%%
#, Read in the first pt_less file
file3_ninetyfive = pd.read_parquet(slash_switch(os.path.join(post_quant_sub_95,second_name)))
file3_ninetyfive["Protein"]= file3_ninetyfive["Protein"].str.replace('-perc_yet','')
file3_ninetyfive = file3_ninetyfive.drop(columns = "Time_post_treatment")
file3_ninetyfive = file3_ninetyfive.set_index(["Protein", "Frames_post_treatment"])
file3_ninetyfive = file3_ninetyfive.add_suffix('_ninetyfive')
file3_ninetyfive.drop(file3_ninetyfive.columns[file3_ninetyfive.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)

#, Read in the second pt_less file
file4_ninetynine = pd.read_parquet(slash_switch(os.path.join(post_quant_sub_99,second_name)))
file4_ninetynine["Protein"]= file4_ninetynine["Protein"].str.replace('-perc_yet','')
file4_ninetynine = file4_ninetynine.drop(columns='Time_post_treatment')
file4_ninetynine.drop(file4_ninetynine.columns[file4_ninetynine.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
file4_ninetynine = file4_ninetynine.set_index(["Protein", "Frames_post_treatment"])
file4_ninetynine= file4_ninetynine.add_suffix('_ninetynine')

#* Overwrite master_temp and do merging with the less files
master_temp = pd.merge(file3_ninetyfive, file4_ninetynine, left_index= True, right_index= True)
master_temp.reset_index(level= 'Protein', inplace= True)
master_temp['Prot_cop'] = master_temp['Protein']
master_temp = master_temp.groupby("Protein").agg('max')

variable_col_95_two = measure_two + "_ninetyfive"
variable_col_99_two = measure_two + "_ninetynine"

#Define the possible conditoins
conditions = [master_temp[variable_col_95_two] > master_temp		[variable_col_99_two],
            master_temp[variable_col_95_two] < master_temp[variable_col_99_two]]

#Define the output choices
choices = ['95th', '99th']

master = master_temp.copy()
#create new column with which has a larger percentage
master['larger(selected)'] = np.select(conditions, choices, default='99th') #* why did I set the default as 99?
#%%
def result_give_pt_less(val1, val2, selector, protein):
	if 'DDC2' in protein:
		return(val2, '99th')
	else:
		if selector == "95th":
			return(val1, selector)
		if selector == "99th":
			return(val2, selector)

def result_give_chosen(val1, val2, selector):
	if selector == "95th":
		return(val1, selector)
	if selector == "99th":
		return(val2, selector)


final_two = master.apply(lambda x: result_give_pt_less(x[variable_col_95_two], x[variable_col_99_two], x["larger(selected)"],x['Prot_cop']), axis = 1, result_type= 'expand').rename(columns= {0: measure_two, 1: "Selected_series"})

#%%
component_one = temp_one.copy()
component_two = final_two.copy()

first_combine = pd.merge(component_one, component_two, left_on= 'Protein', right_on = 'Protein', how = 'left')
final_one =  first_combine.apply(lambda x: result_give_chosen(x[variable_col_95_one], x[variable_col_99_one], x["Selected_series"]), axis = 1, result_type= 'expand').rename(columns= {0: measure_one, 1: "Selected_ser"})


final_combine= pd.merge(final_one, component_two, left_index= True, right_index=True, how = 'left')
final_graph_df = final_combine.sort_values(by = "Percentage_reloc_less").reset_index(drop= False)
final_graph_df = final_graph_df.loc[~(final_graph_df['Protein'].isin(['ATG16', 'NSP1', 'RAD50', 'AIM4', 'YMR031C', 'NOP56', 'YGR122W', 'PIL1']))]
final_graph_df.to_parquet("Final_combined_comparison.parquet")
#%%

final_graph_df["Selected_ser"] = final_graph_df["Selected_ser"] + "perc_f"

#%%
os.chdir(post_path)
final_graph_df = pd.read_parquet("Final_combined_comparison.parquet")
# input('NOTE THAT BELOW WILL REMOVE PROTEINS FROM d0212r1')

figures = "D:/Figures_root/2023-09-22"
# os.chdir(figures)

measure_one = "Percentage_reloc"
measure_two = "Percentage_reloc_less"
import seaborn as sns
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax = sns.scatterplot(x = "Protein", y = measure_one, hue = 'Selected_series', data=final_graph_df)
ax1 = sns.scatterplot(x ="Protein", y = measure_two, hue = 'Selected_series', data=final_graph_df)
fig.savefig(f"Combined_penetrances3.svg")
fig.savefig(f"Combined_penetrances3.pdf")


#%%
import seaborn as sns
fig, ax = plt.subplots()
ax = sns.scatterplot(x = "Protein", y = measure_one, data=final_graph_df)
ax1 = sns.scatterplot(x = "Protein", y = measure_two, data=final_graph_df)
fig.savefig(f"Color_by_meas.svg")
fig.savefig(f'Color_by_meas.pdf')

#%%
final_graph_df = pd.read_parquet("Final_combined_comparison.parquet")
measure_one = "Percentage_reloc"
measure_two = "Percentage_reloc_less"
#%%
Special_prots = pd.read_csv("Protein_list_special.csv", usecols= ["Gene > Standard Name"]).loc[:, "Gene > Standard Name"]

def Special_mark(Protein):
	if Protein in list(Special_prots):
		return(6)
	else:
		return(2)

final_graph_df["Special"] = pd.Series(final_graph_df["Protein"]).apply(Special_mark)
#%%

fig, ax = plt.subplots()
ax = sns.scatterplot(x = "Protein", y = measure_one, data=final_graph_df, size = "Special", sizes = (100, 300))
ax1 = sns.scatterplot(x = "Protein", y = measure_two, data=final_graph_df, size = "Special", sizes = (100, 300))
fig.savefig(f"Color_by_meas_better points.svg")
fig.savefig(f"Color_by_meas_better points.pdf")



#%%
import plotly.express as px
ref_one = px.scatter(x = "Protein", y = measure_one, data_frame=final_graph_df, size = "Special")
ref_two = px.scatter(x = "Protein", y = measure_two, data_frame=final_graph_df, size = "Special")

ref_one.write_html("REFERENCE_Frame_Penetrance.html")
ref_two.write_html("REFERENCE2_Timecourse_Penetrance.html")

#%%
import plotly.express as px
pen_id_refer = px.scatter(data_frame = final_graph_df,
		y= measure_two,
		x= 'Protein')
pen_id_refer.write_html('pen_id_refer.html')
#%%

#%%
fig, ax = plt.subplots()
sns.set_palette('dark')
ax = sns.scatterplot(x = "Protein", y = measure_one, hue = 'Selected_ser', data=final_graph_df, marker = "D")
sns.set_palette('bright')
ax1 = sns.scatterplot(x = "Protein", y = measure_two, hue = 'Selected_series', data=final_graph_df, marker= "o")
fig.savefig(f"Combined_penetrances_dark_light.svg")

#%%
final_graph_df_reloc_frame = final_combine.sort_values(by = "Percentage_reloc").reset_index(drop= False)


import seaborn as sns
fig, ax = plt.subplots()
ax = sns.scatterplot(x = "Protein", y = measure_one, hue = 'Selected_series', data=final_graph_df_reloc_frame, marker="D")
ax1 = sns.scatterplot(x = "Protein", y = measure_two, hue = 'Selected_series', data=final_graph_df_reloc_frame, marker = "o")



#%%
fig.savefig(f"Combined_penetrances_sort_reloc-frame.svg")


#%%
import plotly.express as px
m1 = px.scatter(final_graph_df, x = "Protein", y =measure_one)
m1.write_html("m1.html")
#%%
m2 = px.scatter(final_graph_df, x = "Protein", y =measure_two)
m2.write_html("m2.html")
#%%


# graph.write_image('test.svg')

#view the DataFra
# cols_two = file2.columns()

# sum = file1 - file2
# master = pd.DataFrame([])
# for col in cols_one:
# 	if


#%%
# import matplotlib.pyplot as plt

# x = range(100)
# y = range(100,200)
# fig = plt.figure()
# ax1 = fig.add_subplot(111)

# ax1.scatter(x[:4], y[:4], s=10, c='b', marker="s", label='first')
# ax1.scatter(x[40:],y[40:], s=10, c='r', marker="o", label='second')
# plt.legend(loc='upper left');
# plt.show()






# %%
