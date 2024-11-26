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

#* Derived from "best_perc_from_AUG.py"

#%%

def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
	new = path.replace(os.sep, '/')
	return (new)

#* Made the version versitile so that it does not assume any locations. Also more stable this way
post_quant_root = slash_switch(input("where are the files stored?"))
post_quant_sub_99 = slash_switch(input("Where are the 99th percentile files? Folder name"))
post_quant_sub_95 = slash_switch(input("Where are the 95th percentile files? Folder name"))

#%%
measure_one = input("What are the values being compared? (This is the merged so it should be Percenage_reloc)")
# file_One = input("File One (95)") + ".csv"
# file_Two = input("File Two (99)") + ".csv"
file_One = slash_switch(input("File One (95)")) #* This is the _all file
file_Two = slash_switch(input("File Two (99)"))

#%%
file1_ninetyfive = pd.read_csv(file_One)
file1_ninetyfive = file1_ninetyfive.drop(columns = "Time_post_treatment")
file1_ninetyfive["Protein"]= file1_ninetyfive["Protein"].str.replace('-perc_t','')
file1_ninetyfive = file1_ninetyfive.set_index(["Protein", "Frames_post_treatment"])
file1_ninetyfive = file1_ninetyfive.add_suffix('_ninetyfive')
file1_ninetyfive.drop(file1_ninetyfive.columns[file1_ninetyfive.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)

file2_ninetynine = pd.read_csv(file_Two)
file2_ninetynine["Protein"]= file2_ninetynine["Protein"].str.replace('-perc_t','')
file2_ninetynine = file2_ninetynine.drop(columns='Time_post_treatment')
file2_ninetynine.drop(file2_ninetynine.columns[file2_ninetynine.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
file2_ninetynine = file2_ninetynine.set_index(["Protein", "Frames_post_treatment"])
file2_ninetynine= file2_ninetynine.add_suffix('_ninetynine')

master_temp = pd.merge(file1_ninetyfive, file2_ninetynine, left_index= True, right_index= True, how = 'outer') #* This is a new change 2023-09-04]]
master_temp.reset_index(level= 'Protein', inplace= True)
master_temp = master_temp.groupby("Protein").agg('max')

variable_col_95_one = measure_one + "_ninetyfive"
variable_col_99_one = measure_one + "_ninetynine"

# #Define the possible conditoins
# conditions = [master_temp[variable_col_95] > master_temp[variable_col_99],
#               master_temp[variable_col_95] < master_temp[variable_col_99]]

# #Define the output choices
# choices = ['ninetyfive', 'ninetynine']

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
measure_two = input("What are the values being compared?(This is the merged so it should be Percentage_reloc_less)")
# file_Three = input("File One (95)") + ".csv"
# file_Four = input("File Two (99)") + ".csv"
file_Three = slash_switch(input("File One (95)"))
file_Four = slash_switch(input("File Two (99)"))
#%%
file3_ninetyfive = pd.read_csv(file_Three)
file3_ninetyfive["Protein"]= file3_ninetyfive["Protein"].str.replace('-perc_yet','')
file3_ninetyfive = file3_ninetyfive.drop(columns = "Time_post_treatment")
file3_ninetyfive = file3_ninetyfive.set_index(["Protein", "Frames_post_treatment"])
file3_ninetyfive = file3_ninetyfive.add_suffix('_ninetyfive')
file3_ninetyfive.drop(file3_ninetyfive.columns[file3_ninetyfive.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)

file4_ninetynine = pd.read_csv(file_Four)
file4_ninetynine["Protein"]= file4_ninetynine["Protein"].str.replace('-perc_yet','')
file4_ninetynine = file4_ninetynine.drop(columns='Time_post_treatment')
file4_ninetynine.drop(file4_ninetynine.columns[file4_ninetynine.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
file4_ninetynine = file4_ninetynine.set_index(["Protein", "Frames_post_treatment"])
file4_ninetynine= file4_ninetynine.add_suffix('_ninetynine')

master_temp = pd.merge(file3_ninetyfive, file4_ninetynine, left_index= True, right_index= True)
master_temp.reset_index(level= 'Protein', inplace= True)
master_temp = master_temp.groupby("Protein").agg('max')

variable_col_95_two = measure_two + "_ninetyfive"
variable_col_99_two = measure_two + "_ninetynine"

#Define the possible conditoins
conditions = [master_temp[variable_col_95_two] > master_temp[variable_col_99_two],
              master_temp[variable_col_95_two] < master_temp[variable_col_99_two]]

#Define the output choices
choices = ['ninetyfive', 'ninetynine']

master = master_temp.copy()
#create new column with which has a larger percentage
master['larger(selected)'] = np.select(conditions, choices, default='ninetynine')
#%%
def result_give(val1, val2, selector):
	if selector == "ninetyfive":
		return(val1, selector)
	if selector == "ninetynine":
		return(val2, selector)

final_two = master.apply(lambda x: result_give(x[variable_col_95_two], x[variable_col_99_two], x["larger(selected)"]), axis = 1, result_type= 'expand').rename(columns= {0: measure_two, 1: "Selected_series"})

#%%
component_one = temp_one.copy()
component_two = final_two.copy()

first_combine = pd.merge(component_one, component_two, left_index= True, right_index=True, how = 'left')
final_one =  first_combine.apply(lambda x: result_give(x[variable_col_95_one], x[variable_col_99_one], x["Selected_series"]), axis = 1, result_type= 'expand').rename(columns= {0: measure_one, 1: "Selected_ser"})

final_combine= pd.merge(final_one, component_two, left_index= True, right_index=True, how = 'left')
final_graph_df = final_combine.sort_values(by = "Percentage_reloc_less").reset_index(drop= False)
final_graph_df.to_csv("Final_combined_comparison.csv")
#%%

final_graph_df["Selected_ser"] = final_graph_df["Selected_ser"] + "perc_f"

#%%
final_graph_df = pd.read_csv("Final_combined_comparison.csv")
measure_one = "Percentage_reloc"
measure_two = "Percentage_reloc_less"
import seaborn as sns
fig, ax = plt.subplots()
ax = sns.scatterplot(x = "Protein", y = measure_one, hue = 'Selected_series', data=final_graph_df)
ax1 = sns.scatterplot(x ="Protein", y = measure_two, hue = 'Selected_series', data=final_graph_df)
fig.savefig(f"Combined_penetrances3.svg")


#%%
import seaborn as sns
fig, ax = plt.subplots()
ax = sns.scatterplot(x = "Protein", y = measure_one, data=final_graph_df)
ax1 = sns.scatterplot(x = "Protein", y = measure_two, data=final_graph_df)
fig.savefig(f"Color_by_meas.svg")

#%%
final_graph_df = pd.read_csv("Final_combined_comparison.csv")
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
