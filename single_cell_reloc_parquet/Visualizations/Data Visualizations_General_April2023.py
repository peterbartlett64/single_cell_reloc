# %%
from datetime import date
# from turtle import color
from numpy.lib.function_base import median, quantile
import pandas as pd
import numpy as np
from pandas.core.algorithms import unique
from pandas.core.reshape.merge import merge
from pandas.io.parsers import read_csv
import scipy as scipy
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as scipyio
import glob
import os
import csv
import statistics
from scipy import stats
from sklearn.cluster import KMeans
import matplotlib
# from joblib import Parallel, delayed
from sklearn import preprocessing
import plotly.express as px
from glob import glob
import math
import kaleido
from joblib import Parallel, delayed
import seaborn as sns
import single_cell_reloc_parquet.global_functions.global_variables as glv

# sns.set(rc = {'figure.figsize':(15,8)}) #* Set the figure size to be larger
sns.set(rc = {'figure.figsize':(40,30)}) #* Set the figure size to be the shape of a powerpoint slide

#%%
#TODO: Add the protien/loc grouping for the penetrance plots
#TODO: Further expore the z-score correlation
#? This has partially been done, but I will continue with it

# %%

if __name__ == "__main__":
	microfluidic_results = str(input("Microfluidics results folder"))
	post_path= str(input("Post_quant results folder"))
	# Global_Variables = glv.global_manager()
	microfluidic_results = glv.slash_switch(microfluidic_results)
	post_path = glv.slash_switch(post_path)
	os.chdir(post_path)

else:
	try:
		Global_Variables #* If this is not the main file, then the Global_Variables should already be assigned
	except:
		microfluidic_results = str(input("Microfluidics results folder"))
		post_path= str(input("Post_quant results folder"))
		Global_Variables = glv.Global_variables

figures_root = "C:/Users/pcnba/Dropbox (Grant Brown's Lab)/Peter Bartlett Data/Figures/Direct_output" #! Replace this with a flexible version for final publish
date_today = str(date.today)
figures = os.path.join(figures_root, date_today) #* This was recenlty added to make sure that the figues are not overwriting
time_gap = 7.5

# Global_Variables[timepoint_gap]


#%%
# full_data_name = "Subset_concat.csv"
x=0
while x == 0:
	f_type = input("csv or parquet?")
	full_data_name = input("FileName")

	if f_type == "csv":
		# full_data = pd.read_csv(full_data_name+".csv")
		x=1
	if f_type == "parquet":
		# full_data = pd.read_parquet(full_data_name+".parquet")
		x=1
	else:
		pass
#%%
#, Pearson Correlations
if f_type == "csv":
	Abundance_df = pd.read_csv(full_data_name, usecols= ["z_score_Abund", "z_score_Loc", "Protein", "Pearson_coeff_Loc_Abund"]) # This way only what is needed is read in.
else:
	Abundance_df = pd.read_parquet(full_data_name, usecols= ["z_score_Abund", "z_score_Loc", "Protein", "Pearson_coeff_Loc_Abund"])


#, Abundance-Loc_pearson_r
median_pearson = Abundance_df.groupby(["Protein"])["Pearson_coeff_Loc_Abund"].agg(median)
pearson_coeff_bar = px.bar(median_pearson.sort_values())
os.chdir(figures)
pearson_coeff_bar.write_html(f'Pearson_r_bar.html')
# pearson_coeff_bar =
#%%
#, All_loc_abundance
all_pear_scat = px.scatter(data_frame= full_data, x = "z_score_Loc", y = "z_score_Abund", color = "Protein")

os.chdir(figures) #* Change the path to the figures for direct output
all_pear_scat.write_html(f'Pearson_scat_all.html')

# %%


# test = full_data[full_data["Protein"] == "ARO4"].drop_duplicates()
# px.scatter(data_frame= test, x = "Loc_score", y = "Abundance", color="Frame_x")


# %%
# full_data = pd.read_csv(f"ALL_CHAMBERS_COMBINED.csv")
# full_data = pd.read_csv(f"ALL_CHAMBERS_COMBINED{percentile_read}.csv") #*Read in the percentile
full_data = pd.read_csv(full_data_name)

#%%
#, This is to produce the automated pearson outputs
os.chdir(post_path)
post_treatment_df = full_data.loc[full_data["Frames_post_treatment"] > 0]
Agg_loc = post_treatment_df.groupby(["Protein", "Frame_x"])["Loc_score"].median() #* First get the median/mean of Loc_scores in each protein frame
Agg_loc = pd.DataFrame(Agg_loc)
Agg_loc.reset_index(drop = False, inplace = True)
max_agg_Loc_score = Agg_loc.loc[Agg_loc.groupby(["Protein"])["Loc_score"].idxmax()][["Protein", "Frame_x"]] #*Create a table of maximum aggregated Loc_scores for each protein

#%%
max_agg_Loc_score =  max_agg_Loc_score
#! This is to get rid of shit rows
max_agg_Loc_score = max_agg_Loc_score.iloc[3:].reset_index(drop = True)

#%%
functions_list = ["single_pearson"] #* This is for the function list implementation where they can be selected from a list

func_prot_combo = input(f"Single visualizations. Are there any special protiens that you want to look at? Write the function name from: {functions_list} and the protein of interest. Results will be output the direct ouput folder. If none, write [DONE]")

if func_prot_combo== "DONE":
	pass
else:
	combo_split =  str.split(func_prot_combo)
	func_run = combo_split[0]
	prot_interest = combo_split[1]
	exec(f"{func_run}({prot_interest})") #* This is sus but will work for now until I remember the better way to do

# print("Program is done, but everything is still in memory")

#, Single Prearson Outputs. This is for a dot plot
# prots = full_data["Protein"].unique()
def graph_pearson_single(i, full_data):
	prot = max_agg_Loc_score.iloc[i,:]["Protein"]
	f = max_agg_Loc_score.iloc[i,:]["Frame_x"]
	try:
		sub_data = full_data.loc[(full_data["Protein"] == prot) & (full_data["Frame_x"] == f)].dropna(subset=['z_score_Loc', 'z_score_Abund'])
	except ValueError:
		sub_data["z_score_Abund"] = stats.zscore(sub_data.loc[:,"Abundance"])
		sub_data["z_score_Loc"] = stats.zscore(sub_data.loc[:,"Loc_score"])
	z_Loc_Abun_r = stats.pearsonr(sub_data["z_score_Loc"], sub_data["z_score_Abund"])
	scatter = sns.scatterplot(data= sub_data, x = "z_score_Abund", y = "z_score_Loc")
	fig = scatter.get_figure()
	fig.savefig(f"Pearson{z_Loc_Abun_r}--{prot}_{f}.svg")
	print(f"{prot} is done pearson'n")
	return(prot, z_Loc_Abun_r)

def graph_pearson_plus_minus(i, full_data):
	try:
		prot = max_agg_Loc_score.iloc[i,:]["Protein"]
		f = max_agg_Loc_score.iloc[i,:]["Frame_x"]
		try:
			sub_data = full_data.loc[(full_data["Protein"] == prot) & ((full_data["Frame_x"] == f)| (full_data["Frame_x"] == f+1)| (full_data["Frame_x"] == f-1))].dropna(subset=['z_score_Loc', 'z_score_Abund'])
		except ValueError:
			# sub_data["z_score_Abund"] = stats.zscore(sub_data.loc[:,"Abundance"])
			# sub_data["z_score_Loc"] = stats.zscore(sub_data.loc[:,"Loc_score"])
			sub_data["z_score_Abund"] = sub_data.groupby("Frame_x")["Abundance"].transform(stats.zscore)
			sub_data["z_score_Loc"] = sub_data.groupby("Frame_x")["Loc_score"].transform(stats.zscore)
		z_Loc_Abun_r = stats.pearsonr(sub_data["z_score_Loc"], sub_data["z_score_Abund"])
		scatter = sns.scatterplot(data= sub_data, x = "z_score_Abund", y = "z_score_Loc", hue = "Frame_x")
		fig = scatter.get_figure()
		fig.savefig(f"Pearson{z_Loc_Abun_r}--{prot}_{f}.svg")
		print(f"{prot} is done pearson'n")
		return(prot, z_Loc_Abun_r)
	except Exception as e:
		return(f"{e} on {prot}")

pr = os.cpu_count()
# smaller = full_data.loc[full_data["Protein"] == "ECO1"]
x = Parallel(n_jobs= pr, verbose= 100)(delayed(graph_pearson_plus_minus)(i = i, full_data = smaller) for i in range(len(max_agg_Loc_score)))
print(x)

# %%
# def U_Frame_to_Pos(f):
# 	e = f.find("f")
# 	return(f[:e])

# full_data["Unique_pos"] = pd.Series(full_data["Unique_Frame"]).apply(U_Frame_to_Pos)
# full_data

#%%
# All_proteins = pd.read_csv("ALL_CHAMBERS_FINAL.csv")
# All_proteins = pd.read_csv(f"ALL_CHAMBERS_{percentile_read}.csv")
# All_proteins = pd.read_csv("D:/Microfluidics/RESULTS_ALL/Combined_99_PrePoster/ALL_CHAMBERS.csv")
# All_proteins = pd.read_csv("Chamber_Col_d0222_r1.0_Ch120.0_2022-03-08.csv")
#%%
# All_proteins = pd.read_csv("C:/Users/pcnba/OneDrive/Desktop/Test/ALL_CHAMBERS_95.csv")
All_proteins = full_data
Special_prots = pd.read_csv("Protein_list_special.csv")
Special_prots.rename(columns ={"Gene > Standard Name": "Gene"}, inplace= True)
Special_prots = Special_prots["Gene"]

#%%
Var_ap = All_proteins[["<lambda_0>_y", "Protein", "Frames_post_treatment"]].drop_duplicates().reset_index(drop = True)

##%%
Var_ap = Var_ap[(Var_ap["Protein"] == "RRP5") | (Var_ap["Protein"] == "FLR1")]
Var_ap["CoV_apos"] =  Var_ap["<lambda_0>_y"]
Var_ap = Var_ap.loc[Var_ap["Frames_post_treatment"] >= 0]
Var_ap["Time_post_treatment"] = Var_ap["Frames_post_treatment"] * 7.5
# Var_ap = Var_ap[Var_ap["Frames_post_treatment"] <= 40]
#%%
Var_scatter = sns.scatterplot(data= Var_ap, x = "Time_post_treatment", y = "CoV_apos", hue = "Protein")
grouped_Prot = Var_ap.groupby("Protein")
Var_med = grouped_Prot.agg('median')
Var_med.reset_index(inplace=True)
MRT4_CoV_median = Var_med[Var_med["Protein"] == "RRP5"]["CoV_apos"].values[0]
FLR1_CoV_median = Var_med[Var_med["Protein"] == "FLR1"]["CoV_apos"].values[0]

Var_scatter.axhline(y = MRT4_CoV_median, color = 'blue')
Var_scatter.axhline(y = FLR1_CoV_median, color = 'orange')
Var_scatter_fig = Var_scatter.get_figure()
Var_scatter_fig.savefig("CoV_Prots.svg")

#%%
grouped_Prot = Var_ap.groupby("Protein")
Var_max = grouped_Prot.agg('max')
Var_max.reset_index(inplace = True)

Var_bar  = px.bar(data_frame = Var_max,
		y= "Protein",
		x= 'CoV_apos',
		# color = 'Frames_post_treatment',
		orientation= 'h'
		).update_yaxes(categoryorder='total descending')
Var_bar.show()

#%%
def f_pos_frame(f):
	e = f.find("f")
	pos = f[:e]
	return(pos)

def reloc_yet(x):
    ind = x.idxmax()
    does = x.loc[ind]
    x.loc[:ind] = 0
    x.loc[ind:] = does
    return(x)
def workaround(ind):
    return(Post.loc[ind, "Frame_x"])

Post = full_data.loc[full_data["Frames_post_treatment"] >= 0].copy().reset_index()
# Post = All_proteins.loc[(All_proteins["Protein"] == "ECO1") | (All_proteins["Protein"] == "XRS2")].copy()
# Post = Post.loc[Post["Frames_post_treatment"]>=0].copy().reset_index()
#%%
# Post = Post.dropna(subset=['Relocalized'])

# Post.dropna(subset=["Relocalized"])#! ADDED AUG 11
# Post["n"] = Post.groupby(by = "Cell_Barcode")["Relocalized"].transform(reloc_yet)
# Post["ind"] = Post.groupby(by = "Cell_Barcode")["Relocalized"].transform('idxmax')
# Post["first"] = pd.Series(Post["ind"]).apply(workaround)

Post["Reloc_yet"] = Post.groupby(by= "Cell_Barcode", group_keys=False)["Relocalized"].apply(reloc_yet)  #! changed from transform(''). Check to make srue it


# Post["Position"] = pd.Series(Post["ImageID"]).apply(f_pos_frame)
#%%
Post["Time_post_treatment"] = Post["Frames_post_treatment"] * time_gap

# Post = Post[Post["Frames_post_treatment"] >= 0] # Just values after treatment
# Post.to_csv("ALL_CHAMBERS.csv", index = False)

# def first_time(ser):
#     first = ser.idmax()
#     value = ser.max()
#     return(first, value)

Post["Mover"] = Post.groupby(by="Cell_Barcode")["Relocalized"].transform('max')   # This mover function should be replaced
# does = pd.DataFrame(does)
# does.reset_index(inplace=True)
# does.set_index(["Cell_Barcode", "ImageID"])
# All_proteins.set_index(["Cell_Barcode", "ImageID"])
# All_proteins = pd.merge(All_proteins, does, left_index=True, right_index = True)
# has = All_proteins.groupby(by="Cell_Barcode")
# first = All_proteins.groupby(by=["Cell_Barcode"])["Relocalized"]


# grp = All_proteins.groupby(by = "Cell_Barcode")["Relocalized"].first(1)
# first_values = grp.reset_index()
# first_valuesper

# %%
import seaborn as sns

# sns.set(rc = {'figure.figsize':(15,8)})
Protein = str(input("Protein?"))
# Protein = p
single_protein = Post.loc[Post["Protein"] == Protein]
# times = [0, 165, 315]
times = np.arange(start=0, stop=1000, step=(4*7.5))
single_protein = single_protein.loc[single_protein["Time_post_treatment"].isin(times)]
ax = sns.stripplot(x = "Time_post_treatment", y= "Loc_score", data = single_protein, hue = "n")#hue = "Mover")
# ax = sns.swarmplot(x="Time_post_treatment", y="Loc_score", data=single_protein, color="black", edgecolor="gray")
ax.axhline(y = 1)
fig = ax.get_figure()
fig.savefig(f"{Protein}_strip_testing.svg", dpi = 300)
del fig
del ax



# %%
plot2 = sns.violinplot(x = "Time_post_treatment", y= "Loc_score", data = single_protein, hue= "Mover", innner = None)
plot2.axhline(y = 1)
fig = plot2.get_figure()
fig.savefig(f"{Protein}_third_Time_moversplit.svg", dpi = 300)

# %%
plot3 = sns.swarmplot(x = "Time_post_treatment", y= "Loc_score", data = single_protein)
plot3.axhline(y = 1)
fig = plot3.get_figure()
# fig.savefig(f"{Protein}_third_Time_moversplit.svg", dpi = 300)

# %%
Protein

# %%
# import seaborn as sns
times = np.arange(start=-22.5, stop=1000, step=(4*7.5))
for Protein in All_proteins["Protein"].unique():
	single_protein = All_proteins.loc[(All_proteins["Protein"] == Protein)]
	single_protein = single_protein.loc[single_protein["Time_post_treatment"].isin(times)]
	# ax = sns.violinplot(x = "Time_post_treatment", y= "Loc_score", data = single_protein, inner = None)
	# ax = sns.swarmplot(x="Time_post_treatment", y="Loc_score", data=single_protein, color="black", edgecolor="gray")
	# ax.axhline(y = 1)
	# fig = ax.get_figure()
	# fig.savefig(f"{Protein}_thirdTime_moversplit.svg", dpi = 300)

	plot2 = sns.violinplot(x = "Time_post_treatment", y= "Loc_score", data = single_protein, hue= "Mover", split = True)
	plot2.axhline(y = 1)
	fig = plot2.get_figure()
	fig.savefig(f"{Protein}_fourth_Time_moversplit.png", dpi = 300)


# %%
import plotly.graph_objects as go

# Protein = str(input("Protein?"))

list_p= All_proteins["Protein"].unique()
for Protein in list_p:
	# df = pd.read_csv("d0222r1p070300_movement_pseudo.csv")
	# df = Quant_ALL
	Time_violin_static = go.Figure()

	df = All_proteins[All_proteins["Protein"] == Protein]
	Frames_post_treatment = df["Time_post_treatment"].drop_duplicates()
	# Protein = "HGH1"

	cell_count = len(df["Cell_Barcode"].unique())

	for post_treat_t in Frames_post_treatment:
		Time_violin_static.add_trace(go.Violin(x = df["Time_post_treatment"][df["Time_post_treatment"] == post_treat_t],
																					y = df["Loc_score"][df['Time_post_treatment'] == post_treat_t], #["Loc_score"],
																					name = post_treat_t,
																					box_visible = True,
																					# line_color = 'black',
																					showlegend = False))
	thresh = 1
	Time_violin_static.add_shape(go.layout.Shape(type='line', xref='paper', yref='y',
																				x0=0, y0=thresh, x1=1, y1=thresh, line={'dash': 'dash'}))
	# Time_violin_static_px.add_shape(go.layout.Shape(type='line', xref='x', yref='paper',
	# 																			x0=0, y0=0, x1=0, y1=1, line={'dash': 'dash'}))
	Time_violin_static.update_layout(
		title={
			'text':  f"{Protein} <br><sup>{cell_count} cells </sup>",
			'y':0.9,
			'x':0.5,
			'xanchor': 'center',
			'yanchor': 'top'},
		xaxis_title = "Time post treatment (minutes)",
		yaxis_title = "Loc_score (Proprotion; Horizonal represents baseline)",
		xaxis = dict(
			tickmode = 'linear',
			dtick = 14
    	)
	)
	Time_violin_static.show()
	Time_violin_static.write_html(f'Time_violin_static_{Protein}.html')
	# Time_violin_static.write_image(f"Time_violin_static_{Protein}_2.pdf")


# %%

	Time_violin_static.write_html(f'Time_violin_static_{Protein}+2.html')
	Time_violin_static.write_image(f"Time_violin_static_{Protein}_2.pdf")

# %%
data["Protein"].unique()

# %%
max = 40
min = -1
mid = (max-min)//2

g = [max, min, mid]
g

# %%
import plotly.express as px

# df = pd.read_csv("d0222r1p070300_movement_pseudo.csv")
# df = Quant_ALL

# Frames_post_treatment = df["Frames_post_treatment"].drop_duplicates()
# Protein = "NOP58"
# df = df[df["Protein"] == Protein1]
Protein = str(input("Protein?"))

df = All_proteins[All_proteins["Protein"] == Protein]

#Show max, min, and middle
max = df["Frames_post_treatment"].max()
min = -1
# midmid = (max-min)//2
midlow = (max-min)//3
midhigh  = ((max-min)//3)*2
tmin = min *7.5

g = [max, min, midlow, midhigh]
g

df = df[df["Frames_post_treatment"].isin(g)]

# cell_count = len(df["Cell_Barcode"].unique())
thresh = 1
Time_violin_static_px = px.violin(df,
													x = "Time_post_treatment",
													y = "Loc_score",
													box = True,
													points= 'all',
													title = Protein)
Time_violin_static_px.add_hline(y  = thresh,  line_dash="dash")
# Time_violin_static_px.update_xaxes(tick0 = min*7.5, dtick = 15)
# Time_violin_static_px.add_hline(y = 0.8, line_dash = "dash", line_color="red")
# Time_violin_static_px.add_vline(x = 0)
Time_violin_static_px.show()
# os.chdir(figures)
# Time_violin_static_px.write_html(f"{Protein}_viostatic_wpoints_wbox.html")
Time_violin_static_px.write_image(f"{Protein}_violin_series.svg")

# %%
len(df["Cell_Barcode"].unique())

# %%

# Upgrade to all files.

# # Protein = "NOP58"
df = df[df["Protein"] == Protein1]

thresh = 1.15
Scartter_var_Loc = px.scatter(df, animation_frame = "Frames_post_treatment", x = 'CoV_apos', y = "x95thPercentile_norm_OBJ_Median_GFP", title = Protein1)
Time_violin_static_px.add_hline(y  = thresh,  line_dash="dash")
# Time_violin_static_px.add_vline(x = 0)
Time_violin_static_px.show()

# %%
All_proteins["Protein"]

# %%
All_proteins[All_proteins["Protein"] == 'ECD2']

# %%
Proteins = str(input("Proteins?"))
Proteins = Proteins.split(", ")EDC2


import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'notebook_connected'

# This format takes in a single postion. Coudl change color = "Protein" to color = "Postion" (would need to add)

# df = pd.read_csv("d0222r1p070300_movement_pseudo.csv")
# df = Quant_ALL

# Protein = str(input("Protein?"))

list_p= All_proteins["Protein"].unique()
df = All_proteins[All_proteins["Protein"].isin(Proteins)]
df.sort_values(by = "Frame_x", ascending = True, inplace = True)
# Frames_post_treatment = df["Frames_post_treatment"].drop_duplicates()
# # Protein = "NOP58"
# df = df[df["Protein"] == Protein]

thresh = 1
Time_violin_anim = px.violin(df,
													x = "Protein",
													y = "Loc_score",
													# y = "x95thPercentile_norm_OBJ_Median_GFP",
													points='all',
													animation_frame= "Frames_post_treatment",
													range_y = [0.5,1.5], # Create auto_range in next version of graphing
													color = "Cell_Barcode")
													# hover_data = ["Cell_Barcode", "Protein"]) # This was added to better track the unique cell for the bi-distribution to check that cells are moving/staying in one distribution
 													# color="Postion") # The position is not in the dataframe which is being referecenced right now
Time_violin_anim.add_hline(y  = thresh,  line_dash="dash")
Time_violin_anim.show(renderer = 'notebook_connected')
Time_violin_anim.write_html(f"Time_violin_{Proteins}anim_test.html")

# %%
# df.sort_values(by = 'Frame_x')
df[(df["Cell_Barcode"] == 'd0222r1p340200c0625')]

# %%
All_proteins

# %%
# Proteins = str(input("Proteins?"))
# Proteins = Proteins.split(", ")

import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'notebook_connected'

# This format takes in a single postion. Coudl change color = "Protein" to color = "Postion" (would need to add)

# df = pd.read_csv("d0222r1p070300_movement_pseudo.csv")
# df = Quant_ALL


Protein = str(input("Protein?"))

# list_p= All_proteins["Protein"].unique()
df = All_proteins[All_proteins["Protein"] == Protein]

# Frames_post_treatment = df["Frames_post_treatment"].drop_duplicates()
# # Protein = "NOP58"
# df = df[df["Protein"] == Protein]

df.sort_values(by = "Frames_post_treatment", inplace = True)
thresh = 1
Time_violin_anim = px.violin(df,
													x = "Protein",
													y = "Loc_score",
													# y = "x95thPercentile_norm_OBJ_Median_GFP",
													points='all',
													animation_frame= "Time_post_treatment",
													range_y = [0.5,1.8],
													hover_data = ["Cell_Barcode", "Protein"],
													color = "Mover")
 													# color="Postion") # The position is not in the dataframe which is being referecenced right now
Time_violin_anim.add_hline(y  = thresh,  line_dash="dash")
Time_violin_anim.show(renderer = 'notebook_connected')
Time_violin_anim.write_html(f"Time_violin_new_{Protein}anim.html")

# %%
del figures

# %%
All_proteins[All_proteins["Protein"] == 'MRT4']

# %%
MRT_non = All_proteins[All_proteins["Protein"] == 'MRT4']
# )& (All_proteins[All_proteins["Mover"] == 0])]
MRT_non = MRT_non[(MRT_non["Mover"] == 0) & (MRT_non["Frame_x"] == 69)]["Cell_Barcode"].unique()
ls = pd.DataFrame(MRT_non)
ls.sort_values(by = 0, inplace = True)
ls.set_index(0, drop = False)



# %%
ls.loc["d0222r1p340200c0273"]

# %%

Time_violin_anim.write_html(f"Time_violin_new_{Protein}anim.html")

# %%
# Proteins = str(input("Proteins?"))
# instances = Proteins.split(", ")

import plotly.express as px

df = pd.read_csv("ALL_CHAMBERS.csv")
Time_violin_anim = px.violin(df,
													x = "Protein",
													y = "Loc_score",
													# y = "x95thPercentile_norm_OBJ_Median_GFP",
													points='all',
													facet_col= 4,
													animation_frame= "Frames_post_treatment",
													range_y = [0,2],
													color = "Protein")
													# color="Postion") # The position is not in the dataframe which is being referecenced right now
Time_violin_anim.add_hline(y  = thresh,  line_dash="dash")
Time_violin_anim.show(renderer = 'notebook_connected')
Time_violin_anim.write_html(f"Time_violin_{names}anim.html")

# %%
# names= str(input("Proteins?"))
del Time_scatter_anim

df = pd.read_csv("ALL_CHAMBERS.csv")
Time_scatter_anim = px.scatter(df,
													y = "Loc_score",
													x = "x95thPercentile_norm_OBJ_Median_GFP",
													animation_frame= "Frames_post_treatment",
													range_y = [0,2],
													size = "CoV_apos",
													 hover_data=['Protein'],
													color = "Protein")
													# color="Postion") # The position is not in the dataframe which is being referecenced right now
Time_scatter_anim.show(renderer = 'notebook_connected')
Time_scatter_anim.write_html(f"Scatter_{names}anim.html")

# %%
import plotly.express as px
import plotly.io as pio
pio.renderers.default = 'notebook_connected'

# This format takes in a single postion. Coudl change color = "Protein" to color = "Postion" (would need to add)

df = pd.read_csv("d0222r1p070300_movement_pseudo.csv")

# Frames_post_treatment = df["Frames_post_treatment"].drop_duplicates()
# # Protein = "NOP58"
# # df = df[df["Protein"] == Protein]

thresh = 1.15
Time_violin_anim = px.scatter(df,
													x = "x99thPercentile_Diff_background_mKO", # Want to change it to
													y = "x95thPercentile_norm_OBJ_Median_GFP",
													# size=  "CoV"
													animation_frame= "Frames_post_treatment",
													color = "Protein")
													# color="Postion") # The position is not in the dataframe which is being referecenced right now
# Time_violin_static_px.add_hline(y  = thresh,  line_dash="dash")
Time_violin_anim.show(renderer = 'notebook_connected')
Time_violin_anim.write_html("Time_violin_anim.html")

# %%
os.getcwd()

# %%
#MULTI PLOT


# Quant_prim_index = pd.read_csv("Quant_prim_index.csv")
# pos_list = Quant_prim_index["PositionID"]
# fig, axes = plt.subplots(6,4,figsize=(210,297)) # This is the size of a piece of paper in mm

# # Template for hue overlay of all types for each postion. There are more postions than there are proteins, so it may be worth combining
# col = 0
# row = 0
# for p in pos_list:
# 	Pos_data = pd.read_csv(f"Quant_{p}_primary.csv")
# 	sns.lineplot(ax=axes[row,col],
# 			data = Pos_data,
# 			x='Frames_post_treatment',
# 			y='x95thPercentile_norm_BKGRND_Median_GFP',
# 			hue='Protein',
# 			ci=99).set_title(p)
# 	if col <= 2: # This is 2 because index starts at 0! Don't forget
# 		col += 1
# 	else:
# 		col = 0
# 		row += 1

# # plt.subplots_adjust(wspace=1, hspace=1)

# %%
## MULTIPLOT

# Quant_prim_index = pd.read_csv("Quant_prim_index.csv")
# pos_list = Quant_prim_index["PositionID"].values[80:160]
# fig, axes = plt.subplots(10,8,figsize=(210,297)) # This is the size of a piece of paper in mm

# # Template for hue overlay of all types for each postion. There are more postions than there are proteins, so it may be worth combining
# col = 0
# row = 0
# for p in pos_list:
# 	Pos_data = pd.read_csv(f"Quant_{p}_primary.csv")
# 	sns.lineplot(ax=axes[row,col],
# 			data = Pos_data,
# 			x='Frames_post_treatment',
# 			y='x95thPercentile_norm_BKGRND_Median_GFP',
# 			hue='Myo1Identity',
# 			ci=99).set_title(p)
# 	if col <= 6:
# 		col += 1
# 	else:
# 		col = 0
# 		row += 1

# plt.subplots_adjust(wspace=1, hspace=1)
# plt.savefig(f"8_by_10_AvgTrajectories")

# %%
Quant_prim_index = pd.read_csv("Quant_prim_index.csv")
pos_list = Quant_prim_index["Position"]
fig, axes = plt.subplots(13,13,figsize=(100,100))

for p in pos_list:
	Pos_data = pd.read_csv(f"<Unif_df>{p}.csv")
	Prots = Pos_Data["Protein"].unique().dropna()
	for fig, prot in zip(axes, Prots):
		Prot_Pos_data = Pos_data[Pos_data["Protein"].str.match(prot)]
		sub.lineplot(x=Prot_Pos_data['Frames_post_treatment'],
				y=Prot_Pos_data['x95thPercentile_norm_BKGRND_Median_GFP'],
				hue=Prot_Pos_data['Myo1Identity'],
				ci=99)
		# sub.scatter(x=Prot_Pos_data['Frames_post_treatment'], y=df['x95thPercentile_norm_BKGRND_Median_GFP'], s = 1)
		sub.set_title(prot)


# %%


# %%
os.chdir(microfluidics_results)
Unif_mKa = pd.read_csv(f"unif_mKa_{p}.csv")
Unif_mKO = pd.read_csv(f"unif_mKO_{p}.csv")

MMS_mKa = Unif_mKa
#pd.read_csv(f"MMS_{p}mKa_wReloc.csv")
MMS_mKO = Unif_mKO
#d.read_csv(f"MMS_{p}mKO_wReloc.csv")
#%%
MMS_relocalized_mKa_maxes = pd.DataFrame(MMS_mKa[MMS_mKa["Relocalized"] == 1].groupby(["Cell_Barcode"], sort = False)[["x95thPercentile_norm_BKGRND_Median_GFP", "Protein", "Relocalized"]].max())

MMS_relocalized_mKO_maxes = pd.DataFrame(MMS_mKO[MMS_mKO["Relocalized"] == 1].groupby(["Cell_Barcode"], sort = False)[["x95thPercentile_norm_BKGRND_Median_GFP", "Protein",  "Relocalized"]].max())

MMS_nomove_mKa_maxes = pd.DataFrame(MMS_mKa[MMS_mKa["Relocalized"] == 0].groupby(["Cell_Barcode"], sort = False)[["x95thPercentile_norm_BKGRND_Median_GFP", "Protein",  "Relocalized"]].max())

MMS_nomove_mKO_maxes = pd.DataFrame(MMS_mKO[MMS_mKO["Relocalized"] == 0].groupby(["Cell_Barcode"], sort = False)[["x95thPercentile_norm_BKGRND_Median_GFP", "Protein",  "Relocalized"]].max())


MMS_maxes_conc_mKa = pd.concat([MMS_relocalized_mKa_maxes, MMS_nomove_mKa_maxes]).reset_index()
MMS_maxes_conc_mKO = pd.concat([MMS_relocalized_mKO_maxes, MMS_nomove_mKO_maxes]).reset_index()

# %%
MMS_mKa

# %%
MMS_maxes_conc_mKa

# %%

def Cell_Bar_TO_upos(cb):
	end = cb.find("c")
	return (cb[0:end])



# MMS_maxes_conc_mKa["Unique_pos"] = MMS_maxes_conc_mKa.iloc[0]["Cell_Barcode"][0:
# MMS_maxes_conc_mKa.iloc[0]["Cell_Barcode"].find("c")]
MMS_maxes_conc_mKa["Unique_pos"] = pd.Series(MMS_maxes_conc_mKa["Cell_Barcode"]).apply(Cell_Bar_TO_upos)
MMS_maxes_conc_mKO["Unique_pos"] = pd.Series(MMS_maxes_conc_mKO["Cell_Barcode"]).apply(Cell_Bar_TO_upos)

# MMS_maxes_conc_mKO["Unique_pos"] = MMS_maxes_conc_mKO.iloc[0]["Cell_Barcode"][0:
# MMS_maxes_conc_mKO.iloc[0]["Cell_Barcode"].find("c")]

# %%
MMS_maxes_conc_mKO

# %%
MMS_maxes_conc_mKO

# %%
# mKOSubset_x99median = mKOSubset[["x99thPercentile_norm_BKGRND_Median_GFP", "Frame"]]

# pos_swarm = sns.swarmplot(data = MMS_maxes_conc_mKa,
# x="Relocalized",
# y= "x95thPercentile_norm_BKGRND_Median_GFP",
# hue="Unique_pos").set_title("DDC1")
# plt.savefig("DDC1_0218r2p900200")


# %%



pos_violin = sns.violinplot(data = df,
x="Relocalized",
y= "x95thPercentile_norm_BKGRND_Median_GFP",
hue="Unique_pos").set_title(MMS_maxes_conc_mKO["Protein"].values[0])
plt.savefig(f"x95thbyObjMedian_mKO")


# %%

pos_swarm = sns.violinplot(data = MMS_maxes_conc_mKa,
x="Relocalized",
y= "x95thPercentile_norm_BKGRND_Median_GFP",
hue="Unique_pos").set_title(MMS_maxes_conc_mKa["Protein"].values[0])
plt.savefig(f"x95thbyObjMedian_mKa")

# %%
# sns.boxplot(x= Quant_ALL["x10thPercentile_norm_OBJ_Median_GFP"])
sns.boxplot(x= Quant_ALL["x99thPercentile_norm_OBJ_Median_GFP"], color= 'red')

# %%
Quantification_shift

# %%
info_simple = pd.read_csv(info_simple.csv)
info_simpe.set_index(["Unique_frame", "Unique_cell"])

Quant_ALL = Quant_ALL[["Cell_Barcode", "Unique_frame", "Protein", "Relocalized"]]
# Quant_ALL =
Quant_ALL.set_index(["Unique_frame", "Cell_Barcode"], inplace = True)

subset_feature_move = pd.merge(Quant_ALL, info_simple, left_on=[('Unique_frame', 'Cell_Barcode')], right_on=[('Unique_frame', 'Unique_cell')])
subset_feature_move.to_csv("subset_feature_move.csv")

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import MinMaxScaler

subset_feature_move_copy = subset_feature_move
subset_feature_move_copy.set_index("Cell_Barcode", drop = True, inplace= True)
scalar = preprocessing.StandardScaler().fit(subset_feature_move_copy)
subset_feature_move_copy.reset_index(inplace=True)
le = LabelEncoder()
subset_feature_move_copy["cN"] =le.fit_transform(subset_feature_move_copy['Cell_Barcode'].astype(str))
subset_feature_move_copy.set_index("Cell_Barcode", drop = True, inplace= True)

# %%
import seaborn as sns
subset_feature_move_copy = subset_feature_move_copy.dropna()

kmeans = KMeans(init="random",
                n_clusters=6,
                n_init=10,
                max_iter=500,
                random_state=42)

kmeans.fit(subset_feature_move_copy)

kmeansDF = subset_feature_move_copy.copy()
kmeansDF['labels'] = kmeans.labels_
kmeansDF = kmeansDF.sort_values(by='labels', ascending=1)
del kmeansDF['labels']


colors = ["#f0ece0", "#474f67"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

# kmeansDF = kmeansDF.rename(index=dict(zip(list(kinomeMap['StrainID']),list(kinomeMap['Gene']))))

sns.heatmap(kmeansDF, cmap=cmap, vmax=0.5)

# %%
def myplot(score,coeff,labels=None):
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    plt.scatter(xs * scalex,ys * scaley,s=5)
    for i in range(n):
        plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
        if labels is None:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), color = 'green', ha = 'center', va = 'center')
        else:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')

    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))
    plt.grid()

myplot(components[:,0:2],np.transpose(pca.components_[0:2, :]),list(df.columns))
plt.show()

plt.savefig("pca_vector_contrib")

# %%
subset_shift = Quant_shift[Quant_shift["Cell_Barcode"] == "d0221r1p500200c0001"]
sns.boxplot(x=Quantification_shift["x20thPercentile_norm_OBJ_Median_GFP"], color= 'blue')

# %%
Quant_ALL_copy

# %%
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import MinMaxScaler

Quant_ALLcopy = Quant_ALL.drop(columns = ["Date", "ImageID", "Myo1Identity", "Protein", "TrackID_valid", "mKO_direction", "mKa_direction"])
Quant_ALLcopy.reset_index(inplace=True, drop= True)
le = LabelEncoder()
Quant_ALLcopy["cN"] = le.fit_transform(Quant_ALLcopy['Cell_Barcode'].astype(str))
# Quant_ALLcopy.set_index(["Frame", "Cell_Barcode"], drop = True, inplace= True)
Quant_ALLcopy.set_index("Cell_Barcode", drop = True, inplace= True)
scalar = preprocessing.StandardScaler().fit(Quant_ALLcopy)

# %%
Quant_ALLcopy

# %%
import seaborn as sns

test = Quant_ALLcopy
test_final = test.dropna()

kmeans = KMeans(init="random",
                n_clusters=6,
                n_init=10,
                max_iter=500,
                random_state=42)

kmeans.fit(test_final)

kmeansDF = test_final.copy()
kmeansDF['labels'] = kmeans.labels_
kmeansDF = kmeansDF.sort_values(by='labels', ascending=1)
del kmeansDF['labels']


colors = ["#f0ece0", "#474f67"]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)

# kmeansDF = kmeansDF.rename(index=dict(zip(list(kinomeMap['StrainID']),list(kinomeMap['Gene']))))

sns.heatmap(kmeansDF, cmap=cmap, vmax=0.5)


# %%
Quant_ALLcopy

# %%
#PCA
import plotly.express as px
from sklearn.decomposition import PCA

df = Quant_ALLcopy
n_components = 3

pca = PCA(n_components=n_components)
components = pca.fit_transform(df)

total_var = pca.explained_variance_ratio_.sum() * 100

labels = {str(i): f"PC {i+1}" for i in range(n_components)}
labels['color'] = 'Frame'

fig = px.scatter_matrix(
    components,
    color= df["Frame"],
    dimensions=range(n_components),
    labels=labels,
    title=f'Total Explained Variance: {total_var:.2f}%',
)
fig.update_traces(diagonal_visible=False)
fig.show()

# %%
def myplot(score,coeff,labels=None):
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())
    plt.scatter(xs * scalex,ys * scaley,s=5)
    for i in range(n):
        plt.arrow(0, 0, coeff[i,0], coeff[i,1],color = 'r',alpha = 0.5)
        if labels is None:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, "Var"+str(i+1), color = 'green', ha = 'center', va = 'center')
        else:
            plt.text(coeff[i,0]* 1.15, coeff[i,1] * 1.15, labels[i], color = 'g', ha = 'center', va = 'center')

    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))
    plt.grid()

myplot(components[:,0:2],np.transpose(pca.components_[0:2, :]),list(df.columns))
plt.show()

# plt.savefig("pca_vector_contrib.pdf")

# %%
df.reset_index(inplace=True)

# %%
import plotly.express as px
fig = px.scatter_3d(
    components, x=0, y=1, z=2, color=df['Frame'],
    title=f'Total Explained Variance: {total_var:.2f}%',
    labels={'0': 'PC 1', '1': 'PC 2', '2': 'PC 3', 'color': 'Frame'},
    hover_name=df["Cell_Barcode"]
)
fig.show()
# fig.write_html("PCA_3d.html")

# %%
from sklearn.preprocessing import LabelEncoder

Quant_ALL_copy = Quant_ALL.drop(columns = ["Date", "ImageID", "Myo1", "Protein"])
le = LabelEncoder()
Quant_ALL_copy["cN"] =le.fit_transform(Quant_ALL_copy['Cell_Barcode'].astype(str))
Quant_ALL_copy.set_index("Cell_Barcode", drop = True, inplace = True)

# %%
Quant_ALL_copy

# %%
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import plotly.express as px
# Read data from a csv
# z_data = Quant_ALL_copy
# z = z_data.drop(columns= ["Frame", "cN"])
# sh_0, sh_1 = z.shape
z = Quant_ALL_copy["x99thPercentile_norm_OBJ_Median_GFP"]
x = Quant_ALL_copy["Frame"]
y = Quant_ALL_copy["cN"]

fig = px.scatter_3d(Quant_ALL_copy, x='Frame', y='cN', z=z,
              color='cN')

# fig = go.Figure(data=go.Surface(z=z, x=x, y=y))
# fig.update_layout(title='Test', autosize=True)
fig.show()

# %%
import plotly.express as px
from joblib import Parallel, delayed

os.chdir(microfluidics_results)
if not os.path.exists("Figurs"):
    os.mkdir("Figures")

pos = "testing"

def tdim_metric_auto(metric):
    x = Quant_ALL_copy["Frame"]
    y = Quant_ALL_copy["cN"]
    z = Quant_ALL_copy[metric]

    fig = px.scatter_3d(Quant_ALL_copy, x='Frame', y='cN', z=z,
              color='cN')

    fig.write_html(f"Figures/{pos}_{metric}.html")
    return(1)

cols = Quant_ALL_copy.drop(columns = ["Frame", "cN"]).columns

pr = os.cpu_count() - 1
Parallel(n_jobs=pr)(delayed(tdim_metric_auto)(c) for c in cols)


# def 3dMetric_Frame_cN(metric,):
# 	metric
# 	fig.write_image(f"{}")


# z = Quant_ALL_copy["x90thPercentile_norm_OBJ_Median_GFP"]
# x = Quant_ALL_copy["Frame"]
# y = Quant_ALL_copy["cN"]

# fig = px.scatter_3d(Quant_ALL_copy, x='Frame', y='cN', z=z,
#               color='cN')

# # fig = go.Figure(data=go.Surface(z=z, x=x, y=y))
# # fig.update_layout(title='Test', autosize=True)

# %%


# %%
import plotly.graph_objects as go

fig = go.Figure(go.Surface(
    contours = {
        "x": {"show": True, "start": 1.5, "end": 2, "size": 0.04, "color":"white"},
        "z": {"show": True, "start": 0.5, "end": 0.8, "size": 0.05}
    },
    x = [1,2,3,4,5],
    y = [1,2,3,4,5],
    z = [
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 1],
        [0, 1, 0, 1, 0]
    ]))
fig.update_layout(
        scene = {
            "xaxis": {"nticks": 20},
            "zaxis": {"nticks": 4},
            'camera_eye': {"x": 0, "y": -1, "z": 0.5},
            "aspectratio": {"x": 1, "y": 1, "z": 0.2}
        })
fig.show()


