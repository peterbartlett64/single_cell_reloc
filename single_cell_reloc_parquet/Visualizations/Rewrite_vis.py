
#%%
from datetime import date
import os
import pandas as pd
import math
import os
import statistics
from joblib import Parallel, delayed
from scipy import stats
import plotly.express as px
import single_cell_reloc_parquet.global_functions.global_variables as gv
import seaborn as sns
import matplotlib.pyplot as plt
import plotly
plotly.io.kaleido.scope.mathjax = None #* This is an important line for timely svg output. This is to stop a call to the internet for prettier looking math
import kaleido
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_context('paper')

print(plotly.__version__, kaleido.__version__)

# sns.set(rc = {'figure.figsize':(15,8)}) #* Set the figure size to be larger
# sns.set(rc = {'figure.figsize':(40,30)}) #* Set the figure size to be the shape of a powerpoint slide

#%%
def load_location_info(info_path = "C:\\Users\\pcnba\\Grant Brown's Lab Dropbox\\Peter Bartlett\\Peter Bartlett Data\\Code\\Data_copies\\Tcak_protein_localization.xls"):
	location_information = pd.read_excel(info_path, sheet_name='Calls')
	location_information.dropna(subset = ['Rescreen_MMS'], inplace = True) #* Drop proteins that do not have localizatoin change in MMS

	#, Define a rename dictionary to simplify the localization changes
	rename_dict = {
		'to cyto' : 'Cytoplasm',
		'to nuc' : 'Nucleus',
		#Todo: Add in figure legend
		'nucleus, nuc foci': 'Nucleus',#! MULTI for MCM2
		'cyto foci' : 'Cytoplasm-foci',
		'from budneck/tip': 'Bud neck',
		'from bud tip': 'Bud neck',
		'to pm' : "Plasma Membrane",
		'to bud neck': 'Bud neck',
		'to nuc/periphery': 'Nucleus-peripheral',
		'to nuc (from nuc foci)': 'Nucleus',
		'nuc foci': 'Nucleus-foci',
		'from bud neck': 'Bud neck',
		'to vac' : 'Vacuole',
		'to nuc, nuc foci, cyto foci': 'Nucleus', #!MULTI for YDR132C
		'no cells': 'DROP',
		'to pm foci': 'Plasma Membrabe',
		'to nucleolus, cyto foci' : 'Nucleolus', #* Right now only a single location can be dealt with
		'to vac (from pm)': 'Vacuole',
		'to nucleolus': 'Nucleolus',
		'to nuc periph (cyto?)':'Nucleus-peripheral',
		'nuc foci (weak)': 'Nucleus-foci',
		'to cyto (from pm/endosome)': 'Cytoplasm',
		'to nuc (from nucleolus)': 'Nucleus',
		'er foci': 'ER-foci',
		'to vac (abund inc)': 'Vacuole',
		'to vac (focus)': 'Vacuole-focus',
		'to cyto (from pm)': 'Cytoplasm',
		'shorter microtubules': 'Microtubules',
		'to nuc/nuc periph': 'Nucleus-peripheral',
		'to cyto (daughter)': 'Cytoplasm',
		'to diffuse nuc (really, more numerous less intense foci, from foci)': 'Nucleus-diffuse',
		'to cyto (decreased abund)': 'Cytoplasm',
		'from budneck/tip, to pm' : "Plasma Membrane", #* This is an example of simplying to just desitination
		'to pm foci (endosome)': 'Plamsma Membrane',
		'to nuc (nucleolus)': 'Nucleolus',
		'to nucleolus (weak)': 'Nucleolus',
		'not the same strain': 'DROP',
		'to er' : 'ER',
		'from budneck': 'Bud neck',
		'nuc periph' : 'Nucleus-peripheral',
		'to pm (foci)' : 'Plasma Membrane Foci',
		'to vac (from vac mb)': 'Vacuole',
		'from bud neck (to pm)' : 'Plasma Membrane',
		'to diffuse nuc (from foci)': 'Nucleus-diffuse',
		'to er foci (weak)': 'ER-foci',
		'nuc periph foci':'Nucleus-foci'
	}

	#* Do the renaming
	def rename_local(x, rename_dict):
		modified_location = rename_dict[x]
		return(modified_location)

	#Todo: Finish making the multi location version
	multi_locs = ['nucleus, nuc foci', 'to nuc, nuc foci, cyto foci']
	multi_prots = location_information.loc[location_information['Rescreen_MMS'].isin(multi_locs)]["Standard Name"]

	location_information['Location_simplfied'] = pd.Series(location_information['Rescreen_MMS']).apply(rename_local, rename_dict = rename_dict)
	localization_simple = location_information.copy().set_index('Standard Name')[['Systematic ORF', 'Location_simplfied']]
	return(localization_simple)
#%%
def gen_pearson_r(df):
	post_treatment_df = df.loc[df["Frames_post_treatment"] >= 0] #* Subset to frames after treatment

	Agg_loc = post_treatment_df.groupby(["Protein", "Frame"])["Loc_score"].agg('median') #* First get the median of Loc_scores in each protein frame

	Agg_loc = pd.DataFrame(Agg_loc)
	Agg_loc.reset_index(drop = False, inplace = True) #* Reassociate the Protein and Frame to median Loc_score

	max_agg_Loc_score = Agg_loc.loc[Agg_loc.groupby(["Protein"])["Loc_score"].idxmax()][["Protein", "Frame"]] #*Create a table of maximum aggregated Loc_scores for each protein

	# for p in max_agg_Loc_score["Protein"]:
	# 	pearson_r = df.apply(lambda x: stats.pearsonr(x["z_score_Abund"], x["z_score_Loc"]))

	return(max_agg_Loc_score)

#%%
def graph_pearson_single(i, full_data, max_agg_Loc_score):
	prot = max_agg_Loc_score.iloc[i,:]["Protein"]
	f = max_agg_Loc_score.iloc[i,:]["Frame"]
	try:
		sub_data = full_data.loc[(full_data["Protein"] == prot) & (full_data["Frame"] == f)].dropna(subset=['z_score_Loc', 'z_score_Abund'])
	except ValueError:
		sub_data["z_score_Abund"] = stats.zscore(sub_data.loc[:,"Abundance"])
		sub_data["z_score_Loc"] = stats.zscore(sub_data.loc[:,"Loc_score"])
	z_Loc_Abun_r = stats.pearsonr(sub_data["z_score_Loc"], sub_data["z_score_Abund"])
	scatter = sns.scatterplot(data= sub_data, x = "z_score_Abund", y = "z_score_Loc")
	fig = scatter.get_figure()
	fig.savefig(f"Pearson{z_Loc_Abun_r}--{prot}_{f}.pdf")
	print(f"{prot} is done pearson'n")
	return(prot, z_Loc_Abun_r)

def graph_pearson_plus_minus(i, full_data, max_agg_Loc_score):
	try:
		prot = max_agg_Loc_score.iloc[i,:]["Protein"]
		f = max_agg_Loc_score.iloc[i,:]["Frame"]
		try:
			sub_data = full_data.loc[(full_data["Protein"] == prot) & ((full_data["Frame"] == f)| (full_data["Frame"] == f+1)| (full_data["Frame"] == f-1))].dropna(subset=['z_score_Loc', 'z_score_Abund'])
		except ValueError:
			# sub_data["z_score_Abund"] = stats.zscore(sub_data.loc[:,"Abundance"])
			# sub_data["z_score_Loc"] = stats.zscore(sub_data.loc[:,"Loc_score"])
			sub_data["z_score_Abund"] = sub_data.groupby("Frame")["Abundance"].transform(stats.zscore)
			sub_data["z_score_Loc"] = sub_data.groupby("Frame")["Loc_score"].transform(stats.zscore)
		z_Loc_Abun_r = stats.pearsonr(sub_data["z_score_Loc"], sub_data["z_score_Abund"])
		scatter = sns.scatterplot(data= sub_data, x = "z_score_Abund", y = "z_score_Loc", hue = "Frame")
		fig = scatter.get_figure()
		fig.savefig(f"Pearson{z_Loc_Abun_r}--{prot}_{f}.pdf")
		print(f"{prot} is done pearson'n")
		return(prot, z_Loc_Abun_r)
	except Exception as e:
		return(f"{e} on {prot}")

#%%# %%
#! Correct this to actually run
if __name__ == "__main__":
	# Global_Variables = gv.global_manager()
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': 'D:/TRY AGAIN', #. gv.slash_switch(input("Post quant path?")) , #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': int(math.floor(os.cpu_count()*0.7)),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True,
	'figures_root': 'D:/Figures_root'}

	#, Make directory with figures, but if it exists skip
	figures_root = Global_variables['figures_root']
	date_today = str(date.today())
	figures = os.path.join(figures_root, date_today) #* This was recenlty added to make sure that the figues are not overwriting
	try:
		os.mkdir(figures)
	except FileExistsError:
		pass

	time_gap = Global_variables['timepoint_gap']
	post_path = Global_variables["post_path"]
	post_path = os.path.join(post_path, 'Combined_by_perc')
	os.chdir(post_path)

	#* Set final file name
	output_file = 'Final_w_Abund.parquet'
	#, Load in file if present, otherwise create. If fail exit
	#* This is good when this script has already been run so it can save time, but also go auto if not.
	try:
		full_data = pd.read_parquet(output_file, columns = ['Cell_Barcode','Date', 'Frame', 'Unique_Frame', 'Protein', 'Is_treated', 'Frames_post_treatment', 'Unique_pos', 'Loc_score', 'Relocalized', 'CoV_spos','CoV_apos', 'Abundance', 'z_score_Loc', 'z_score_Abund']) #* This speeds it up a TON by not loding extra information
	except FileNotFoundError:
		# Create an empty DataFrame to store the merged data
		merged_data = pd.DataFrame([])
		# Iterate over each file in the directory
		c = 0
		for filename in os.listdir(post_path):
			c += 1
			if filename.endswith('.parquet'):
				# Load the Parquet file into a DataFrame
				file_path = os.path.join(post_path, filename)
				df = pd.read_parquet(file_path, columns = ['Cell_Barcode','Date', 'Frame', 'Unique_Frame', 'Protein', 'Is_treated', 'Frames_post_treatment', 'Unique_pos', 'Loc_score', 'Relocalized', '<lambda_0>_x','<lambda_0>_y', 'Abundance', 'z_score_Loc', 'z_score_Abund']).reset_index(drop = False)

				# Append the data to the merged_data DataFrame
				merged_data = pd.concat([merged_data, df], ignore_index=True)
				print(c)

		# Specify the path for the output merged Parquet file
		merged_data = merged_data.rename(columns={'<lambda_0>_x': "CoV_spos",'<lambda_0>_y': "CoV_apos"})
		output_file = os.path.join(post_path, output_file)

		# Write the merged data to a single Parquet file
		merged_data.to_parquet(output_file)

		print(f'Merged data saved to {output_file}')
		full_data = merged_data
	except:
		print('Could not find full_data file in post_quant dir. Exiting...')
		exit()

	localization_simple = load_location_info() #? This does take in a variable info path, but for right now just using the default Tkach verison
	full_data_wLocal = pd.merge(full_data, localization_simple, left_on= 'Protein', right_index=True, how = 'outer')
	full_data_wLocal = full_data_wLocal.dropna(subset='Date')
	full_data['Run'] = pd.Series(full_data['Unique_Frame']).apply(lambda x: x[:7])
	full_data = full_data.loc[~(full_data["Run"] == 'd0212r1')]

	os.chdir(post_path) #* It likely already is this directory but no harm doing it anyway
	full_data_wLocal.to_parquet('Final_wAbundLocal.parquet')

	input("pause")
	# Abundance_df = full_data_wLocal.reset_index(inplace = True, drop = True)

	max_agg_Loc_score, full_data_wLocal = gen_pearson_r(full_data_wLocal)
	#* This requires the pearson corr to alread by attached/performed

	os.chdir(figures)
	x = Parallel(n_jobs= Global_variables['cpu_se'], verbose= 100)(delayed(graph_pearson_single)(i = i, full_data = full_data,  max_agg_Loc_score = max_agg_Loc_score) for i in range(len(max_agg_Loc_score)))
	print(x)

	# x = Parallel(n_jobs= pr, verbose= 100)(delayed(graph_pearson_plus_minus)(i = i, full_data = smaller,  max_agg_Loc_score =  max_agg_Loc_score) for i in range(len(max_agg_Loc_score)))
	# print(x)

	agg_way =  ['median', 'max']
	os.chdir(figures)

	for ag in agg_way: #* It might be a better idea to zip execute with color
		agg_pearson = full_data_wLocal.groupby(['Protein'])['Pearson_coeff_Loc_Abund'].agg(ag)
		pearson_coeff_bar = px.bar(agg_pearson.sort_values())
		pearson_coeff_bar.write_html(f'Pearson_r_bar_mono_{ag}.html')
		pearson_coeff_bar.write_image(f'Pearson_r_bar_mono_{ag}.pdf')

	# for p in full_data['Proteins']:
	Var_ap = full_data[["CoV_apos", "Protein", "Frames_post_treatment", "Run"]].copy().drop_duplicates().reset_index(drop = True)
	Var_ap = Var_ap.loc[Var_ap["Frames_post_treatment"] >= 0]
	Var_ap["Time_post_treatment"] = Var_ap["Frames_post_treatment"] * Global_variables['timepoint_gap']
	Var_ap['Median_CoV_apos'] = Var_ap.groupby(['Protein']).CoV_apos.transform('median') #* Generate median CV for each protien

	# Generate a median
	def f_median_day_graph():
		Var_ap['Median_Day_CoV_apos'] = Var_ap.groupby(['Run']).Median_CoV_apos.transform('median')
		Var_ap_day_check = Var_ap.loc[:, ['Run', "Median_Day_CoV_apos"]].drop_duplicates().sort_values(by = 'Run')
		Run_CoV_bar = px.bar(Var_ap_day_check, x = 'Run', y = 'Median_Day_CoV_apos')
		Run_CoV_bar.write_html('Median_Day_CoV_apos.html')
		Run_CoV_bar.write_image('Median_Day_CoV_apos.pdf')

	f_median_day_graph()
	def f_bin_CV_Median():
		Var_ap = Var_ap.sort_values(by = ["Median_CoV_apos", "Frames_post_treatment"])
		n_bin_plot = 5

		#* This a manual bin strategy
		Var_prot_list = Var_ap['Protein'].unique()
		Var_ap_prot_len = len(Var_prot_list)
		Var_ap_first = list(Var_prot_list[:5])
		Var_ap_last = list(Var_prot_list[-5:])
		middle = math.floor(Var_ap_prot_len/2)
		Var_ap_middle = list(Var_prot_list[middle - 2: middle + 4])
		Var_ap_sel =   Var_ap_first + Var_ap_middle + Var_ap_last

		Var_ap_sub_df = Var_ap.loc[Var_ap['Protein'].isin(Var_ap_sel)].copy()
		Var_line = px.scatter( Var_ap_sub_df, x = "Time_post_treatment", y = "CoV_apos", color = "Protein")
		#. This did not generate very intersting results
		return(Var_line)


	#, Figure to see if there is an assocaition
	pentrance_data = pd.read_parquet("D:/TRY AGAIN/Final_combined_comparison.parquet") #Todo: Change this from hard coding
	merged_CV_pen = pd.merge(pentrance_data, Var_ap, left_on = 'Protein', right_on = 'Protein')
	merged_CV_pen = merged_CV_pen[['Protein', 'Percentage_reloc_less', 'Percentage_reloc', 'Median_CoV_apos']].copy().drop_duplicates()
	Var_pen_compare_scatter = sns.scatterplot(data = merged_CV_pen, x = "Percentage_reloc_less", y= 'Median_CoV_apos', hue='Protein')
	# Var_commpare_scatter = px.scatter(test, x = 'Percentage_reloc_less', y = 'Median_CoV_apos', color = 'Protein') #* There does not seem to be an association. Should explore this further


	#, Write the Comparison between RRP5 and FLR1
	Var_ap_RRP5_FLR1 = Var_ap.loc[Var_ap['Protein'].isin(['RRP5', 'FLR1', 'FGV2'])]
	ax = sns.lineplot(data= Var_ap_RRP5_FLR1, x = "Time_post_treatment", y = "CoV_apos", hue = "Protein")
	ax.set(xlabel='Time post treatment (minutes)', ylabel='Coefficient of Variation (CV)')
	plt.savefig('RRP5-FLR1_CV.pdf', format = 'pdf')

	after_copy = full_data.loc[full_data['Frames_post_treatment']>= 0].copy()
	def reloc_yet(x):
		ind = x.idxmax()
		does = x.loc[ind]
		x.loc[:ind] = 0
		x.loc[ind:] = does
		return(x)
	def workaround(ind):
		return(after_copy.loc[ind, "Unique_Frame"])
	def does_workaround(ind):
		return(after_copy.loc[ind,"Relocalized"])

	after_copy["Yet"] = after_copy.groupby(by = "Cell_Barcode")["Relocalized"].transform(reloc_yet) # This repesents wether there has been relocaliztion yet
	after_copy["ind"] = after_copy.groupby(by = "Cell_Barcode")["Relocalized"].transform('idxmax')
	after_copy["Does"] = pd.Series(after_copy["ind"]).apply(does_workaround) #this will work for now but should make it come out of one of the other functions
	after_copy["When"] = pd.Series(after_copy["ind"]).apply(workaround)

	after_copy['Time_post_treatment'] = after_copy['Frames_post_treatment'] * Global_variables['timepoint_gap']

	subset_FLR1 = after_copy.loc[(after_copy['Protein'] == 'FLR1') & (after_copy['Frames_post_treatment'].isin([0,8,16,24,32,40]))]
	plot_FLR1 = sns.stripplot(x = "Time_post_treatment", y= "Loc_score", data = subset_FLR1, hue = 'Yet')
	plot_FLR1.axhline(y = 1)
	plot_FLR1.set(xlabel= 'Time post treatment (minutes)', ylabel = 'LocScore')
	fig = plot_FLR1.get_figure()
	fig.savefig("FLR1_strip.pdf")

	subset_RRP5 = after_copy.loc[(after_copy['Protein'] == 'RRP5') & (after_copy['Frames_post_treatment'].isin([0,8,16,24,32,40]))]
	plot_RRP5 = sns.violinplot(data = subset_RRP5, x = "Time_post_treatment", y= "Loc_score")#,hue = 'Yet', split = True)
	plot_RRP5.axhline(y = 1)
	plot_RRP5.set(xlabel= 'Time post treatment (minutes)', ylabel = 'LocScore')
	fig = plot_RRP5.get_figure()
	fig.savefig(f"RRP5_Violin_unsplit.pdf", format = 'pdf')
	fig.savefig(f"RRRP5_Violin_unsplit.svg", format = 'svg')

	for p in after_copy['Protein'].unique():
		subset = after_copy.loc[(after_copy['Protein'] == p) & (after_copy['Frames_post_treatment'].isin([0,8,16,24,32,40]))].copy()
		plot = sns.violinplot(data = subset, x = "Time_post_treatment", y= "Loc_score", inner = 'quartile') #,hue = 'Yet, split = True)
		plot.axhline(y = 1)
		plot.set(xlabel= 'Time post treatment (minutes)', ylabel = 'LocScore')
		fig = plot.get_figure()
		fig.savefig(f"{p}_Violin.pdf", format = 'pdf')
		del subset
		del plot
		del fig
#%%
barcs_keep = subset_RRP5.groupby(['Unique_pos'])['Cell_Barcode'].sample(n = 10)
yet_yes = pd.DataFrame(subset_RRP5.loc[(subset_RRP5['Yet'] == 1) & (subset_RRP5['Frames_post_treatment'] == 0)]['Cell_Barcode'])
yes_keep = yet_yes.sample(n =5)['Cell_Barcode']
barcs_keep = pd.concat([barcs_keep, yes_keep])

subset_RRP5_smaller = subset_RRP5.loc[subset_RRP5["Cell_Barcode"].isin((barcs_keep))]
plot_RRP5 = sns.lineplot(data = subset_RRP5_smaller, x = "Time_post_treatment", y= "Loc_score")
plot_RRP5.axhline(y = 1)
plot_RRP5.set(xlabel= 'Time post treatment (minutes)', ylabel = 'LocScore')
fig = plot_RRP5.get_figure()
#%%
different_subset = after_copy.loc[(after_copy['Protein'] == 'RRP5')] #*This is not subsetted by time because the result will be a line. There are no groups so this will generate much smoother lines
primed = different_subset.loc[(different_subset['Relocalized'] == 1) & (different_subset['Frames_post_treatment'] == 0)]['Cell_Barcode']
different_subset['Group'] = 0 #* temp set the column to 1
different_subset.loc[(different_subset['Cell_Barcode'].isin(primed)), 'Group'] = 1 #* Change all barcodes were it was positive for reloc at treatment time
#%%
line_RRP5 = sns.lineplot(data = different_subset, x = 'Time_post_treatment', y = 'Loc_score', hue = 'Group')
line_RRP5.axhline(y = 1)
line_RRP5.set(xlabel= 'Time post treatment (minutes)', ylabel = 'LocScore')
fig = line_RRP5.get_figure()
fig.savefig("lineplot_RRR5.pdf", format = 'pdf')
#%%
test = pd.DataFrame(primed.apply(lambda x: x[:x.find("c")]))
test['t'] = True #* This is just for pandas so it has a column to count elements
test.groupby('Cell_Barcode')['t'].count()