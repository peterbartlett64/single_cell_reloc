#%%
import pandas as pd
import os
import single_cell_reloc_parquet.global_functions.global_variables as glv
from icecream import ic
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

#%%
#. Dictionaries were removed because no longer combining reloc paterns across lit. datasets
#! SHould this information be required again, it is included in the original script

#, I don't know
class dataset_desc:
	def __init__(self, merge_on, information) -> None:
		self.merge_on = merge_on
		self.information = information
		self.func_match = f"{self}_f"

# #This is callable funciton to get information on the microscope
# class microcope_info():
# 	def __init__(self) -> None:
# 		self.pixel_ratio_microns = 0.1081

# class GO_term:
# 	def __init__(self, go_list) -> None:
# 		self.listed = str.split(go_list, ", ")

#. Leave in case of future use
def f_GO_map(micro_map): #! This was left out because decided to display based on a network rather than colored scatter
	go = pd.read_excel(os.path.join(Global_Variables['information_path'], "GO_Proteins.xlsx"))
	# go = pd.read_excel("C:\\Users\\pcnba\\Grant Brown's Lab Dropbox\\Peter Bartlett\\Peter Bartlett Data\\Code\\Data_copies\\Information_files\\Localization_merging\\GO_Proteins.xlsx")
	terms = go['TERM']

	# def add_GOs(prot:str, go_group_prot:str, go_group_name:str, go_group_list:list): #! This is not in use because deicided to use a network graph to visualize interaction enrichment rather than as scatter with variable grouped GO color
	#. This below code has not been micro_maped
	# 	if prot in go_group_list:
	# 		go_group_prot + "," + go_group_name
	# 	else:
	# 		pass
	# 	return(go_group_prot)

	# micro_map['GO_group_collected'] = ''
	# for r in terms:
	# 	go_matches = go.loc[r, 'ANNOTATED_GENES']
	# 	micro_map['GO_group_collected'] = micro_map.apply(lambda x: add_GOs(x['Protein'], x['GO_group_collected'], r, go_matches), axis = 1)
	return(go) #. , micro_map)

def sgd_map_f():
	# sgd = pd.read_csv("results_best.csv")
	sgd = pd.read_csv("results.tsv", sep = '\t')
	# sgd.rename(columns={'input':'Gene_Standard_Name'}, inplace=True)
	# sgd.rename(columns={'length':'Gene_Length'}, inplace=True)
	sgd.columns = sgd.columns.str.replace(" > ", "_")
	sgd.columns = sgd.columns.str.replace(" ", "_")
	# sgd = sgd.loc[:, ["Gene_Systematic_Name", "Gene_Standard_Name", "Gene_Name", "Gene_Length", "Gene_Phenotype_Summary", "Gene_Length"]]
	sgd['Gene_Standard_Name'] = sgd['Gene_Standard_Name'].fillna(sgd['Gene_Systematic_Name'])
	sgd = sgd.set_index("Gene_Systematic_Name")
	#* global info_sgd
	#* info_sgd = dataset_desc('Standard_Name','information')
	return(sgd)

def tkach_f():
	tkach = pd.read_excel("Tcak_protein_localization.xlsx", sheet_name='Calls')
	tkach.columns = tkach.columns.str.replace(" ", "_")
	tkach['Standard_Name'] = tkach['Standard_Name'].fillna(tkach['Systematic_ORF'])
	tkach = tkach.set_index("Systematic_ORF")

	tkach['EndLOC_Rescreen_MMS_Tcak'] = tkach['EndLOC_Rescreen_MMS_Tcak'].map(Tkach_dictionary).fillna(tkach['EndLOC_Rescreen_MMS_Tcak'])
	tkach['EndLOC_Rescreen_HU_Tcak'] = tkach['EndLOC_Rescreen_HU_Tcak'].map(Tkach_dictionary).fillna(tkach['EndLOC_Rescreen_HU_Tcak'])
	#* global info_tkach
	#* info_tkach = dataset_desc('Standard_Name','localization')
	return(tkach)

def microfluidics_map_f():
	microfluidics_map = pd.read_excel("MicrofluidicsMap_wCol_USE.xlsx", sheet_name="ProteinLocations")
	microfluidics_map.dropna(subset='Protein', inplace=True)
	microfluidics_map.columns = microfluidics_map.columns.str.replace(' ', '_')
	microfluidics_map['Protein'] = microfluidics_map['Protein'].str.upper()
	microfluidics_map = microfluidics_map.set_index('Protein')
	#* global info_microfluidcs
	#* info_microfluidcs = dataset_desc('Standard_Name/Mix','map')
	return(microfluidics_map)

def denervaud_ycd_f():
	denervaud_ycd = pd.read_excel("Den_data_bestgood.xlsx", sheet_name='Sheet1')
	# denervaud_ycd.rename(columns={'denervaud_ycd_Call':'Call'}, inplace=True)
	denervaud_ycd.drop(denervaud_ycd.columns[denervaud_ycd.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
	# denervaud_ycd = denervaud_ycd.add_suffix("_denervaud_ycd")
	# denervaud_ycd['Standard_Name_denervaud_ycd'] = denervaud_ycd['Standard_Name_denervaud_ycd'].str.upper()
	denervaud_ycd.fillna('unclassified', inplace = True)

	denervaud_ycd['initial_localization'] = denervaud_ycd['initial_localization'].map(Den_ycd_map_dict).fillna(denervaud_ycd['initial_localization'])
	denervaud_ycd['end_localization'] = denervaud_ycd['end_localization'].map(Den_ycd_map_dict).fillna(denervaud_ycd['end_localization'])
	denervaud_ycd.sort_values(by = "movieTag", inplace = True) #* Put in order so that the best movie is first before taking the first instance of an ORF label
	denervaud_ycd['geneName'] = denervaud_ycd['geneName'].replace({'-': np.nan})
	denervaud_ycd['geneName'] = denervaud_ycd['geneName'].fillna(denervaud_ycd['yORF'])
	denervaud_ycd = denervaud_ycd.groupby('yORF').aggregate(lambda x: x.iloc[0])
	denervaud_ycd = denervaud_ycd.drop(columns=['geneName', 'exp_cond', 'movieTag'])
	return(denervaud_ycd)

# def denervaud_f():
# 	denervaud = pd.read_excel("Denervaud Calls.xlsx", sheet_name='Sheet1')
# 	denervaud.columns = denervaud.columns.str.replace(' ', '_')
# 	denervaud.rename(columns={'Denervaud_Call':'Call'}, inplace=True)
# 	denervaud.drop(denervaud.columns[denervaud.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
# 	denervaud = denervaud.add_suffix("_Denervaud")
# 	denervaud['Standard_Name_Denervaud'] = denervaud['Standard_Name_Denervaud'].str.upper()

# 	denervaud['Call_Denervaud'] = denervaud['Call_Denervaud'].map(Denervaud_dictionary).fillna(denervaud['Call_Denervaud'])
# 	denervaud = denervaud.set_index('Standard_Name_Denervaud')
# 	#* global info_Denervaud
# 	#* info_Denervaud = dataset_desc('Standard_Name','Localization')
# 	return(denervaud)

def Ho_loc_pen_f():
	Ho_loc = pd.read_excel("Loc_Ho_SuppT5.xlsx", sheet_name="20210809_HUMMS_penetrance")
	Ho_loc.drop(Ho_loc.columns[Ho_loc.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
	Ho_loc.drop(Ho_loc.columns[Ho_loc.columns.str.contains('HU',case = False)],axis = 1, inplace = True)
	Ho_loc.columns = Ho_loc.columns.str.replace(" ", "_")
	Ho_loc = Ho_loc.add_suffix("_Ho")
	Ho_loc = Ho_loc.set_index('Gene_Ho')
	#* global info_Ho
	#* info_Ho = dataset_desc('Standard Name','LocPen')
	return(Ho_loc)

def Mazumder_f():
	mazumder = pd.read_excel("Mazumder_ver2.xlsx", sheet_name="Mod_dest")
	# mazumder.columns = mazumder.columns.str.replace(" ", "_")
	mazumder = mazumder.add_suffix("_Mazumder")
	mazumder['Localization_Mazumder'] = mazumder['Localization_Mazumder'].map(Mazumder_dictionary).fillna(mazumder['Localization_Mazumder'])
	mazumder['CommName_Mazumder'] = mazumder['CommName_Mazumder'].fillna(mazumder['ORF_Mazumder'])
	mazumder = mazumder.set_index('ORF_Mazumder')
	#* global info_Mazumder
	#* info_Mazumder = dataset_desc('ORF','Localization')
	return(mazumder)

def Huh_f():
	Huh = pd.read_csv('Huh2003.txt', sep = '\t').set_index('yORF')
	Huh['localization summary'] = Huh['localization summary'].map(Huh_dict).fillna(Huh['localization summary'])
	return(Huh)


#%%
# def origin_destination(directional:str):
# 	orig_dest = directional.str.split(by = ' -> ')
# 	return(orig_dest) #! must be unpacked

# def merge_manager():
# 	os.chdir(Global_Variables['information_path'])
# 	files = input("What are the files to be micro_map based on protein? [Comma deliminate]").split(', ')
# 	for f in files: #TODO: Add multi-read_functionality
# 		pd.read_csv(f).set_index()
# 	micro_map = pd.DataFrame([])
# 	file_n = 0
# 	for f in files:
# 		exec(f'Global {f}')
# 		exec(f"temp = {f}_f({f})") #! This isn't good form, but will work for now
# 		if file_n > 0:
# 			micro_map = pd.merge(micro_map, temp, left_index = True, right_index = True)
# 			file_n += 1
# 		else:
# 			micro_map = temp.copy()
# 			file_n += 1

#%%
def control_replace(x,search):
	ic(x, search)
	if search.upper() in x.upper():
		return(search.upper())
	else:
		return(x.upper())
#%%
#, Desigante the main function call
if __name__ == "__main__":
	# Global_Variables = glv.global_manager()
	Global_Variables = {
	# 	"analyze": "F:/Microfluidics/Missing_Analyze2",
		"microfluidics_results": "F:/Microfluidics/RES_N_ULTS",
		"information_path": "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging",
		"post_path": "D:/ALL_FINAL"} # * ,
	# 	"subset": False,
	# 	'subset_by': 'range',
	# 	'subset_collection': '',
	# 	"cpu_se": 16,
	# 	"timepoint_gap": 7.5,
	# 	"percentiles": [95, 99],
	# 	"multiplex": True}
	os.chdir(Global_Variables['information_path'])
	# files = input("What are the files to be micro_map based on protein? [Comma deliminate]").split(', ')
	# for f in files: #TODO: Add multi-read_functionality
		# pd.read_csv(f).set_index()

	# micro_map = microfluidics_map_f()
	# # .drop(columns=['Predicted_localization_Change', 'Notes', 'Current_Stage', 'Location', 'Fullmicro_map'])
	# map_drop_columns = ['Predicted_localization_Change', 'Notes', 'Current_Stage', 'Location', 'Fullmicro_map']
	# micro_map.drop(micro_map.columns[micro_map.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
	# micro_map.drop(columns=[col for col in micro_map if col in map_drop_columns], inplace=True)

	sgd = sgd_map_f()
	# micro_map = pd.merge(sgd, micro_map, left_index=True, right_index=True, how = 'left')
	#. Decided that will just use the list of proteins from SGD to merge with the other datasets

	#< tkach = tkach_f()
	#< # denervaud = denervaud_f()
	#< mazumder = Mazumder_f()
	#< Den_ycd_map_df = denervaud_ycd_f()
	#< Brandons_map = pd.read_excel("Brandons_Paper.xlsx", sheet_name="Sheet1").set_index('ORF')
	#< Brandons_map = Brandons_map.drop(columns=['Protein']).rename(columns={'Subcellular Compartment Re-localization': 'Dest_Call'}).add_suffix('_Brandons')
	#< Huh = Huh_f()

	#, Read in the dataframes to be joined
	tkach = pd.read_excel("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging/Tkach_refined.xls", sheet_name='Localization scoring')
	tkach = tkach.rename(columns={'MMS localization change class': 'MMS_localization_class', 'HU Localication change class': 'HU_localization_class'})
	micro_map = sgd.merge(tkach, left_on = "Gene_Standard_Name", right_on = "Standard Name", how = 'left')
	micro_map['MMS_HU_merged_class'] = micro_map['MMS_localization_class'].fillna(micro_map['HU_localization_class']) #* This did not end of getting used



	#< loqate = pd.read_excel('proteomesummarylamicro_mapversion.xlsx', sheet_name='Sheet1', usecols=['ORF', 'control Localization']).set_index('ORF').replace('below threshold', np.nan)

	#< micro_map = sgd.merge(Den_ycd_map_df, left_index= True, right_index= True, how = 'left')
	#< micro_map = micro_map.merge(tkach, left_index=True, right_index = True, how= "left")

	#< #Either should work below
	#< micro_map = micro_map.merge(mazumder, left_index = True, right_index = True, how = "left")
	#< # micro_map = micro_map.merge(mazumder, left_on = 'Gene_Standard_Name', right_on = 'CommName_Mazumder', how = "left")

	#< micro_map = pd.merge(micro_map, Brandons_map, right_index=True, left_index=True, how= 'left')
	#< #Artifact of removed micorfluidics map
	#< # micro_map = micro_map.sort_values(by = ['Date', 'Run_Number', 'MapID_(Col_Range)'])

	#< micro_map = micro_map.merge(Huh, left_index=True, right_index=True, how='left')
	#< micro_map = micro_map.merge(loqate, left_index=True, right_index=True, how='left')

#%%
#. These are all the protiens which do not have matching localization information
#Check the rows where Na for 'localization_change', 'initial_localization', 'end_localization', 'EndLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_HU_Tcak', Localization_Mazumder', 'Dest_Mazumder', 'Function_Mazumder'

# >micro_map.loc[micro_map[['localization_change', 'initial_localization', 'end_localization','EndLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_HU_Tcak', 'Localization_Mazumder', 'Dest_Mazumder', 'Dest_Call_Brandons']].isna().all(axis = 1),['Gene_Standard_Name', 'localization_change', 'initial_localization', 'end_localization','EndLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_HU_Tcak', 'Localization_Mazumder', 'Dest_Mazumder']]

#%%
#< micro_map[['StartLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_MMS_Tcak']] = micro_map['EndLOC_Rescreen_MMS_Tcak'].str.split(' -> ', expand=True).replace('', np.nan, regex=True)
#< micro_map[['StartLOC_Rescreen_HU_Tcak', 'EndLOC_Rescreen_HU_Tcak']] = micro_map['EndLOC_Rescreen_HU_Tcak'].str.split(' -> ', expand=True).replace('', np.nan, regex=True)

#%%
#< micro_map['EndLOC_Rescreen_HU_Tcak'] = micro_map['EndLOC_Rescreen_HU_Tcak'].str.replace(' -> ', '').replace('', np.nan, regex=True)
#< micro_map['EndLOC_Rescreen_HU_Tcak'] = micro_map['EndLOC_Rescreen_HU_Tcak'].str.replace('-> ', '').replace('', np.nan, regex=True)

# micro_map['EndLOC_Rescreen_MMS_Tcak'] = micro_map['EndLOC_Rescreen_MMS_Tcak'].str.replace(' -> ', '').replace('', np.nan, regex=True)
# micro_map['EndLOC_Rescreen_MMS_Tcak'] = micro_map['EndLOC_Rescreen_MMS_Tcak'].str.replace('-> ', '').replace('', np.nan, regex=True)

#%%
#< def simp_pref(*args):
#< 	for arg in args:
#< 		if pd.isna(arg):
#< 			continue
#< 		elif arg == 'unclassified':
#< 			continue
#< 		elif arg == 'nothing':
#< 			continue
#< 		else:
#< 			return arg
#< 	return np.nan

#< micro_map['Combined_destination'] = micro_map.apply(lambda x: simp_pref(x['Dest_Call_Brandons'], x['end_localization'], x['EndLOC_Rescreen_MMS_Tcak'], x['EndLOC_Rescreen_HU_Tcak'], x['Dest_Mazumder']), axis=1)

#< micro_map['Combined_origin'] = micro_map.apply(lambda x: simp_pref(x['StartLOC_Rescreen_MMS_Tcak'], x['StartLOC_Rescreen_HU_Tcak'], x['localization summary']), axis=1)

#< micro_map['Combined_destination'] = micro_map['Combined_destination'].fillna('unclassified').str.replace(',', ' and ').str.title().str.strip()
#< micro_map['Combined_origin'] = micro_map['Combined_origin'].fillna('unclassified').str.replace(',', ' and ').str.title().str.strip()

#%%
#* Split the combined origins and destinations into their respective levels of granularity
#< micro_map[['Single_origin', 'Double_origin', 'Triple_origin', 'Quad_origin']] = micro_map['Combined_origin'].str.split(' And ', expand=True)
#< micro_map[['Single_destination', 'Double_destination', 'Triple_destination', 'Quad_destination']] = micro_map['Combined_destination'].str.split(' And ', expand=True)

# micro_map = micro_map.fillna(value=np.nan)
#%%
#* Define a dictionary to rename the origin and destination to a more general term
#< ori_dest_dict = {
# 	'Nucleolus Irregular': 'Nucleolus',
# 	'Cytoplasmic Foci': 'Cytoplasm',
# 	'Cytoplasm Irreg.': 'Cytoplasm',
# 	'Punctate': 'Cytoplasm foci',
# 	'Nuclear Periphery': 'Nucleus',
# 	'Cyto Foci': 'Cytoplasm foci',
# 	'Cell Periphery': 'Cytoplasm',
# 	'Plasma Membrane': 'Cytoplasm',
# 	'Nucleus P': 'Nucleus',
# 	'Er Foci': 'ER',
#     'Microtubules': 'Cytoplasm',
# 	'Vacuole Foci' : 'Vacuole',
# 	'Bud': 'Bud',
# 	'Er': 'ER',
# 	'Nothing': 'Unclassified',
#     'Valuole Foci': 'Vacuole',
# 	'Mitochondrion': 'Mitochondria',
# 	'Pm (Punctate)': 'Plasma Membrane',
# 	'Nuclear / Cyto Foci': 'Nuclear and Cyto Foci',
# 	'Bud Neck': 'Bud Neck'}

#< #* Rename each origin level
#< micro_map['Single_origin'] = micro_map['Single_origin'].map(ori_dest_dict).fillna(micro_map['Single_origin'])
#< micro_map['Double_origin'] = micro_map['Double_origin'].map(ori_dest_dict).fillna(micro_map['Double_origin'])
#< micro_map['Triple_origin'] = micro_map['Triple_origin'].map(ori_dest_dict).fillna(micro_map['Triple_origin'])
#< micro_map['Quad_origin'] = micro_map['Quad_origin'].map(ori_dest_dict).fillna(micro_map['Quad_origin'])
#%%
#* Recombined the levels to have layered granularity instead of indivual terms
#< micro_map['Double_origin'] = micro_map['Single_origin'] + " and " + micro_map['Double_origin']
#< micro_map['Double_origin'] = micro_map['Double_origin'].fillna(micro_map['Single_origin'])
#%%
#< micro_map['Triple_origin'] = micro_map['Double_origin'] + " and " + micro_map['Triple_origin']
#< micro_map['Triple_origin'] = micro_map['Triple_origin'].fillna(micro_map['Double_origin'])
#%%
#< micro_map['Quad_origin'] = micro_map['Triple_origin'] + " and " + micro_map['Quad_origin']
#< micro_map['Quad_origin'] = micro_map['Quad_origin'].fillna(micro_map['Triple_origin'])

#%%
#* Rename each destination level
#< micro_map['Single_destination'] = micro_map['Single_destination'].map(ori_dest_dict).fillna(micro_map['Single_destination'])
#< micro_map['Double_destination'] = micro_map['Double_destination'].map(ori_dest_dict).fillna(micro_map['Double_destination'])
#< micro_map['Triple_destination'] = micro_map['Triple_destination'].map(ori_dest_dict).fillna(micro_map['Triple_destination'])
#< micro_map['Quad_destination'] = micro_map['Quad_destination'].map(ori_dest_dict).fillna(micro_map['Quad_destination'])

#* Recombined the destination levels to have layered granularity instead of indivual terms
#< micro_map['Double_destination'] = micro_map['Single_destination'] + " and " + micro_map['Double_destination']
#< micro_map['Double_destination'] = micro_map['Double_destination'].fillna(micro_map['Single_destination'])

#< micro_map['Triple_destination'] = micro_map['Double_destination'] + " and " + micro_map['Triple_destination']
#< micro_map['Triple_destination'] = micro_map['Triple_destination'].fillna(micro_map['Double_destination'])

#< micro_map['Quad_destination'] = micro_map['Triple_destination'] + " and " + micro_map['Quad_destination']
#< micro_map['Quad_destination'] = micro_map['Quad_destination'].fillna(micro_map['Triple_destination'])

#%%%
# micro_map['Single_origin'] = micro_map['Single_origin'].map(ori_dest_dict).fillna(micro_map['Single_origin'])
# micro_map['Double_origin'] = micro_map['Double_origin'].map(ori_dest_dict).fillna(micro_map['Double_origin'])
# micro_map['Triple_origin'] = micro_map['Triple_origin'].map(ori_dest_dict).fillna(micro_map['Triple_origin'])
# micro_map['Quad_origin'] = micro_map['Quad_origin'].map(ori_dest_dict).fillna(micro_map['Quad_origin'])
# # reference a pandas column and split into required columns by delimiter
# micro_map['Single_destination'] = micro_map['Single_destination'].map(ori_dest_dict).fillna(micro_map['Single_destination'])
# micro_map['Double_destination'] = micro_map['Double_destination'].map(ori_dest_dict).fillna(micro_map['Double_destination'])
# micro_map['Triple_destination'] = micro_map['Triple_destination'].map(ori_dest_dict).fillna(micro_map['Triple_destination'])
# micro_map['Quad_destination'] = micro_map['Quad_destination'].map(ori_dest_dict).fillna(micro_map['Quad_destination'])


# micro_map[['Primary']] = micro_map['Combined_origin'].str.split(' And ', expand=True).replace('', np.nan, regex=True)

#%%
micro_map.to_parquet('protein_origin_dest_new.parquet')

#%%
# Read in the percentages and merge with the micro_map
penetrances = pd.read_parquet("D:\ALL_FINAL\Final_combined_comparison.parquet")
#%%
penetrances = penetrances.merge(micro_map, left_on = 'Protein', right_on='Gene_Standard_Name', how='left')
#%%
yet_percentiles = pd.read_parquet("D:/ALL_FINAL/Combined_by_perc/new_percs.parquet")

#ADD the trimming to make sure that proteins are less than 32 frames (4 hours)
yet_percentiles = yet_percentiles.loc[yet_percentiles['Frames_post_treatment'] <= 32]



#%%
updated_yet_perc = yet_percentiles.groupby('Protein').Yet_perc.agg('max').rename('updated_yet_perc')

penetrances = penetrances.merge(updated_yet_perc, left_on = 'Protein', right_on='Protein', how='left') #* This has extra information that is not needed

#%%
Ho = Ho_loc_pen_f()
Ho_agg = Ho.agg('max', axis= 1)*100 #*Convert the decimal to a percentage
Ho_agg = Ho_agg.rename('Ho_max').to_frame()
merged_frame_pens = penetrances.merge(Ho_agg, left_on = 'Protein', right_index = True, how = 'left')
merged_frame_pens.to_parquet('D:/ALL_FINAL/Combined_by_perc/penetrance_updated_trimmed.parquet')

#%%
#* Reading in this file is a time limiting step
Loc_data = pd.read_parquet("D:\ALL_FINAL\Combined_by_perc\merged_data_final.parquet", columns=["Cell_Barcode", "Date", "Frame", "Unique_Frame", "Protein", "Is_treated", "Frames_post_treatment", "Unique_pos", "Loc_score", "Relocalized", "Abundance", "log_Abundance", "log_Loc_score", "z_score_Loc", "z_score_Abund", "z_score_logLoc", "z_score_logAbund", "ImageID", "Progen_bud", "Yet", "Does", "When", "Yes_yet", "No_yet", "pres_end", "c_count_end_nl", "c_count_end_p5"])
#%%
Loc_data_merged = micro_map.merge(Loc_data, left_on = 'Gene_Standard_Name', right_on = 'Protein', how = 'right')
#%%
Loc_data["LogAbundance"] = pd.Series(Loc_data["Abundance"]).apply(lambda x: np.log2(x))

import scipy.stats as stats
grouped_pearson = Loc_data.groupby(["Protein", "Frame"]).apply(lambda x:pd.Series(stats.pearsonr(x.z_score_Loc, x.z_score_logAbund), index=["corr", "pval"]))
medians_pearson = grouped_pearson.groupby(level=0).agg('median')
medians_pearson.to_csv('median_pearson.csv')

grouped_spearman = Loc_data.groupby(["Protein", "Frame"]).apply(lambda x:pd.Series(stats.spearmanr(x.z_score_Loc, x.z_score_logAbund), index=["corr", "pval"]))
medians_spearmans = grouped_spearman.groupby(level=0).agg('median')
medians_spearmans.to_csv('median_spearman.csv')

grouped_kendall = Loc_data.groupby(["Protein", "Frame"]).apply(lambda x:pd.Series(stats.kendalltau(x.z_score_Loc, x.z_score_logAbund), index=["corr", "pval"]))
medians_kendlall = grouped_kendall.groupby(level=0).agg('median')
medians_kendlall.to_csv('median_kendall.csv')



def gen_pearson_r_w_p(Loc_data):





	ProtFrameCorr_df = grouped_correlations.unstack().iloc[:,1].rename('ProtFrameCorrs').to_frame().reset_index(drop = False)
	med_corr =  ProtFrameCorr_df.copy().reset_index(drop = False).groupby("Protein")['ProtFrameCorrs'].agg('median').rename('MedianProtCorr').to_frame().reset_index(drop = False)

	med_corr.sort_values(by = ['MedianProtCorr'], inplace = True)
	return(ProtFrameCorr_df, med_corr)

Loc_data = Loc_data
ProtFrameCorr_df, med_corr = gen_pearson_r(Loc_data)
Loc_data = Loc_data.merge(ProtFrameCorr_df, right_on=['Protein','Frame'], left_on=['Protein', 'Frame'], how = 'left')
Loc_data = Loc_data.merge(med_corr, right_on="Protein", left_on= "Protein", how = 'left')

#%%
Loc_data_merged["LogAbundance"] = pd.Series(Loc_data_merged["Abundance"]).apply(lambda x: np.log2(x))
def gen_pearson_r(Loc_data):
	grouped_correlations = Loc_data.groupby(["Protein", "Frame"])[['Loc_score', 'LogAbundance']].corr(method = 'spearman')
	ProtFrameCorr_df = grouped_correlations.unstack().iloc[:,1].rename('ProtFrameCorrs').to_frame().reset_index(drop = False)
	med_corr =  ProtFrameCorr_df.copy().reset_index(drop = False).groupby("Protein")['ProtFrameCorrs'].agg('median').rename('MedianProtCorr').to_frame().reset_index(drop = False)
	med_corr.sort_values(by = ['MedianProtCorr'], inplace = True)
	return(ProtFrameCorr_df, med_corr)

Loc_data = Loc_data_merged
ProtFrameCorr_df, med_corr = gen_pearson_r(Loc_data)
Loc_data = Loc_data.merge(ProtFrameCorr_df, right_on=['Protein','Frame'], left_on=['Protein', 'Frame'], how = 'left')
Loc_data = Loc_data.merge(med_corr, right_on="Protein", left_on= "Protein", how = 'left')
#%%

# cols = ['Gene_Standard_Name', 'Combined_destination', 'Combined_origin',
#        'Single_origin', 'Double_origin', 'Triple_origin', 'Quad_origin',
#        'Single_destination', 'Double_destination', 'Triple_destination',
#        'Quad_destination', 'Cell_Barcode', 'Date', 'Unique_Frame', 'Protein', 'Unique_pos', 'ImageID', 'When']
# Loc_data[cols] = Loc_data[cols].astype('category')
table = pa.Table.from_pandas(Loc_data)
pq.write_table(table, 'D:\ALL_FINAL\Combined_by_perc\Loc_data_comp_merged_Tkach.parquet')

#%%
#Requested addions
Loc_data = pd.read_parquet('D:\ALL_FINAL\Combined_by_perc\Loc_data_comp_merged_Tkach.parquet')
#%%
Reloc_count = Loc_data.groupby(['Protein', 'Frames_post_treatment', 'Relocalized']).Cell_Barcode.agg('count')
Reloc_count = Reloc_count.rename('Reloc_state_count').unstack().rename(columns={0: 'CurrNot', 1: 'CurrYes'}).reset_index(drop = False)


Reloc_count.sort_values(by = ['Protein', 'Frames_post_treatment'], inplace = True)
Reloc_count['currProportion'] = Reloc_count['CurrYes']/(Reloc_count['CurrYes'] + Reloc_count['CurrNot'])
Reloc_count['RelocVelocity'] = Reloc_count.groupby(['Protein'])['currProportion'].diff() #* This is is just to see if there is a variable rate of relocalization because if cells are being relocalized at a constant rate or not at all
Reloc_count['RelocAcceleration'] = Reloc_count.groupby(['Protein'])['RelocVelocity'].diff() #* This is is just to see if there is a time where cells are relocalizing faster

#%%
Yet_count = Loc_data.groupby(['Protein', 'Frames_post_treatment', 'Yet']).Cell_Barcode.agg('count')
Yet_count = Yet_count.rename('Yet_state_count').unstack().rename(columns={0: 'CurrNotYet', 1: 'CurrYet'}).reset_index(drop = False)
Yet_count['currYetProportion'] = Yet_count['CurrYet']/(Yet_count['CurrNotYet'] + Yet_count['CurrYet'])

Yet_count['YetVelocity'] = Yet_count.groupby(['Protein'])['currYetProportion'].diff() #* This is is just to see if there is a variable rate of relocalization because if cells are being relocalized at a constant rate or not at all
Yet_count['YetAcceleration'] = Yet_count.groupby(['Protein'])['YetVelocity'].diff() #* This is is just to see if there is a time where cells are relocalizing faster


#%%
#, Calculate the proportion of cells that have relocalized
Does_count = Loc_data.groupby(['Protein', 'Frames_post_treatment', 'Does']).Cell_Barcode.agg('count')
Does_count = Does_count.rename('Does_state_count').unstack().rename(columns={0: 'countDoesNot', 1: 'countDoes'}).reset_index(drop = False)
Does_count.sort_values(by = ['Protein', 'Frames_post_treatment'], inplace = True)
Does_count['DoesProportion'] = Does_count['countDoes']/(Does_count['countDoesNot'] + Does_count['countDoes'])
Does_count['DoesVelocity'] = Does_count.groupby(['Protein'])['DoesProportion'].diff() #* This is is just to see if there is a variable rate of cell loss because if cells are being lost at a constant rate or not at all
#%%
def get_range(s):
 return s.max() - s.min()

Does_count['Does_FinDiff'] = Does_count.groupby(['Protein'])['countDoes'].transform(get_range) #* This is to see if there is a net loss of cells
#%%
Loc_data = Loc_data.merge(Reloc_count, left_on = ['Protein', 'Frames_post_treatment'], right_on = ['Protein', 'Frames_post_treatment'], how = 'left')
Loc_data = Loc_data.merge(Yet_count, left_on = ['Protein', 'Frames_post_treatment'], right_on = ['Protein', 'Frames_post_treatment'], how = 'left')
Loc_data = Loc_data.merge(Does_count, left_on = ['Protein', 'Frames_post_treatment'], right_on = ['Protein', 'Frames_post_treatment'], how = 'left')
#%%

table = pa.Table.from_pandas(Loc_data)
pq.write_table(table, 'D:\ALL_FINAL\Combined_by_perc\Loc_data_comp_merged_everything.parquet')

#%%




# Loc_data_merged = Loc_data_merged.loc[:,["Gene_Standard_Name", "Combined_destination", "Combined_origin", "Single_origin", "Double_origin", "Triple_origin", "Quad_origin", "Single_destination", "Double_destination", "Triple_destination", "Quad_destination", "Cell_Barcode", "Date", "Frame", "Unique_Frame", "Protein", "Is_treated", "Frames_post_treatment", "Unique_pos", "Loc_score", "Relocalized", "Abundance", "log_Abundance", "log_Loc_score", "z_score_Loc", "z_score_Abund", "z_score_logLoc", "z_score_logAbund", "ImageID", "Progen_bud", "Yet", "Does", "When", "Yes_yet", "No_yet", "pres_end", "c_count_end_nl", "c_count_end_p5", "LogAbundance"]]
# table = pa.Table.from_pandas(Loc_data_merged)
# pq.write_table(table, 'D:\ALL_FINAL\Combined_by_perc\Loc_data_comp_merged.parquet')
#%%
# Real_mols = pd.read_excel("Abundance_TableS8_MOD.xlsx", sheet_name="Table S8")
# Real_foldChange = pd.read_excel("foldAbundanceTableS9_MOD.xlsx", sheet_name="Table S9")


#%%
#> This is an older version that was changed to the above
# small = micro_map[['EndLOC_Rescreen_MMS_Tcak','Call_Denervaud','Localization_Mazumder']]
# def combine_logic(tkach, denervaud, mazumder):
# 	if len(tkach) > 0:
# 		t = True

# 	else: t = False

# 	if len(denervaud)>0:
# 		d = True
# 	else: d = False

# 	if len(mazumder) > 0:
# 		m = True
# 	else: m = False

# 	if t == True and d == True and m == True:
# 		if tkach == denervaud == mazumder:
# 			combo = tkach
# 		else:
# 			combo = tkach + denervaud + mazumder

# 	elif t == True:
# 		if d == True and m == False:
# 			if tkach == denervaud:
# 				combo = tkach
# 			else:
# 				combo = tkach + denervau
# 		elif d == False and m == True:
# 			if tkach == mazumder:
# 				combo = tkach
# 			else:
# 				combo = tkach + denervaud
# 		else:
# 			combo = tkach

# 	elif d == True and t == False:
# 		if m == True:
# 			if denervaud == mazumder:
# 				combo = mazumder
# 			else:
# 				combo = denervaud + mazumder
# 		else:
# 			combo = denervaud

# 	elif m == True and t == False and d == False:
# 		if mazumder == denervaud:
# 			combo = mazumder
# 		else:
# 			combo = denervaud + mazumder

# 	if t == False and d == False and m == False:
# 		return(None)

# 	return(combo)


# small['combined'] = small.apply(lambda x: combine_logic(tkach=x['EndLOC_Rescreen_MMS_Tcak'], denervaud= x['Call_Denervaud'], mazumder= x['Localization_Mazumder']), axis = 0)

# def join_strings(s):
	# return ''.join(str(i) for i in s if not pd.isna(i))

# %%
# #> This is to match the DDC2 to FLR1 in library
# #* The search and join for control should be done last once all the merging has been performed
# search = input("What is the name of the control protein? [Esc to skip]")
# if search == '':
# 	pass
# else:
# 	micro_map.reset_index(inplace=True, drop= False)
# 	micro_map['Protein'] = pd.Series(micro_map['Protein']).apply(lambda x: control_replace(x = x, search = search))
# 	micro_map.set_index('Protien')
