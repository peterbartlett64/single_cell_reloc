#%%
import pandas as pd
import os
import single_cell_reloc_parquet.global_functions.global_variables as glv
#%%
#TODO Finish the automated version for upload

class dataset_desc:
	def __init__(self, merge_on, information) -> None:
		self.merge_on = merge_on
		self.information = information
		self.func_match = f"{self}_f"

class microcope_info():
	def __init__(self) -> None:
		self.pixel_ratio_microns = 0.1081

def sgd_map_f():
	sgd = pd.read_csv("results_best.csv")
	sgd.columns = sgd.columns.str.replace(" > ", "_")
	sgd.columns = sgd.columns.str.replace(" ", "_")
	sgd = sgd.loc[:, ["Gene_Systematic_Name", "Gene_Standard_Name", "Gene_Name", "Gene_Length", "Gene_Phenotype_Summary", "Gene_Length"]]
	sgd['Gene_Standard_Name'] = sgd['Gene_Standard_Name'].fillna(sgd['Gene_Systematic_Name'])
	sgd = sgd.set_index("Gene_Standard_Name")
	#* global info_sgd
	#* info_sgd = dataset_desc('Standard_Name','information')
	return(sgd)

def tkach_f():
	tkach = pd.read_excel("Tcak_protein_localization.xlsx", sheet_name='Calls')
	tkach.columns = tkach.columns.str.replace(" ", "_")
	tkach = tkach.set_index("Standard_Name")
	#* global info_tkach
	#* info_tkach = dataset_desc('Standard_Name','localization')
	return(tkach)

def microfluidics_map_f():
	microfluidics_map = pd.read_excel("MicrofluidicsMap_wCol (version 1).xlsx", sheet_name="ProteinLocations")
	microfluidics_map.columns = microfluidics_map.columns.str.replace(' ', '_')
	microfluidics_map['Protein'] = microfluidics_map['Protein'].str.upper()
	microfluidics_map = microfluidics_map.set_index('Protein')
	#* global info_microfluidcs
	#* info_microfluidcs = dataset_desc('Standard_Name/Mix','map')
	return(microfluidics_map)

def denervaud_f():
	denervaud = pd.read_excel("Denervaud Calls.xlsx", sheet_name='Sheet1')
	denervaud.columns = denervaud.columns.str.replace(' ', '_')
	denervaud = denervaud.add_suffix("_Denervaud")
	denervaud['Standard_Name_Denervaud'] = denervaud['Standard_Name_Denervaud'].str.upper()
	denervaud = denervaud.set_index('Standard_Name_Denervaud')
	#* global info_Denervaud
	#* info_Denervaud = dataset_desc('Standard_Name','Localization')
	return(denervaud)

def Ho_loc_f():
	Ho_loc = pd.read_excel("Loc_Ho_SuppT5.xlsx", sheet_name="20210809_HUMMS_penetrance")
	Ho_loc.columns = Ho_loc.columns.str.replace(" ", "_")
	Ho_loc = Ho_loc.add_suffix("_Ho")
	Ho_loc = Ho_loc.set_index('Gene_Ho')
	#* global info_Ho
	#* info_Ho = dataset_desc('Standard Name','LocPen')
	return(Ho_loc)

def Mazumder_f():
	mazumder = pd.read_excel("Mazumder_localization.xlsx", sheet_name="All GFP intensities")
	mazumder.columns = mazumder.columns.str.replace(" ", "_")
	mazumder = mazumder.add_suffix("_Mazumder")
	mazumder = mazumder.set_index('ORF_Mazumder')
	#* global info_Mazumder
	#* info_Mazumder = dataset_desc('ORF','Localization')
	return(mazumder)

# def merge_manager():
# 	os.chdir(Global_Variables['post_path'])
# 	files = input("What are the files to be merged based on protein? [Comma deliminate]").split(', ')
# 	for f in files: #TODO: Add multi-read_functionality
# 		pd.read_csv(f).set_index()
# 	merged = pd.DataFrame([])
# 	file_n = 0
# 	for f in files:
# 		exec(f'Global {f}')
# 		exec(f"temp = {f}_f({f})") #! This isn't good form, but will work for now
# 		if file_n > 0:
# 			merged = pd.merge(merged, temp, left_index = True, right_index = True)
# 			file_n += 1
# 		else:
# 			merged = temp.copy()
# 			file_n += 1


def control_replace(x,search):
	if search.upper() in x.upper():
		return(search.upper())
	else:
		return(x.upper())

if __name__ == "__main__":
    # Global_Variables = glv.global_manager()
	Global_Variables = {
	# 	"analyze": "F:/Microfluidics/Missing_Analyze2",
		"microfluidics_results": "F:/Microfluidics/RES_N_ULTS",
		"post_path": "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files"} # * ,
	# 	"subset": False,
	# 	'subset_by': 'range',
	# 	'subset_collection': '',
	# 	"cpu_se": 16,
	# 	"timepoint_gap": 7.5,
	# 	"percentiles": [95, 99],
	# 	"multiplex": True}
	os.chdir(Global_Variables['post_path'])
	# files = input("What are the files to be merged based on protein? [Comma deliminate]").split(', ')
	# for f in files: #TODO: Add multi-read_functionality
		# pd.read_csv(f).set_index()

	merged = microfluidics_map_f()
	tkach = tkach_f()
	sgd = sgd_map_f()
	denervaud = denervaud_f()
	mazumder = Mazumder_f()

	search = input("What is the name of the control protein? [Esc to skip]")
	if search == '':
		pass
	else:
		merged['Protein'] = pd.Series(merged['Protein']).apply(lambda x: control_replace(x = x, search = search))
	merged = merged.merge(tkach, left_index=True, right_index = True, how= "outer")
	merged = merged.merge(denervaud, left_index = True, right_index = True, how = "outer")
	merged = merged.merge(sgd, left_index = True, right_index = True, how = "outer")
	merged = merged.merge(mazumder, left_on = "Gene_Systematic_Name", right_index = True, how = "outer")

	merged = merged.sort_values(by = ['Date', 'Run_Number', 'MapID_(Col_Range)'])
	print("DONE merging")
	merged.dropna(subset='Date', inplace=True)

