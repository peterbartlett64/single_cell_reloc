#%%
#, Import packages
import pandas as pd
import os

#, Define functions
#*Calculate the proportion of cells that are not relocalized currently but have been relocalized in the past
def deloc_calc(row):
	if row["Yet"] == 1 and row["Relocalized"] == 1:
		deloc = 1
	# else: #* This is shorter but possibly less readable and reliable
	# 	return(0)
	elif row["Yet"] == 1 and row["Relocalized"] == 0:
		deloc = 0
	elif row["Yet"] == 0 and row["Relocalized"] == 1:
		deloc = 0
	elif row["Yet"] == 0 and row["Relocalized"] == 0:
		deloc = 0
	else:
		return(9)
	return(deloc)
#%%
#, Main function call
if __name__ == '__main__':
	#. Load temporary functions
	Global_Variables
	Global_Variables = {
	# 	"analyze": "F:/Microfluidics/Missing_Analyze2",
		"microfluidics_results": "F:/Microfluidics/RES_N_ULTS",
		"information_path": "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging",
		"post_path": "D:/ALL_FINAL"}

	df = pd.read_parquet("D:\ALL_FINAL\Combined_by_perc\Loc_data_comp_merged_complex.parquet", engine='pyarrow')
	#* Temporary subset of the data
	df = df.loc[(df['Protein'] == 'FLR1')]
#%%
#, Calcuate the number of cell that have delocalized
sorted_df = df.sort_values(by = ['Protein', 'Frames_post_treatment'], ascending=True)
sorted_df['Delocalized'] = sorted_df.apply(deloc_calc, axis = 1)
#%%

# sorted_df['Delocalized_count'] = sorted_df.groupby(['Protein', 'Frames_post_treatment']).transform('count')
sorted_df['State_count'] = sorted_df['Delocalized'] #* Define a temporary column that will be counted later
deloc_info = sorted_df.groupby(['Protein', 'Frames_post_treatment', 'Delocalized'])['State_count'].agg('count')

#%%
deloc_info = deloc_info.rename(columns = {'Delocalized': 'Delocalized_count'})
deloc_info.set_index(['Protein', 'Frames_post_treatment'], inplace = True)


#%%
df = df.groupby(['Protein', 'Frames_post_treatment', 'D']).agg({'Does': 'max', 'Yet': 'max'}).reset_index(drop = False)

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