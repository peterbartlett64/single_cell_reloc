#%%
import pandas as pd
import os
import parquet

#%%
if __name__ == '__main__':
	# ddc2_confirmation = pd.read_csv('C:/Users/pcnba/OneDrive - University of Toronto/Desktop/Ddc2_confirmation.csv')
	# test = ddc2_confirmation.groupby(['Protein', 'Frames_post_treatment', 'ImageID']).cell_count.agg('count')
	lib_confirmation = pd.read_parquet("D:\ALL_FINAL\Combined_by_perc\Quant_ALL.parquet", columns= ['Protein', 'Frames_post_treatment', 'ImageID', 'cell_count'])
	# lib_confirmation = pd.read_csv('C:/Users/pcnba/OneDrive - University of Toronto/Desktop/lib_confirmation.csv')
	#%%
	test = lib_confirmation.groupby(['Protein', 'Frames_post_treatment', 'ImageID']).cell_count.agg('count')
	test = test.reset_index(drop = False)
	final = test.loc[(test['Frames_post_treatment'] < 0) & (test['Frames_post_treatment'] >= -10) & (test['cell_count'] < 50)].copy()
	final['Pos'] = pd.Series(final['ImageID']).apply(lambda x: x[:-5])
	final2 = final.groupby(['Pos']).cell_count.agg('sum').reset_index(drop = False)
	#%%
	final2 = final2.loc[final2['cell_count'] < 100]['Pos']
	final3 = final.loc[final['Pos'].isin(final2)]['Pos'].unique()
	final4 = final.loc[final['Pos'].isin(final3)]['Protein'].unique()
#%%
	final = test.loc[(test['Frames_post_treatment'] <= 32) & (test['Frames_post_treatment'] >= -10) & (test['cell_count'] < 50)].copy()
	final['Pos'] = pd.Series(final['ImageID']).apply(lambda x: x[:-5])
	final5 = final['Protein'].unique()

# %%
