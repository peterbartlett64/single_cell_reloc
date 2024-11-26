#, f_prepare for STRING network
#%%
import pandas as pd
import os
import arrow
import math

#%%
if __name__ == "__main__":
	# Global_Variables = gv.global_manager()
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': 'D:/ALL_FINAL', #. gv.slash_switch(input("Post quant path?")) , #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': int(math.floor(os.cpu_count()*0.7)),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True,
	'figures_root': 'D:/Figures_root'}

	post_combined_path = os.path.join(Global_variables['post_path'], 'Combined_by_perc')
	os.chdir(post_combined_path)

df = pd.read_parquet("D:\\ALL_FINAL\\Final_combined_comparison.parquet", columns=['Protein', 'Percentage_reloc_less'])
# %%
df.astype({'Percentage_reloc_less': 'int16'})
#%%
df.to_csv('for_network.csv', index= False)
