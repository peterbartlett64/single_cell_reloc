#%%
import pandas as pd
from glob import glob
import os
import plotly.express as px

#%% Do execution
if __name__ == "__main__":
	#, Define the varaibles
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
		'microfluidics_results': 'E:/Microfluidics/RESULTS',
		'post_path': 'E:/Microfluidics/MOST FINAL', #Todo: This needs to be changed to a input call
		'subset': False,
		'subset_by': '',
		'subset_collection': '',
		'cpu_se': os.cpu_count(),
		'timepoint_gap': 7.5,
		'percentiles': [95, 99],
		'multiplex': True}

	path_ten = "E:/Microfluidics/MOST FINAL"
	path_three = "E:/Microfluidics/TRY THREE"
	os.path.join(path_ten,"Final_combined_comparison.parquet")
	file_ten = pd.read_parquet(os.path.join(path_ten,"Final_combined_comparison.parquet")).loc[:, ["Percentage_reloc", "Percentage_reloc_less", "Protein"]].set_index("Protein", drop = True)
	file_three = pd.read_parquet(os.path.join(path_three,"Final_combined_comparison.parquet")).loc[:, ["Percentage_reloc", "Percentage_reloc_less", "Protein"]]

	combine_compare_time = pd.merge(file_ten, file_three, left_on = 'Protein', right_on='Protein', how = 'outer', suffixes=("_10", "_3"))
	combine_compare_time['reloc_Diff'] = combine_compare_time["Percentage_reloc_10"] - combine_compare_time["Percentage_reloc_3"]
	combine_compare_time['less_Diff'] = combine_compare_time["Percentage_reloc_less_10"] - combine_compare_time["Percentage_reloc_less_3"]

#%%
if __name__ == '__main__':
	compare_time_less = px.bar(combine_compare_time, x = 'Protein', y = 'reloc_Diff')
	compare_time_less.write_html('compare_times_Reloc_less.html')
	compare_time_less.write_image('compare_times_Reloc_less.pdf')

	compare_time_cross = px.bar(combine_compare_time, x = 'Protein', y = 'reloc_Diff')
	compare_time_cross.write_html('compare_times_Reloc_cross.html')
	compare_time_cross.write_image('compare_times_Reloc_cross.pdf')
