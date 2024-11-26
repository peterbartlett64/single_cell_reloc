#%%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
#%%
regular =pd.read_parquet("D:/ALL_FINAL/Combined_by_perc/Quant_ALL.parquet", columns=['x95thPercentile_norm_BKGRND_Median_GFP', 'x99thPercentile_norm_BKGRND_Median_GFP', 'Frames_post_treatment', 'Cell_Barcode', 'Myo1Identity','x95thPercentile_GFP_RAW', 'x99thPercentile_GFP_RAW', 'x95thPercentile_norm_OBJ_Median_GFP', 'x99thPercentile_norm_OBJ_Median_GFP', 'Unique_Frame'])
regular.drop_duplicates(inplace=True)
#%%

agg_median = regular.groupby(['Frames_post_treatment', 'Unique_Frame', 'Myo1Identity']).aggregate('median').reset_index(drop = False)
agg_median = agg_median.loc[agg_median['Frames_post_treatment'] < 32]
#%%
os.chdir("D:/Second Mind/Academic/Project Stuff/Figures")
#%%

sns.lineplot(data = agg_median,
			#  hue='Myo1Identity',
			x='Frames_post_treatment',
			y='x95thPercentile_norm_OBJ_Median_GFP',
			errorbar=('ci', 95),
			alpha = 0.75)
#%%
plt.savefig("Frames_pre_treatment_95_Obj.svg", format='svg')
#%%
sns.lineplot(data = agg_median,
			x='Frames_post_treatment',
			y='x99thPercentile_norm_OBJ_Median_GFP',
			errorbar=('ci', 95),
			alpha = 0.75)
plt.savefig("Frames_pre_treatment_99_Obj.svg", format='svg')