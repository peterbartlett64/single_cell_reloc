#%%
import pandas as pd
import os
import single_cell_reloc_parquet.global_functions.global_variables as gv
from datetime import datetime
from glob import glob
import seaborn as sns
import plotly.express as px
import single_cell_reloc_parquet.Post_quant.abundance_gen as cr_abund
import math
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import statistics

sns.set(style = 'whitegrid', palette = 'pastel', color_codes = True)
#%%
root = gv.slash_switch(input('Path'))
os.chdir(root)
#%%
year = str(datetime.today().year)
df_DDC2_all = pd.DataFrame([])
for path in glob('*.parquet', recursive = True):
	df = pd.read_parquet(path).reset_index(drop = False)
	col_name = path[path.find("Col"):path.find(year)-1]
	df = df.loc[df['Protein'].str.startswith("DDC2")] #* Filter protein names staring with "DDC2"
	df_DDC2_all = pd.concat([df_DDC2_all, df],ignore_index=True)
	print(path)
df_DDC2_all.to_parquet('DDC_all.parquet')
#%%
df_DDC2_all = pd.read_parquet('DDC_all.parquet')
df_DDC2_all['Run'] = pd.Series(df_DDC2_all['ImageID']).apply(lambda x: x[:x.find("p")])
df_DDC2_all = cr_abund.Abundance_calc_manager(df_DDC2_all)
#%%
barcs_keep = df_DDC2_all.groupby(['Run'])['Cell_Barcode'].sample(frac = 0.9)
df_DDC2_keep = pd.DataFrame(df_DDC2_all.loc[df_DDC2_all["Cell_Barcode"].isin((barcs_keep))].groupby(['Run','Cell_Barcode'])[['Loc_score', 'Abundance']].median()).reset_index(drop =False)
#%%
compare_Loc = px.box(df_DDC2_keep, x = 'Run', y = 'Loc_score', hover_data=['Cell_Barcode'])
compare_Loc.update_layout(showlegend=False)
compare_Loc.write_image(f"DDC2_compare_Loc.pdf")
compare_Loc.write_html(f"DDC2_compare_Loc.html")
compare_Loc
#%%
compare_Abundance = px.box(df_DDC2_keep, x = 'Run', y = 'Abundance', color = 'Run')
compare_Abundance.update_layout(showlegend = False)
compare_Abundance.write_image(f"DDC2_compare_abundance.pdf")
compare_Abundance.write_html(f"DDC2_compare_abundance.html")
compare_Abundance
#%%
fig, ax = plt.subplots(1, figsize = (12, 6))
sns.violinplot(data = df_DDC2_keep, x = 'Run', y = 'Abundance', hue= 'Run')
sns.despine()

plt.xticks(rotation = 35, ha = 'right')
plt.yticks(np.arange(4.5, 5.25, 0.25))
plt.show()

#%%
test = sns.stripplot(data = df_DDC2_keep, x = 'Run', y = 'Loc_score', hue
='Run', legend = False)

# compare_Abundance = sns.stripplot(data = df_DDC2_keep, x = 'Run', y = 'Abundance', hue ='Run', legend = False)
#%%




# %%

#%%
#, This needs to be added to the Individual Pearson plots
fig = px.scatter(df, x="", y="", color="Cell_Barcode", facet_col="Frame")
#%%
pen_frame = [25.76271186440678, 24.641833810888254, 20.17010935601458, 21.125143513203216, 27.86885245901639, 41.38297872340426, 19.187675070028014, 28.767942583732058, 19.906542056074766, 20.789473684210527, 27.106227106227106, 18.457101658255226, 22.90836653386454, 19.510135135135133, 23.227383863080682, 27.020648967551626, 22.85237698081735, 19.577464788732392, 30.278422273781903, 26.650755767700872, 21.960297766749378, 24.456521739130434, 22.49422632794457, 51.515151515151516, 30.14705882352941, 55.27950310559007, 52.33751425313569]

mean = statistics.mean(pen_frame)
for p in pen_frame:
	print(p - mean)