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
# plotly.io.kaleido.scope.mathjax = None #* This is an important line for timely svg output. This is to stop a call to the internet for prettier looking math
import kaleido
sns.set(rc={'figure.figsize':(11.7,8.27)})
sns.set_context('paper')

print(plotly.__version__, kaleido.__version__)
#%%
# df_cross = pd.read_csv()
# df_longitude = pd.read_csv()

# def f_conv_name(p):
# 	return(p[:p.find("-")])

# df_cross['Protien'] = pd.Series(df_cross['Protien']).apply(f_conv_name)
# df_longitude['Protien'] = pd.Series(df_longitude['Protien']).apply(f_conv_name)

# df_merged = pd.merge(df_cross, df_longitude, by = "Protein")


df_merged = pd.read_parquet("D:\ALL_FINAL\95th_percentile\All_pos_pt_percentage_sync.parquet")
#%%
df_merged = df_merged.loc[df_merged['Frames_post_treatment'] <42]
px.line(df_merged, x = 'Frames_post_treatment', y=["FLR1-perc_t", "FLR1-perc_yet"])


