#%%
from joblib import Parallel, delayed
import pandas as pd
import math
from scipy import stats
from joblib import Parallel, delayed
from datetime import date
import os
import psutil as p
import math
from glob import glob
import single_cell_reloc_parquet.global_functions.global_variables as gv
from plotnine import ggplot, geom_line, aes, scale_alpha_discrete
#%%

os.chdir("E:/ALL_FINAL/99th_percentile")
tc_pen =pd.read_parquet("All_pos_pt_t_less_percentages_melt.parquet")
frame_pen =pd.read_parquet("All_pos_pt_percentages_melt.parquet")
#%%
tc_pen["Protein"] = pd.Series(tc_pen["Protein"]).apply(lambda x: x[:x.find("-")])
frame_pen["Protein"] = pd.Series(frame_pen["Protein"]).apply(lambda x: x[:x.find("-")])
#%%
tc_pen.set_index(['Protein', 'Frames_post_treatment'], inplace=True)
frame_pen.set_index(['Protein', 'Frames_post_treatment'], inplace=True)
#%%
merged_pens_t = pd.merge(tc_pen, frame_pen, left_index=True, right_index=True, how = 'inner', suffixes=('', '_y')).reset_index(drop=False)
# merged_pens_t = merged_pens_t[merged_pens_t.columns.drop(list(merged_pens_t.filter(regex='_y')))] #* Drop cols containing _y (ie repeat)
#%%
(ggplot(merged_pens_t, aes(x = "Time_post_treatment	", y= "Percentage_reloc_less", color='Protein'))+
 geom_line()+
 scale_alpha_discrete(guide=False))
#  geom_line(aes = ("Time_post_treamtent", "Perecentage_reloc")))
