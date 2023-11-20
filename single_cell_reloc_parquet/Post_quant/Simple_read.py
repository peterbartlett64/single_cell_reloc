#%%
import pandas as pd
import math
from scipy import stats
from joblib import Parallel, delayed
from datetime import date
import os
import psutil as p
import math
from datetime import datetime
from joblib import Parallel, delayed
from glob import glob
import numpy as np

#%%
today = str(date.today())

def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
	new = path.replace(os.sep, '/')
	return (new)

#* This was for comparing between the flattened array vs index = 0 percentiles. Confirmed that it was not an issue and values were indeed indentical
#. Tested with >>
#! D:\ALL_FINAL\Raw_quant\2023-09-14\Quantification_d0224r1p060300f0059_2023-09-14.parquet
#! E:\Microfluidics\RESULTS\2023-10-08\Quantification_d0224r1p060300f0059_2023-10-08.parquet


path_read1 = slash_switch(input('Path to read'))
path_read2 = slash_switch(input('Next Path to read'))
#%%
file1 = pd.read_parquet(path_read1).set_index(['Cell_Barcode', 'ImageID'])
#%%
file2 = pd.read_parquet(path_read2).set_index(['Cell_Barcode', 'ImageID'])

#%%
combined = pd.merge(file1, file2, left_on= 'Cell_Barcode', right_on='Cell_Barcode', suffixes=('file1', 'file2'))
columns = pd.Series(combined.columns)
columns = columns.apply(lambda x: x[:-5]).unique()
#%%
difference_map = pd.DataFrame([])
for c in columns:
	difference_map[c] = np.where((combined[f"{c}file1"] <= combined[f"{c}file2"]), True, False)
	not_l = difference_map.loc[difference_map[c] == False]

	if len(not_l) > 0:
		print(c)
	else:
		pass