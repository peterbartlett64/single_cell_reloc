import pandas as pd
import os




for f in dfs:
	exec(f"{f} = pd.read_parquet({f}.parquet)")


def check_present(dataframes, pseudo_index , pseudo_col == 'Unique_frame'):
	res_df = pd.DataFrame([])
	for df in dataframes:
		df = pd.read_parquet(f"{df}.parquet") #* convert from a df name to the actual dataframe
		t = df.loc[df[pseudo_col] == pseudo_index]
		if len(t >= 1):
			return(True)
		else:
			return(False)
	return()