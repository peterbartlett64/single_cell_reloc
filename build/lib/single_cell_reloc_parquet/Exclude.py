#%%
import pandas as pd
import os
import global_variables.slash_swtich as slash_switch

def Exclude():
	name = input("What file name?")
	excl_type = input("On what level is exclusion?")
	exclude = input("What runs to exclude?")
	exclude = str.split(exclude, ", ")

	df = pd.read_parquet(f"{name}.parquet")
	df_subset = df.loc[~(df[excl_type].isin(exclude))]
	subset_df = df.to_parquet(f"{name}_subset.parquet")
	return("Subset complete!")

if __name__ == "__main__":
	path = input("What is the path?")
	path = slash_switch(path)
	os.chidr(path)
	Exclude()

