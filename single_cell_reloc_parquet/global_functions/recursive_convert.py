#%%
from glob import glob
import pandas as pd
import os
import single_cell_reloc_parquet.global_functions.global_variables as gv

#, Function to recursively convert csv files to parquet files
def recursive_convert(root):
	os.chdir(root)
	for path in glob(root + '**', recursive=True):
		if path.endswith('.csv'):
			df = pd.read_csv(path)
			df.to_parquet(path.replace('.csv', '.parquet'))
			print(path)

if __name__ == '__main__':
	# Global_variables = gv.global_manager()
	# microfluidics_results = Global_variables.microfluidics_results
	microfluidics_results = gv.slash_switch(input("Where is microfludics_results?"))
	print(f'Converting csv files to parquet files recursively within {microfluidics_results}')
	recursive_convert(microfluidics_results)
else:
	recursive_convert(Global_variables.microfluics_results)
	print(f'Converting csv files to parquet files recursively within {Global_variables.microfluidics_results}')
