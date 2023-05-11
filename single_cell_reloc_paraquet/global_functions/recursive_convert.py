from glob import glob
import pandas as pd
import single_cell_reloc_parquet.global_functions.global_variables as gv

#, Function to recursively convert csv files to parquet files
def recursive(root):
	for path in glob(root + '**', recursive=True):
		if path.endswith('.csv'):
			df = pd.read_csv(path)
			df.to_parquet(path.replace('.csv', '.parquet'))
			print(path)

if __name__ == '__main__':
	Global_variables = gv.global_manager()
	print(f'Converting csv files to parquet files recursively within {Global_variables.microfluidics_results}')
	recursive(Global_variables.microfluidics_results)
else:
	recursive(Global_variables.microfluics_results)
	print(f'Converting csv files to parquet files recursively within {Global_variables.microfluidics_results}')
