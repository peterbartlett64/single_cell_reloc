import os
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import psutil as p

def rename_files(directory, mem_avail):
	for filename in os.listdir(directory):
		if filename.endswith(".parquet"):
			file_path = os.path.join(directory, filename)
			os.system
			if os.path.getsize(file_path) < (mem_avail/5):
				smaller_files_rename(file_path)
			else:
				big_file_rename(file_path)

def smaller_files_rename(file_path):
	df = pd.read_parquet(file_path)
	if 'Loc_score' in df.columns:
		df.rename(columns={'Loc_score': 'ReLocScore'}, inplace=True)
		df.to_parquet(file_path)
		print(f"Renamed column in {file_path}")

def big_file_rename(file_path):
	pq_raw = pq.read_table(source=file_path)
	print(f'loaded {file_path}')
	pq_raw.rename_columns(['ReLocScore' if x == 'Loc_score' else x for x in pq_raw.column_names])
	print(f'Renamed column in {file_path}')
	pq.write_table(pq_raw, file_path)


if __name__ == '__main__':
	mem = p.virtual_memory()
	mem_avail = mem.available
	rename_files(input("Enter the directory path: "), mem_avail)
	big_file_rename('D:\\Renaming\\Combined_by_perc\\Quant_ALL.parquet')