import pandas as pd
import os
import pathlib
import subprocess
import single_cell_reloc_parquet.Pre_seg.xml_generator as xml_generator
import single_cell_reloc_parquet.Pre_track.orgAllmask_er as orgAllmask_er
import single_cell_reloc_parquet.global_functions.global_variables as gv


#, Create the conda environment. This is a one-time operation.
conda_env = 'pycellxworkflow'
subprocess.run(['conda', 'env', 'create', '-f', 'segmentation_manager.yml', '-n', conda_env])

# def segmentation_helper(Global_variables):
# 	fileseries = Global_variables['fileseries']
# 	parameters = Global_variables['parameters']
# 	cellx_location = Global_variables['cellx_location']
# 	cpu_count = Global_variables['cpu_se']
# 	return(fileseries, parameters, cellx_location, cpu_count)

def segmentation_manager_parralllel(fileseries, parameters, cellx_location, cpu_count):
	os.chdir(pathlib.Path(__file__).parent.absolute())
	#* Run the segmentation manager
	os.chdir(os.path.join(pathlib.Path(__file__).parent.absolute(), 'RUN_SegProgLib'))
	subprocess.run(['conda', 'run', '-n', conda_env, 'snakemake', '--config', f'fileseries={fileseries}', f'calibration={parameters}', f'cellx={cellx_location}', '--cores', cpu_count, '--latency-wait 60', '--keep-going'])

def segmentation_manager_single_run(fileseries_small, parameters, cellx_location):
	os.chdir(os.path.join(pathlib.Path(__file__).parent.absolute(), 'RUN_SegProgLib'))

	#! This is a single call that will run the remainder of segmentation single core all RAM
	subprocess.run([f'{cellx_location}', f'{parameters}', '-s', f'{fileseries_small}'])

def seg_manager(Global_variables)->None:
	#* Moved away from the helper function to just have the global variables passed in
	# fileseries, parameters, cellx_location, cpu_count = segmentation_helper(Global_variables = Global_variables)

	#* Run the segmentation manager
	segmentation_manager_parralllel(fileseries = Global_variables['fileseries'], parameters = Global_variables['parameters'], cellx_location = Global_variables['cellx_location'], cpu_count = Global_variables['cpu_se'])

	org_mask_completed_file = orgAllmask_er.orgAllmask_er(Global_variables = Global_variables, pre = True)
	orgAllmasks_completed = pd.read_csv(os.path.join(Global_variables['microfluidics_results'], 'orgAllmasks_completed.csv'))
	imgIndex = pd.read_csv(os.path.join(Global_variables['microfluidics_results'], 'imgIndex.csv'))

	#* Check what is missing
	missing = imgIndex[~imgIndex['imgIndex'].isin(orgAllmasks_completed['imgIndex'])]

	xml_file_name = xml_generator(Global_variables = Global_variables, imgIndex = missing, retry = True)
	fileseries_small = Global_variables['microfluidics_results'] + xml_file_name

	segmentation_manager_single_run(fileseries_small = fileseries_small, parameters = Global_variables['parameters'], cellx_location = Global_variables['cellx_location'])
	print('Segmentation Complete')

if __name__ == '__main__':
	Global_variables = gv.global_manager()
	seg_manager(Global_variables = Global_variables)