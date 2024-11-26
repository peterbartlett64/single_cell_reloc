# pylint: disable=bad-indentation
#%%
#, Read in the required functions
import os
import matlab.engine
import single_cell_reloc_parquet.global_functions.global_variables as gv
import single_cell_reloc_parquet.ImageJ as ImageJ
from joblib import Parallel, delayed
import pandas as pd
from datetime import date
import single_cell_reloc_parquet.Pre_seg.img_indexer as img_er
import single_cell_reloc_parquet.Pre_track.orgAllmask_er as org_er
import logging
#%%/
logger_tracking = logging.getLogger('tracking_logger')

#%%
def switch_slash(path):
	new = path.replace(os.sep, '/')
	return (new)


#TODO: Modify this so that it can be used in the global variables calls
# def locate_tracking_file(name = __name__, file = __file__): #* This should work regardless of where the function is called from to find the "MATLAB" folder
# 	if name == '__main__':
# 		parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(file))))
# 		relative_path = "MATLAB/test.m"
# 		matlab_file = os.path.join(parent_dir, relative_path)
# 		test = os.path.exists(matlab_file)
# 		return(matlab_file, test)
# 	elif name != '__main__': #* This is for when the function is being called from the master file, so it should only go up two levels instead of 3.
# 		parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(file)))
# 		relative_path = "MATLAB/test.m"
# 		matlab_file = os.path.join(parent_dir, relative_path)
# 		test = os.path.exists(matlab_file)
# 		return(matlab_file, test)

def pre_Tracking(microfluidics_results, subset = False, subset_by = '', subset_collection = ''): #* Prepare trracking information
	#// try: #* Test if the tracker exists and if doesn't create it. Avoid overwriting the files
	#// 	logger_tracking
	#// except NameError:

	#! Fix this
	#Todo: ...
	# global logger_tracking
	# logger_tracking = gv.setup_logging(logger_name='TrackLogger', microfluidics_results= microfluidics_results)

	imgIndex = pd.read_parquet('imgIndex.parquet').reset_index(drop = False)
	# imgIndex = img_er() #Todo: Modify so that it takes in the input variables for image locatoin and microfluids_results
	if subset == True:
		imgIndex = imgIndex.loc[(imgIndex[subset_by].isin(subset_collection))]
	else:
		pass

	Pos_frame_list = imgIndex[["Unique_pos", "Unique_frame",
							"Date"]].drop_duplicates().set_index(["Unique_pos"])

	try:
		orgAllmasks = pd.read_parquet("orgAllmasks.parquet").reset_index(drop = False)
	except FileNotFoundError:
		orgAllmasks = org_er()

	seg_dirDF = orgAllmasks[["Unique_pos", "Path", "Date"]].set_index(["Unique_pos"]) #* This is a call to the global variable

	raw_dirDF = imgIndex[(imgIndex["Channel"] == "Sub") | (
	imgIndex["Channel"] == "out")][["Unique_pos", "Path", "Channel"]].set_index(["Unique_pos"])
	pos_index = raw_dirDF.reset_index()["Unique_pos"].unique()

	global current_date
	current_date = str(date.today())

	return(Pos_frame_list, seg_dirDF, pos_index, raw_dirDF)


#! This is the older version which does not have the stage tracking
def tracking_all_para(p, MATLAB_files, par = True, quant = False):#? , Pos_frame_list, seg_dirDF, raw_dirDF, current_date):
	try:
		eng = matlab.engine.start_matlab()
		eng.cd(MATLAB_files, nargout=0)
		pos = p
		posFingerprint = "position" + p[p.find("p")+1:]
		length = len(Pos_frame_list.loc[p, ])
		try:
			g = seg_dirDF.loc[p, ["Path"]].iloc[0] #Todo: Make this actually the directory instead of a random file it can be derived from
		except KeyError:
			# logger_tracking.error(f"KeyError at {p}")
			return(f"KeyError at {p}")
		else:
			try:
				g = seg_dirDF.loc[p, ["Path"]].iloc[0].values[0]
			except AttributeError:
				# logger_tracking.error(f"Attribute error on {p}")
				return(f"Attribute error on {p}")
			# print(p)
			path_seg = os.path.dirname(g)
			path_raw = os.path.dirname(
				raw_dirDF.loc[p, ["Path"]].iloc[0].values[0])
			BFchan = raw_dirDF.loc[p, ["Channel"]].iloc[0].values[0]
			print(pos, current_date, posFingerprint, length, path_seg, path_raw, BFchan)
			if BFchan == 'Sub':
				eng.TrackALL(pos, current_date, posFingerprint, length, path_seg, path_raw, par, quant, nargout=0) #* Pass the [position, current_date, full position name ("position_dMMDDr[0-9]{1}p[0-9]{6,7}, length, segmentation directory, image directory, parallel boolean, segmentation quant boolean]"
			elif BFchan == 'out':
				eng.TrackALL_out(pos, current_date, posFingerprint, length, path_seg, path_raw, par, quant, nargout=0)
		eng.quit()
	except:
		return(p)

#%%
if __name__ == "__main__":
	# Global_variables = glv.global_manager()
	# cpu_se = Global_variables['cpu_se']
	Global_variables = {'microfluidics_results': gv.slash_switch(input("microfluidics_results?")),
						'subset': True,
						'subset_by':'Unique_pos',
						'subset_collection': ['d0214r1p010200', 'd0214r1p010300', 'd0214r1p030200', 'd0214r1p030300', 'd0214r1p040200', 'd0214r1p040300', 'd0214r1p050200', 'd0214r1p050300', 'd0214r1p070200', 'd0214r1p070300', 'd0214r1p080200', 'd0214r1p080300', 'd0214r1p090200', 'd0214r1p090300', 'd0214r1p110200', 'd0214r1p110300', 'd0214r1p120200', 'd0214r1p120300', 'd0214r1p130200', 'd0214r1p130300', 'd0214r1p150200', 'd0214r1p150300', 'd0214r1p160200', 'd0214r1p160300', 'd0214r2p010200', 'd0214r2p010300', 'd0214r2p030200', 'd0214r2p030300', 'd0214r2p040200', 'd0214r2p040300', 'd0214r2p050200', 'd0214r2p050300', 'd0214r2p070200', 'd0214r2p070300', 'd0214r2p080200', 'd0214r2p080300', 'd0214r2p090200', 'd0214r2p090300', 'd0214r2p110200', 'd0214r2p110300', 'd0214r2p120200', 'd0214r2p120300', 'd0214r2p130200', 'd0214r2p130300', 'd0214r2p150200', 'd0214r2p150300', 'd0214r2p160200', 'd0214r2p160300', 'd0215r1p010200', 'd0215r1p010300', 'd0215r1p030200', 'd0215r1p030300', 'd0215r1p040200', 'd0215r1p040300', 'd0215r1p050200', 'd0215r1p050300', 'd0215r1p070200', 'd0215r1p070300', 'd0215r1p080200', 'd0215r1p080300', 'd0215r1p090200', 'd0215r1p090300', 'd0215r1p110200', 'd0215r1p110300', 'd0215r1p120200', 'd0215r1p120300', 'd0215r1p130200', 'd0215r1p130300', 'd0215r1p150200', 'd0215r1p150300', 'd0215r1p160200', 'd0215r1p160300']
	}
	# matlab_file, t = locate_tracking_file()
	microfluidics_results = Global_variables['microfluidics_results']
	os.chdir(microfluidics_results)
	# try:
	# 	os.chdir(microfluidics_results)
	# except:
	# 	microfluidics_results = gv.slash_switch(input("microfluidics_results?"))
	# 	os.chdir(microfluidics_results)

	matlab_file = "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/MATLAB"
	# matlab_file = "C:/Users/Nikon/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/MATLAB"

	# r = 1
	# while t != True:
	# 	if r <= 3:
	# 		exit()
	# 	matlab_file, t = locate_tracking_file()

	Pos_frame_list, seg_dirDF, pos_index, raw_dirDF = pre_Tracking(microfluidics_results = Global_variables['microfluidics_results'], subset=Global_variables['subset'], subset_by=Global_variables['subset_by'], subset_collection= Global_variables['subset_collection']) #* These will now be global variables which can be accessed
	# Pos_frame_list, seg_dirDF, pos_index, raw_dirDF = pre_Tracking() #* These will now be global variables which can be accessed

	# imgIndex = pd.read_parquet("imgIndex.parquet").reset_index(drop = False)


	# cpu_se = Global_variables['cpu_se']
	cpu_se = os.cpu_count()

	eng = matlab.engine.start_matlab() #* Start a temporary MATLAB server to check that the unit test works

	eng.cd(matlab_file, nargout=0)
	answer = eng.Test()
	if answer == 2:
		print('Found MATLAB files!')
	else:
		print('Did not find MATLAB files')
	eng.quit()


	# def test_func():
	# 	r = random.randint(0,100)
	# 	if r> 80:
	# 		logger_tracking.error('test')
	# 	else:
	# 		pass
	# 	return(r)

	# ImageJ.composite_manager() #! This functionality is if the functions must be rerun with altered quantification (on denoised)
	if len(pos_index) <= cpu_se*2: #* This is an choice based on experience
			for p in pos_index:
				tracking_all_para(p, MATLAB_files = matlab_file)
	else:
		fails_captures = Parallel(n_jobs=cpu_se//2, verbose = 100)(delayed(tracking_all_para)(p, MATLAB_files = matlab_file, ) for p in pos_index)
		fails = pd.Series(fails_captures)
		fails_key = fails[~fails.str.contains("KeyError")]
		fails_key.to_csv('Tracking_fails_key')

		for p in fails: #for the failures
			tracking_all_para(p, MATLAB_files = matlab_file)
		del Pos_frame_list
		del seg_dirDF
		del pos_index
		del raw_dirDF
