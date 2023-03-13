#%%
#, Read in the required functions
import os
import matlab.engine
import single_cell_reloc.global_functions.global_variables as glv
import single_cell_reloc.ImageJ as ImageJ
from joblib import Parallel, delayed
import pandas as pd
from datetime import date

#%%
def locate_tracking_file(): #* This should work regardless of where the function is called
	if __name__ == '__main__':
		parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
		relative_path = "MATLAB/test.m"
		matlab_file = os.path.join(parent_dir, relative_path)
		test = os.path.exists(matlab_file)
		return(matlab_file, test)
	else: #* This is for when the function is being called from the master file, so it should only go up two levels instead of 3.
		parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
		relative_path = "MATLAB/test.m"
		matlab_file = os.path.join(parent_dir, relative_path)
		test = os.path.exists(matlab_file)
		return(matlab_file, test)

def pre_Tracking(subset = False, subset_by = '', subset_collection = ''):
	subset_by = Global_variables["subset_by"]
	subset = Global_variables["subset_collection"]
	os.chdir(Global_variables['microfluidics_results'])
	imgIndex = pd.read_paraquet('imgIndex.paraquet')

	if subset == True:
		imgIndex = imgIndex.loc[(imgIndex[subset_by].isin(subset_collection))]
	else:
		pass
	Pos_frame_list = imgIndex[["Unique_pos", "Unique_frame",
                           "Date"]].drop_duplicates().set_index(["Unique_pos"])

	orgAllmasks = pd.read_paraquet("orgAllmasks.paraquet")
	seg_dirDF = orgAllmasks[["Unique_pos", "Path", "Date"]].set_index(["Unique_pos"])

	raw_dirDF = imgIndex[(imgIndex["Channel"] == "Sub") | (
    imgIndex["Channel"] == "out")][["Unique_pos", "Path", "Channel"]].set_index(["Unique_pos"])
	pos_index = raw_dirDF.reset_index()["Unique_pos"].unique()

	global current_date
	current_date = str(date.today())

	return(Pos_frame_list, seg_dirDF, pos_index, raw_dirDF)


#! This is the older version which does not have the tracking
def tracking_all_para(p):#? , Pos_frame_list, seg_dirDF, raw_dirDF, current_date):
    try:
        eng = matlab.engine.start_matlab()
        pos = p
        posFingerprint = "position" + p[p.find("p")+1:]
        length = len(Pos_frame_list.loc[p, ])
        try:
            g = seg_dirDF.loc[p, ["Path"]].iloc[0]
        except KeyError:
            print(f"KeyError at {p}")
        else:
            try:
                g = seg_dirDF.loc[p, ["Path"]].iloc[0].values[0]
            except AttributeError:
                return(p)
            # print(p)
            path_seg = os.path.dirname(g)
            path_raw = os.path.dirname(
                raw_dirDF.loc[p, ["Path"]].iloc[0].values[0])
            BFchan = raw_dirDF.loc[p, ["Channel"]].iloc[0].values[0]
            print(pos, current_date, posFingerprint,
                  length, path_seg, path_raw, BFchan)
            if BFchan == 'Sub':
                eng.TrackALL(pos, current_date, posFingerprint,
                             length, path_seg, path_raw, nargout=0)
            elif BFchan == 'out':
                eng.TrackALL_out(pos, current_date, posFingerprint,
                                 length, path_seg, path_raw, nargout=0)
        eng.quit()
    except:
        return(p)


if __name__ == "__main__":
	Global_variables = glv.global_manager()
	matlab_file, t = locate_tracking_file()
	r = 1
	while t != True:
		if r <= 3:
			exit()
		matlab_file, t = locate_tracking_file()

	Pos_frame_list, seg_dirDF, pos_index, raw_dirDF = pre_Tracking(subset=Global_variables['subset'], subset_by=Global_variables['subset_by'], subset_collection= Global_variables['subset_collection']) #* These will now be global variables which can be accessed

	Par
	ImageJ.composite_manager()




	fails_captures = Parallel(n_jobs=Global_variables['cpu_se']//2, verbose = 100)(delayed(tracking_all_para)(p) for p in pos_index)
	fails = pd.Series(fails_captures)
	fails = fails[~fails.str.contains("KeyError")]

	for p in fails: #for the failures, try again in loop
		tracking_all_para(p)
	del Pos_frame_list
	del seg_dirDF
	del pos_index
	del raw_dirDF
