#%%
#, Read in the required functions
import os
import matlab.engine
import single_cell_reloc.global_functions.global_variables as glv
#%%
def locate_tracking_file(): #* This should work regardless of where it grabbed from
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



if __name__ = "__main__":
	Global_variables = glv.global_manager()
	matlab_file, t = locate_tracking_file()
	r = 1
	while t != True:
		if r <= 3:
			exit()
		matlab_file, t = locate_tracking_file()
	



# #, This is an implementation which while elegant, may not be worth it.
# # print(os.path.abspath(__file__)) #, NOTE: This line will not run in the interactive window
# package_root =  os.path.abspath(os.path.join(os.path.abspath(__file__),"../.."))
# print(package_root)

# # absolute_path = os.path.dirname(__file__)
# relative_path = "\MATLAB"
# script_dir_path = os.path.join(package_root, relative_path)
# print(script_dir_path)

# #%%
# eng = matlab.engine.start_matlab()
# # path = "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc/MATLAB"
# eng.cd(script_dir_path, nargout=0)
# returned =eng.Test(nargout = 1)
# print(returned)
# eng.quit()

# #%%
# eng.cd(r'../..single_cell_reloc/single_cell_reloc/MATLAB', nargout=0)
# #%%
# eng.Test(nargout=0)
# eng.quit()