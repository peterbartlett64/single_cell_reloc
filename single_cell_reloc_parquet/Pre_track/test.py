#%%
import os
import matlab.engine
# from pathlib import Path
# import sys

# #%%
# import os

# parent_directory = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
# matlab_file = os.path.join(parent_directory, 'matlab_file.m')
# print(matlab_file)
def locate_tracking_file()
	parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	relative_path = "MATLAB/test.m"
	matlab_file = os.path.join(parent_dir, relative_path)
	t = os.path.exists(matlab_file)

print(matlab_file, t)
# Do something with the matlab_file
# for example, you can use the scipy library to read the matlab file

# from scipy.io import loadmat
# matlab_data = loadmat(matlab_file)



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