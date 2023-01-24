#%%
import matlab.engine
import os
from pathlib import Path
import sys

#%%
#, This is an implementation which while elegant, may not be worth it.
# print(os.path.abspath(__file__)) #, NOTE: This line will not run in the interactive window
package_root =  os.path.abspath(os.path.join(os.path.abspath(__file__),"../.."))
print(package_root)

# absolute_path = os.path.dirname(__file__)
relative_path = "\MATLAB"
script_dir_path = os.path.join(package_root, relative_path)
print(script_dir_path)

#%%
eng = matlab.engine.start_matlab()
# path = "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/single_cell_reloc/single_cell_reloc/MATLAB"
eng.cd(script_dir_path, nargout=0)
returned =eng.Test(nargout = 1)
print(returned)
eng.quit()

#%%
eng.cd(r'../..single_cell_reloc/single_cell_reloc/MATLAB', nargout=0)
#%%
eng.Test(nargout=0)
eng.quit()