#%%
#, This will only run with the Microfluidics_Pipe_highest environment
#%%
import pandas as pd
import os
import imagej
from joblib import Parallel, delayed
import single_cell_reloc.global_functions.global_variables as glv
# import scyjava
#%%
#! Correct this to prompt user for the local instalation of imageJ
ij = imagej.init('C:/Users/pcnba/Desktop/Fiji.app') #* Load the local instalation of imageJ.

#%%
def composite_manger(position):
	os.chdir(Global_Variables["microfluidics_results"])
	imgIndex = pd.read_csv('imgIndex.csv')
	Parallel()


def composite_imageJ():
	ij = imagej.init('C:/Users/pcnba/Desktop/Fiji.app') #* Load the local instalation of imageJ.


if __name__ = "__main__":
	Global_Variables = glv.global_manager()
	composite_manger()