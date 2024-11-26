#%%

#Import requried Packages
# import h5py
import subprocess
#import dask.dataframe as daf
# from numba import njit, vectorize, cuda
import time
from tkinter import E
from scipy.fft import fft, ifft
import tables
import itertools
import scipy.sparse as spd
import pandas as pd
import datetime
import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
from matplotlib import image as mpimg
# import seaborn as sns
import scipy.io as scipio
from PIL import Image
import glob
import os
import ntpath as ntp
import csv
import statistics
import pymongo
import json
import multiprocessing as mp
from scipy.io import loadmat
import time
from multiprocessing import Pool

def slash_switch(path):
    new = path.replace(os.sep, '/')
    return (new)

#REFER to README.txt for setup protocols
pn = os.cpu_count()# to get the cpu core count
print ("Pipeline will run with", pn, "process nodes (number of system cores)")

# temp_store = "C:/Users/Peter/Desktop/DATA_TEST"

# code = "C:/Users/Peter/Dropbox (Grant Brown's Lab)/Peter Bartlett Data/Code"

# #Set the location of the experimental flow file/condtion information
# locationC = "F:/Microfluidics/DATA_TEST/20210112_BaselMicrofluidics_GeneMap.xlsx" 

#Fill in the fluorescent information
# Fluorescent = {'GFP': 'PROT', 'mKa': 'BHY175_Myo1', 'mKO': 'BHY131_Myo1'} 

#Set the image path to the location of the raw image files (Largest directory)


# Define experimental information
Experimental_info = {
    "timepoint_space": 7.5,
	"TracX_version": "TracX",
    "cell_tracking": True,
    "flourescent": True,
    "data_subset": False,
    "windows_strings": True,
	"seg_rerun": False #This is set so that if there is a segmenation frame found, it is not rerun. If parameters or paths have changed, set to True
    }

roots = {
	"F_drive_root" : "F:/",
	"D_drive_root" : "D:/",
	"xps_usr" : "C:/Users/pcnba",
	"nikon_usr" : "C:/User/Nikon"
}

results_dir = "RESULTS"
analyze_dir = "Analyze"
RUN_seg_dir = "RUN_seg_Prog_Lib"

ask_loc = True
while ask_loc == True:
	Comp_sel = str(input("Which computer? [XPS or Nikon]"))
	if Comp_sel == "XPS":
		current_code = "C:/Users/pcnba/Dropbox (Grant Brown's Lab)/Peter Bartlett Data/Code/Current"
		microfluidics_results = "D:/Microfluidics/RES_N_ULTS"
		image_path = "D:/Microfluidics/Analyze"
		RUN_seg_lib = os.path.join(current_code, RUN_seg_dir)
		TracX_dir = f"/{Experimental_info['TracX_version']}/src/"
		TracX_path = os.path.join(RUN_seg_lib, TracX_dir)

		# microfluidics_results = os.path.join(roots["D_drive_root"], results_dir) # This is another method of implementation whcih I will do in the future
		ask_loc = False
	elif Comp_sel == "Nikon":
		current_code = "C:/User/Nikon/Desktop/Current_Both"
		microfluidics_results = "F:/Microfluidics/RESULTS"
		image_path = "F:/Microfluidics/Analyze"
		RUN_seg_lib = os.path.join(current_code, RUN_seg_dir)
		TracX_dir = f"/{Experimental_info['TracX_version']}/src/"
		TracX_path = os.path.join(RUN_seg_lib, TracX_dir)

		ask_loc = False
	else:
		ask_loc = True

#%%
image_path = str(input("What is the root folder containing all source images?"))
if Experimental_info["windows_strings"] == True:
    image_path = slash_switch(image_path)
else:
    pass
#"F:/Microfluidics/Analyze" # This would normally be a folder containing all date folders with the images aquired at those days {[All]/[2018_x]/[images...]}

#Set the Results directory. This is where the CellX will output to and where the Master hdf5 file will be stores

# microfluidics_results = str(input("Where are the results being stored? This folder should also contain the condition informaton file that relates postion to: protein, type of localization change, time of treatment, and any notes. Shape of file is outlined in documentation"))
# if Experimental_info["windows_strings"] == True:
#     microfluidics_results = slash_switch(microfluidics_results)
# 	try:
# 		os.mkdir(microfluidics_results)
# 	except FileExistsError:
# 		pass
# else:
#     pass

post_quant = str(input("Where would you like the quantificaiton to be outputted to? (By default, it will go to todays date). The folder will be created if it does not exist"))

if Experimental_info["windows_strings"] == True:
    post_quant = slash_switch(post_quant)
else:
    pass

# "F:/Microfluidics/RESULTS" #This is a folder containing all image masks and other output from CellX. It should be of shape ALL/Date/[GFP],[RFP].../  AND the HDF5 file in my case

os.chdir(current_code)
from Pre_seg import imgIndex_er, xml_generator #Single core, single core
from Pre_seg import segment_manager # Multi-core

from Pre_track import mask_indexer.orgmask_er as orgmask_er # Single core
from Pre_track import mask_indexer.trackmasker as trackmaster # Single core
if Comp_sel == "XPS":
	from Pre_track import tracker.tracker_xps as tracker
elif Comp_sel == "Nikon":
	from Pre_track import tracker.tracker_nikon as tracker #Multi-core on number of positions

from Post_track import mask_expand_er #Multi-core

from Quantification_funcs import quantification #Multi-core

from Post_quant import post_quantification #Multi-core - mixed on: number of frames, number of positions  
from Post_quant import primary_visualizer #Multi-core on number of postions

#Set the location of the CellX program and the MCR version folder. On Linux, the CellX will be an .sh file
# CellXc = "C:/Users/Peter/CellX/Windows_7_g1.12_c1.12/CellX.exe"
#MCRc = "C:/Program Files (x86)/MATLAB/MATLAB Compiler Runtime/v716# %%

#Get the ask_loc fucntion to work better. Right now the ask_loc function does not properly deliminate between requirments. Manually enter requirement in the pandas selection parameters later on for now. 
if Experimental_info["data_subset"] == True:
    subdat = str(input("Enter position barcodes to be analyzed (comma seperated and in the form dMMDDpPosition)"))

    instances = subdat.split(", ")
    print(f"Analysis will be performed for positions which correspond to the following pos_barcodes: {instances}")
else:
    pass

# try: subdat
# except NameError: pass
# else: del subdat

#Run the image indexing of available files
imgIndex_er()

#Create the targetting file to be handed to the segmentation program
orgmask_er()
xml_generator() #This xml_generator will check which files have already been segmented. From the above function run, and skip these files to save time. 
#(This is for the case that the pipeline is run in parts)

#Run the segmentation. The manager takes the input from xml_generator and
segment_manager()

#Primary generated masks are used in tracking targeting
orgmask_er()

#Tracker program is automaticallly run and then the output masks are indexed
track_manager(instances, imgIndex, orgAllmasks)

#Tracked masks are expanded to capture the cell membrane and budneck more effectively. Masks are epanded by 1/2 watershed distance (4 pixels)
mask_expand_er()

#Parallel quantification of microscope images by frame
quantification()

#Run the post_quantifcation pipeline
post_quantification()

#Output figures from the data collected
primary_visualizer() #This will prepare some of the default visualizations

# %%
