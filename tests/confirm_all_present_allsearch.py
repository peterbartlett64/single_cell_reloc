#%%
####################################################################################         MAKE SURE TO COPY THIS IS MARKED SECTION   ####################################################################################################
from xml.dom import xmlbuilder
from numpy.lib.utils import info
import pandas as pd
import os
from scipy.io import loadmat
import datetime
import csv
from operator import pos
from xml.etree.ElementTree import Element, SubElement, Comment, TreeBuilder, tostring
import datetime
from xml.dom import minidom
import os
import pandas as pd
import numpy as np
import datetime
import warnings
from pandas.core.indexing import IndexingError
from datetime import date
import os
from re import sub
import matlab.engine
import pandas as pd
from joblib import Parallel, delayed
import time
#%%

date_today = datetime.date.today()
current_date = date_today


def slash_switch(path): ##k This function is currently unused but could be usefull in the future for the cwd setting
    new = path.replace(os.sep, '/')
    return (new)

# current_code = str(input("Current code"))
# current_code = slash_switch(current_code)
# current_code = "C:/Users/pcnba/Dropbox (Grant Brown's Lab)/Peter Bartlett Data/Code/Current"
current_code = "C:/Users/Nikon/Documents/Peter/Current_Both"

os.chdir(current_code)
# import track_manager_MAY as track_manager
# import xmlgen_MAY as xmlbuilder
try:
    prompt = input(f"{microfluidic_results}? y/n")
    if prompt == "n":
        microfluidics_results = str(input("Microfluidics results folder"))
    else:
        pass
except NameError:
    microfluidics_results = str(input("Microfluidics results folder"))

# microfluidics_results = "F:/Microfluidics/RESULTS" ######################################## Change this. Must be hard coded as there is no memmory once anaconda reactivates
microfluidics_results = slash_switch(microfluidics_results)
microfluidic_results = microfluidics_results

try:
    fluor
except NameError:
    fluor = str(input("Run CellX fluorescence? (True or False)"))

file_mod = str(input("Is this for a specific case? What should the file modifier be?"))

# Run_segProgLib = str(input("Where is Run_segProgLib?"))

Run_segProgLib = str(input("where is RUN_segProgLib"))
Run_segProgLib = slash_switch(Run_segProgLib)
# Run_segProgLib = "C:/Users/pcnba/Dropbox (Grant Brown's Lab)/Peter Bartlett Data/Code/Current/RUN_segProgLib"
#%%
# M_instances = str(input("What the instances that you wnat to test within? This is to avoid segmenting and tracing files that are not immediatly needed"))
seg = True
Track = True
#Todo Fix the input subset
if GlobalVariables['Subset'] == True:
    instances = GlobalVariables['Subset_instances'] #! Fix this
else:
    pass

#################################################################################################
#%%

os.chdir(microfluidics_results)
orgmaskpaths = []
count = 0
fluors = "GFP_mKO_mKa"
exclude = list(['No_fluor', 'No_flour', 'mKO_mKa'])
for root, dirs, files, in os.walk(os.getcwd(), topdown= True):
    [dirs.remove(d) for d in list(dirs) if d in exclude]
    for name in files:
        if name.endswith(".mat") and name.startswith("mask"):
            orgmaskpaths.append({'Path': os.path.join(root, name)})
            count = count + 1
            print(count, end="\r")
        else:
            pass

def f_Pos_mask(z):
    start = z.find('p') #Note: The shift forward is dependant upon how you write out position
    end = z.find('GFP')-1  #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
    pos = z[start:end]
    # if len(pos) <= 8:
    #     return(pos)
    # else:
    #     pass
    return(pos)


def mask_inf(m):
    num1 = m.find('ask_')+5
    num2 = m.find('.mat')
    num_ = ("f"+ m[num1:num2])
    return(num_)

def mask_hour(m):
    hourS = m.find("_")+2
    hourE = m.find("_p")
    return(m[hourS:hourE])

def mask_barcode(m):
    start =m.find('d02')
    end = m.find('GFP') -1
    return (m[start:end])

def f_maskdate(x):
    dstart = x.find('d0')
    dend = dstart + 5
    expdate = (x[dstart:dend])
    return(expdate)

orgAllmasks = pd.DataFrame(orgmaskpaths)

def f_non_sync(p):
    substring = ".sync"
    if substring in p:
        return(None)
    else:
        return(p)

orgAllmasks["Path"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_non_sync)
orgAllmasks.dropna(inplace=True)

orgAllmasks["Date"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_maskdate)
orgAllmasks["Position_ID"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_Pos_mask)
orgAllmasks["MaskNum"] = pd.Series(orgAllmasks.iloc[:,0]).apply(mask_inf)
# Allmasks["Hour/Run_t"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_hour)
# Allmasks["Unique_pos_hour"] = Allmasks['Date'] + "h" + Allmasks["Hour/Run_t"] + Allmasks['Position_ID']

orgAllmasks["Unique_pos"] = pd.Series(orgAllmasks.iloc[:,0]).apply(mask_barcode)
orgAllmasks["Unique_frame"] = orgAllmasks["Unique_pos"] + orgAllmasks['MaskNum'] #By performing this way, can avoid one extra step per index

del orgmaskpaths

orgAllmasks.dropna(inplace=True) # temp drop because there are some extra segmentations
orgAllmasks.reset_index(inplace=True, drop = True)

orgAllmasks.set_index(["Unique_frame"], inplace = True)
orgAllmasks

maskpaths = []
count = 0
for root, dirs, files, in os.walk(os.getcwd()):
    for name in files:
        if name.endswith(".mat") and name.startswith("track_mask"):
            maskpaths.append({'Path': os.path.join(root, name)})
            count = count + 1
            print(count, end="\r")
        else:
            pass

# def f_Pos_mask(z):
#     start = z.find('RESULTS') + 8 #Note: The shift forward is dependant upon how you write out position
#     end = z.find('GFP')-1  #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
#     pos = z[start:end]
#     if len(pos) <= 16:
#         return(pos)
#     else:
#         pass

def mask_inf(m):
    num1 = m.find('ask_')+5
    num2 = m.find('.mat')
    num_ = ("f"+ m[num1:num2])
    return(num_)

def mask_hour(m):
    hourS = m.find("_")+2
    hourE = m.find("_p")
    return(m[hourS:hourE])

def mask_barcode(m):
    start =m.find('d02')
    end = m.find('GFP') -1
    return (m[start:end])

def f_maskdate(x):
    dstart = x.find('d0')
    dend = dstart + 5
    expdate = (x[dstart:dend])
    return(expdate)

Allmasks = pd.DataFrame(maskpaths)

def f_non_sync(p):
    substring = ".sync"
    if substring in p:
        return(None)
    else:
        return(p)

Allmasks["Path"] = pd.Series(Allmasks.iloc[:,0]).apply(f_non_sync)
Allmasks.dropna(inplace=True)

Allmasks["Date"] = pd.Series(Allmasks.iloc[:,0]).apply(f_maskdate)
# Allmasks["Position_ID"] = pd.Series(Allmasks.iloc[:,0]).apply(f_Pos_mask)
Allmasks["MaskNum"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_inf)
# Allmasks["Hour/Run_t"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_hour)
# Allmasks["Unique_pos_hour"] = Allmasks['Date'] + "h" + Allmasks["Hour/Run_t"] + Allmasks['Position_ID']

Allmasks["Unique_pos"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_barcode)
Allmasks["Unique_frame"] = Allmasks["Unique_pos"] + Allmasks['MaskNum'] #By performing this way, can avoid one extra step per index

# del maskpaths

# Allmasks.dropna(inplace=True) # temp drop because there are some extra segmentations
Allmasks.reset_index(inplace=True, drop = True)

Allmasks.set_index(["Unique_frame"], inplace = True)
Allmasks

orgAllmasks.to_csv('orgAllmasks.csv')
Allmasks.to_csv('Allmasks.csv')

#%%
os.chdir(microfluidic_results)
imgIndex = pd.read_csv('imgIndex.csv')
imgIndex_temp = imgIndex.rename(columns={"Path" : "image Path"})
imgIndex_temp = imgIndex_temp.loc[(imgIndex_temp["Channel"] == "Sub") | (imgIndex_temp["Channel"] == "out")]
frame_list = imgIndex_temp["Unique_frame"].unique()
imgIndex_temp.set_index("Unique_frame", inplace = True, drop=False)

orgAllmasks = pd.read_csv("orgAllmasks.csv")
orgAllmasks_temp = orgAllmasks[["Unique_frame", "Path"]].rename(columns={"Path": "org_Mask Path", "Unique_frame" : "Unique_frame_org"})
orgAllmasks_temp.set_index(["Unique_frame_org"], inplace= True)

Allmasks = pd.read_csv('Allmasks.csv')
Allmasks_temp = Allmasks[["Unique_frame", "Path"]].rename(columns={"Path": "Track_Mask Path", "Unique_frame" : "Unique_frame_Track"})
Allmasks_temp.set_index(["Unique_frame_Track"], inplace = True)
#%%
merged = pd.merge(imgIndex_temp, orgAllmasks_temp, left_index= True, right_index= True, how = 'left') # This assumes that no images were lost since segmentation
#%%
# if M_instances == "All"
# 	pass

# else:
# merged = merged.loc[merged["Unique_pos"].isin(M_instances)] # Substet to just the files that are immediately needed
merged = pd.merge(merged, Allmasks_temp, left_index= True, right_index= True, how = 'left')# Again, this assumes that no images were lost since tracking

#%%
merged_na = merged.loc[(merged["org_Mask Path"]).isna() | (merged["Track_Mask Path"].isna())] # If there is a missing mask (either from seg or tracking), it will be labelled here
#%%
merged_na = merged_na.drop(columns="Unique_frame")
merged_na.to_csv("Files_missing_df.csv") #output the files that are missing
#%%

org_Mask_missing = merged_na.loc[merged_na["org_Mask Path"].isna()].reset_index(drop = False)
org_Mask_missing_pos = org_Mask_missing["Unique_pos"].unique()

Track_Mask_missing = merged_na.loc[(merged_na["Track_Mask Path"].isna()) & ~(merged_na["Unique_pos"].isin(org_Mask_missing_pos))].reset_index(drop = False)
#Make sure to remove the positions which have not had segmentation completed
Track_Mask_missing_pos = Track_Mask_missing["Unique_pos"].unique()

###########################################################################################################################################################################################################################
#%%
l_missing_orgMask = len(org_Mask_missing["Unique_pos"].unique()) # get the lengths of the missing at each stage
l_missing_TrackMask = len(Track_Mask_missing)

ret_fir = f"There are {l_missing_orgMask} missing segmentation masks, and {l_missing_TrackMask} tracking masks"
print(ret_fir)
with open('Missing_out', 'a+') as posit:
	posit.write(ret_fir)
posit.close()

time.sleep(10) # Pause for n seconds to allow for the number of postions which must be calculatd to be read
# seg = bool(input("Would you like to run the segmentation for the missing files? [True or False]"))
# Track = bool(input("Would you like to run the tracking for the missing files?"))


#%%

if seg == True:
    print("Running the segmenation")
    instances = org_Mask_missing["Unique_frame"].unique().tolist()
    cwd = os.getcwd()
    print (f"Placing results in {cwd}")
    uPos_list = imgIndex["Unique_pos"].drop_duplicates()
    for p in uPos_list:
        try:os.mkdir(p)
        except FileExistsError: pass
    del uPos_list

    imgIndex = imgIndex.loc[imgIndex["Unique_frame"].isin(instances)]

    now = datetime.datetime.now()

    positions = imgIndex.reset_index()
    imgIndex_T = imgIndex.reset_index()[
        ["Path", "Unique_frame", "Unique_pos", "Channel", "PositionID", "Frame"]]
    positions = positions["Unique_pos"]
    positions = positions.drop_duplicates()

    positions
    imgIndex_T

    date_today = datetime.date.today()


    os.chdir(Run_segProgLib)
    warnings.filterwarnings('ignore')
    id = 1

    # First, initiate the minidom root
    root = minidom.Document()

    # The first child is CellXFiles, create it's element,
    # and then set the attributes, then append to the XML document
    CellXFiles_root = root.createElement('CellXFiles')
    CellXFiles_root.setAttribute('Creator', 'CellXGui')
    CellXFiles_root.setAttribute('timestamp', str(now))
    root.appendChild(CellXFiles_root)

    for p in positions:
        # print to console to keep track which position we are working with
        print(p)

        idstr = str(id)

        # Next child is timeSeries, child of CellXFiles
        # Set attributes, append to CellXFiles
        timeSeries = root.createElement('CellXTimeSeries')
        timeSeries.setAttribute('fluotypes', 'GFP_mKO_mKa')
        # timeSeries.setAttribute('fluotypes', 'No_fluor')#no fluor
        timeSeries.setAttribute('id', idstr)
        timeSeries.setAttribute('position', '')
        timeSeries.setAttribute('tracking', '1')
        CellXFiles_root.appendChild(timeSeries)

        id = id + 1
        # Next is directory location for the output data
        # child of time series. Set attributes, append to
        # time series

        strp = str(p)
        directory = root.createElement('CellXResultDir')
        dirT = microfluidic_results+ "/" + strp
        directoryText = root.createTextNode(dirT)
        directory.appendChild(directoryText)
        timeSeries.appendChild(directory)

        # Next is the file set. This is where the for-loops will
        # probably come into play, since it looks for EACH time point
        # and this file set is a child of the timeSeries element

        pos_subset = imgIndex_T[imgIndex_T["Unique_pos"] == p]

        timepoint = np.unique(pd.to_numeric(pos_subset["Frame"], errors='coerce'))

        for i in timepoint:
            # The following will subset the dataframe to include only the
            # image paths for each channel, at the specified time point in the loopo
            pos_subset['Frame'] = pd.to_numeric(pos_subset['Frame'], errors='coerce')
            timeSubset = pos_subset[pos_subset['Frame'] == i]


            try: subImagePath = str(timeSubset[timeSubset['Channel'] == 'Sub'].iloc[0,0])
            except IndexError:
                try: subImagePath = str(timeSubset[timeSubset['Channel'] == 'out'].iloc[0,0])
                except IndexError:
                    continue ## If a frame does not exist, then exit and do not create a element

            fileSet = root.createElement('CellXFileSet')
            fileSet.setAttribute('frame', str(i))
            timeSeries.appendChild(fileSet)

            # Next up is the individual channels to be appended to the fileset
            # elements, this includes all 4 channels
            brightfield = root.createElement('oofImage')

            brightfieldText = root.createTextNode(subImagePath)
            brightfield.appendChild(brightfieldText)
            fileSet.appendChild(brightfield)
            if fluor == "True":
                channels = ['GFP','mKO', 'mKa']

                for j in channels:
                    try: channelImagePath = timeSubset[timeSubset['Channel'] == j].iloc[0,0]
                    except IndexError:
                        continue

                    fluorophoreSet = root.createElement('fluoSet')
                    fluorophoreSet.setAttribute('type', j)
                    fileSet.appendChild(fluorophoreSet)


                    fluorophoreImage = root.createElement('fluoImage')
                    fluorophoreText = root.createTextNode(channelImagePath)
                    fluorophoreImage.appendChild(fluorophoreText)
                    fluorophoreSet.appendChild(fluorophoreImage)
            else:
                pass


    xml_str = root.toprettyxml()
    save_path_file = "Target_" + str(date_today) + file_mod + ".xml"

    with open(save_path_file, "w") as f:
        f.write(xml_str)

    # r = str(os.cpu_count() - 1)

    # missing_imgIndex = pd.read_csv('imgIndex.csv')
    # missing_imgIndex = missing_imgIndex[(missing_imgIndex["Unique_pos"].isin(subset))]
    # missing_Pos_frame_list = missing_imgIndex[["Unique_pos", "Unique_frame"]].drop_duplicates().set_index(["Unique_pos"])

#%%
#%%
os.chdir(microfluidics_results)
orgAllmasks = pd.read_csv("orgAllmasks.csv")
seg_dirDF = orgAllmasks[["Unique_pos" , "Path"]].set_index(["Unique_pos"])


# Pos_frame_list = pd.concat(Pos_frame_list, missing_Pos_frame_list)

missing_Pos_frame_list = Track_Mask_missing["Unique_pos"].unique()
missing_imgIndex = imgIndex
Pos_frame_list = missing_Pos_frame_list
del missing_Pos_frame_list

raw_dirDF = imgIndex[imgIndex["Channel"] == "Sub"][["Unique_pos" , "Path"]].set_index(["Unique_pos"])
# Missing_raw_dirDF = missing_imgIndex[(missing_imgIndex["Channel"] == "Sub") | (missing_imgIndex["Channel"] == "out")][["Unique_pos" , "Path", "Channel"]].set_index(["Unique_pos"])
# del imgIndex
# del missing_imgIndex

# raw_dirDF = pd.concat(raw_dirDF, Missing_raw_dirDF)
# raw_dirDF = Missing_raw_dirDF
# del Missing_raw_dirDF


# seg_dirDF = seg_dirDF
# del seg_dirDF
# del missing_Allmasks
# del Allmasks

pos_index =  raw_dirDF.reset_index()["Unique_pos"].unique()

os.chdir(current_code)
#need to add auto BF_out vs sub function

def tracking_all_para(p):
    try:
        eng = matlab.engine.start_matlab()
        pos = p
        posFingerprint = "position" + p[p.find("p")+1:]
        length = len(Pos_frame_list.loc[p,])
        try: g = seg_dirDF.loc[p, ["Path"]].iloc[0]
        except KeyError: print(f"Failure at {p}")
        else:
            try: g = seg_dirDF.loc[p, ["Path"]].iloc[0].values[0]
            except AttributeError:
                return(f"Failure on {p}")
            # print(p)
            path_seg = os.path.dirname(g)
            path_raw = os.path.dirname(raw_dirDF.loc[p, ["Path"]].iloc[0].values[0])
            BFchan = raw_dirDF.loc[p, ["Channel"]].iloc[0].values[0]
            print(pos, current_date, posFingerprint, length, path_seg, path_raw, BFchan)
            if BFchan == 'Sub':
                eng.TrackALL(pos, current_date, posFingerprint, length, path_seg, path_raw, nargout=0)
            elif BFchan == 'out':
                eng.TrackALL_out(pos, current_date, posFingerprint, length, path_seg, path_raw, nargout=0)
    except:
        return(f"Error on {p}")
    eng.quit()


pn = os.cpu_count()
l = len(pos_index)

if l >= pn:
    pr = pn-1
else:
    pr = l # this is an L not a 1!
#This is just to make sure that not too many cores are reserved
for p in pos_index:
    tracking_all_para(p)




# Parallel(n_jobs=pr, verbose = 100)(delayed(tracking_all_para)(p) for p in pos_index)

#%%

    #need to add auto BF_out vs sub function

    # def tracking_all_para(p):
    #     try:
    #         eng = matlab.engine.start_matlab()
    #         pos = p
    #         posFingerprint = "position" + p[p.find("p")+1:]
    #         length = len(Pos_frame_list.loc[p,])
    #         try: g = seg_dirDF.loc[p, ["Path"]].iloc[0]
    #         except KeyError: print(f"Failure at {p}")
    #         else:
    #             try: g = seg_dirDF.loc[p, ["Path"]].iloc[0].values[0]
    #             except AttributeError:
    #                 return(f"Failure on {p}")
    #             # print(p)
    #             path_seg = os.path.dirname(g)
    #             path_raw = os.path.dirname(raw_dirDF.loc[p, ["Path"]].iloc[0].values[0])
    #             BFchan = raw_dirDF.loc[p, ["Channel"]].iloc[0].values[0]
    #             print(pos, current_date, posFingerprint, length, path_seg, path_raw, BFchan)
    #             if BFchan == 'Sub':
    #                 eng.TrackALL(pos, current_date, posFingerprint, length, path_seg, path_raw, nargout=0)
    #             elif BFchan == 'out':
    #                 eng.TrackALL_out(pos, current_date, posFingerprint, length, path_seg, path_raw, nargout=0)
    #         eng.quit()
    #     except:
    #         return(f"Error on {p}")

    #This is just to make sure that not too many cores are reserved


    # Parallel(n_jobs=pr, verbose = 100)(delayed(tracking_all_para)(p) for p in pos_index)
# #%%

# print("Re-testing... ")
# os.chdir(microfluidics_results)

# orgmaskpaths = []
# count = 0
# fluors = "GFP_mKO_mKa"
# exclude = list(['No_fluor', 'No_flour', 'mKO_mKa'])
# for root, dirs, files, in os.walk(os.getcwd(), topdown= True):
#     [dirs.remove(d) for d in list(dirs) if d in exclude]
#     for name in files:
#         if name.endswith(".mat") and name.startswith("mask"):
#             orgmaskpaths.append({'Path': os.path.join(root, name)})
#             count = count + 1
#             print(count, end="\r")
#         else:
#             pass

# def f_Pos_mask(z):
#     start = z.find('p') #Note: The shift forward is dependant upon how you write out position
#     end = z.find('GFP')-1  #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
#     pos = z[start:end]
#     # if len(pos) <= 8:
#     #     return(pos)
#     # else:
#     #     pass
#     return(pos)


# def mask_inf(m):
#     num1 = m.find('ask_')+5
#     num2 = m.find('.mat')
#     num_ = ("f"+ m[num1:num2])
#     return(num_)

# def mask_hour(m):
#     hourS = m.find("_")+2
#     hourE = m.find("_p")
#     return(m[hourS:hourE])

# def mask_barcode(m):
#     start =m.find('d02')
#     end = m.find('GFP') -1
#     return (m[start:end])

# def f_maskdate(x):
#     dstart = x.find('d0')
#     dend = dstart + 5
#     expdate = (x[dstart:dend])
#     return(expdate)

# orgAllmasks = pd.DataFrame(orgmaskpaths)

# def f_non_sync(p):
#     substring = ".sync"
#     if substring in p:
#         return(None)
#     else:
#         return(p)
# orgAllmasks["Path"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_non_sync)
# orgAllmasks.dropna(inplace=True)

# orgAllmasks["Date"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_maskdate)
# orgAllmasks["Position_ID"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_Pos_mask)
# orgAllmasks["MaskNum"] = pd.Series(orgAllmasks.iloc[:,0]).apply(mask_inf)
# # Allmasks["Hour/Run_t"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_hour)
# # Allmasks["Unique_pos_hour"] = Allmasks['Date'] + "h" + Allmasks["Hour/Run_t"] + Allmasks['Position_ID']

# orgAllmasks["Unique_pos"] = pd.Series(orgAllmasks.iloc[:,0]).apply(mask_barcode)
# orgAllmasks["Unique_frame"] = orgAllmasks["Unique_pos"] + orgAllmasks['MaskNum'] #By performing this way, can avoid one extra step per index

# del orgmaskpaths

# orgAllmasks.dropna(inplace=True) # temp drop because there are some extra segmentations
# orgAllmasks.reset_index(inplace=True, drop = True)

# orgAllmasks.set_index(["Unique_frame"], inplace = True)
# orgAllmasks

# maskpaths = []
# count = 0
# for root, dirs, files, in os.walk(os.getcwd()):
#     for name in files:
#         if name.endswith(".mat") and name.startswith("track_mask"):
#             maskpaths.append({'Path': os.path.join(root, name)})
#             count = count + 1
#             print(count, end="\r")
#         else:
#             pass

# # def f_Pos_mask(z):
# #     start = z.find('RESULTS') + 8 #Note: The shift forward is dependant upon how you write out position
# #     end = z.find('GFP')-1  #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
# #     pos = z[start:end]
# #     if len(pos) <= 16:
# #         return(pos)
# #     else:
# #         pass

# def mask_inf(m):
#     num1 = m.find('ask_')+5
#     num2 = m.find('.mat')
#     num_ = ("f"+ m[num1:num2])
#     return(num_)

# def mask_hour(m):
#     hourS = m.find("_")+2
#     hourE = m.find("_p")
#     return(m[hourS:hourE])

# def mask_barcode(m):
#     start =m.find('d02')
#     end = m.find('GFP') -1
#     return (m[start:end])

# def f_maskdate(x):
#     dstart = x.find('d0')
#     dend = dstart + 5
#     expdate = (x[dstart:dend])
#     return(expdate)

# Allmasks = pd.DataFrame(maskpaths)

# def f_non_sync(p):
#     substring = ".sync"
#     if substring in p:
#         return(None)
#     else:
#         return(p)
# Allmasks["Path"] = pd.Series(Allmasks.iloc[:,0]).apply(f_non_sync)
# Allmasks.dropna(inplace=True)

# Allmasks["Date"] = pd.Series(Allmasks.iloc[:,0]).apply(f_maskdate)
# # Allmasks["Position_ID"] = pd.Series(Allmasks.iloc[:,0]).apply(f_Pos_mask)
# Allmasks["MaskNum"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_inf)
# # Allmasks["Hour/Run_t"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_hour)
# # Allmasks["Unique_pos_hour"] = Allmasks['Date'] + "h" + Allmasks["Hour/Run_t"] + Allmasks['Position_ID']

# Allmasks["Unique_pos"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_barcode)
# Allmasks["Unique_frame"] = Allmasks["Unique_pos"] + Allmasks['MaskNum'] #By performing this way, can avoid one extra step per index

# # del maskpaths

# # Allmasks.dropna(inplace=True) # temp drop because there are some extra segmentations
# Allmasks.reset_index(inplace=True, drop = True)

# Allmasks.set_index(["Unique_frame"], inplace = True)
# Allmasks

# orgAllmasks.to_csv('orgAllmasks.csv')
# Allmasks.to_csv('Allmasks.csv')

# # imgIndex = pd.read_csv('imgIndex.csv')[["Unique_frame", "Path"]].rename(columns={"Path" : "image_Path"})
# # frame_list = imgIndex["Unique_frame"].unique()

# orgAllmasks = pd.read_csv("orgAllmasks.csv")[["Unique_frame", "Path"]].rename(columns={"Path": "org_Mask Path"})

# Allmasks = pd.read_csv('Allmasks.csv')[["Unique_frame", "Path"]].rename(columns={"Path": "Track_Mask Path"})
# Allmasks = Allmasks.set_index(["Unique_frame"])

# merged = pd.merge(imgIndex, orgAllmasks, left_index= True, right_index= True, how = 'outer') # This is by default a left merge

# # if M_instances == "All"
# # 	pass

# # else:
# merged = merged.loc[merged["Unique_pos"].isin(M_instances)] # Substet to just the files that are immediately needed

# merged = pd.merge(merged, Allmasks, left_index= True, right_index= True, how = 'left')
# merged_na = merged.loc[(merged["org_Mask Path"]).isna() | (merged["Track_Mask Path"].isna())] # If there is a missing mask (either from seg or tracking), it will be labelled here

# merged_na.to_csv("Files_STILL_missing_df.csv") #output the files that are missing

# org_Mask_missing = merged_na.loc[merged_na["org_Mask Path"].isna()].reset_index(drop = False)
# imgIndex = merged_na.loc[merged_na["Track_Mask Path"].isna()].reset_index(drop = False)

# l_missing_orgMask = len(org_Mask_missing) # get the lengths of the missing at each stage
# l_missing_TrackMask = len(imgIndex)


# ret_sec =f"STILL MISSING ----- There are {l_missing_orgMask} missing segmentation masks, and {l_missing_TrackMask} tracking masks"
# print(ret_sec)

# with open('Missing_out', 'a+') as posit:
# 	posit.write(ret_sec)
# posit.close()




# def mask_load(m):
#     mask = loadmat(m)
#     a = mask['data']
#     return (a)

# f = "d0221r1p500200f0002"
# info_concat_s = info_concat[info_concat["cell_frame"] == 2]
# info_concat_s = info_concat_s[info_concat_s["Budneck"] == "mKa"]
# info_concat_s.reset_index(inplace=True, drop = True)


# mask = pd.DataFrame(mask_load(Allmasks.loc[f,["Path"]].values[0]))

# with open('Coupling_errors.txt', 'a+') as couple:
# 	for f in frame_list:
# 		for i in range(len(info_concat_s)):
# 			# f = info_concat.iloc[i,0]
# 			vt = info_concat_s.iloc[i,1]
# 			x = info_concat_s.iloc[i,2]
# 			y = info_concat_s.iloc[i,3]

# 			vm = mask.iloc[y,
# 			x]
# 			if vt == vm:
# 				continue
# 			elif vt != vm:
# 				print(f"Error at frame {f} with info {vt} --> mask {vm}")
# 				couple.write(f"Error at frame {f} with info {vt} --> mask {vm}\n")
# 			else:
# 				print("ERROR?")
# 	couple.close()


# count = len(open("Coupling_errors.txt").readlines())
# if count > 0:
# 	# print(f"There are {count} uncoupling errors. Would you like to continue?")
# 	con = input(f"There are {count} uncoupling errors. Would you like to continue? <Y/N>")
# 	if con == "Y":
# 		print("Continuing with analysis")
# 		cont = True
# 	elif con == "N":
# 		print("Exiting analysis")
# 		cont = False
# else:
# 	cont = True

# %%
