#%%
import csv
from operator import pos
from xml.etree.ElementTree import Element, SubElement, Comment, TreeBuilder, tostring
import datetime
from xml.dom import minidom
import os
import pandas as pd
import numpy as np
# import datetime
import warnings
from pandas.core.indexing import IndexingError
#%%
# subset_by = 'Position'
# subset_first = ['d0214r1p110300']
# # subset_second = ['d0214r2p150300', 'd0214r2p160200', 'd0214r2p160300']
#%%
microfluidics_results = 'D:/Microfluidics/RESULTS'
# microfluidics_results = 'F:/Microfluidics/RESULTS'
os.chdir(microfluidics_results)

# %%
# The following generates the basis of the XML structure
# which will be populated with the metadata from the image
# file name
# os.chdir(microfluidics_results)
imgIndex = pd.read_parquet('imgIndex.parquet')
#// Subset for filling May23
imgIndex = imgIndex.astype({'Frame':'int', 'Run': 'int'})
# imgIndex = imgIndex.loc[(imgIndex['Date'] == 'd0214') & (imgIndex['Run'] == 1) & (imgIndex['Frame'] <= 45) & (imgIndex['Unique_pos'].isin(subset_first))] #* Fill specific frames for specific positions for d0214r1 and fill all frames for specific pos d0214r2
# imgIndex = imgIndex.loc[(imgIndex['Unique_pos'].isin(subset_first))] #* Fill specific frames for specific positions for d0214r1 and fill all frames for specific pos d0214r2


#%%
# THIS SHOULD ALREADY BE DONE IN THE PREVIOUS PIPELINE
print(f"Placing results in {os.getcwd()}")
uPos_list = imgIndex["Unique_pos"].drop_duplicates()
for p in uPos_list:
    try:os.mkdir(p)
    except FileExistsError: pass
del uPos_list
#%%

now = datetime.datetime.now()
date_today = datetime.date.today()
# if Experimental_info["data_subset"] == True:
#     imgIndex = imgIndex[imgIndex["Unique_pos"].isin(instances)]
# else:
#     pass

positions = imgIndex.reset_index()
imgIndex_T = imgIndex.reset_index()[
    ["Path", "Unique_frame", "Unique_pos", "Channel", "PositionID", "Frame"]]
positions = positions["Unique_pos"]
positions = positions.drop_duplicates()

# positions
# imgIndex_T

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
    dirT = microfluidics_results+ "/" + strp
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

        channels = ['GFP','mKO', 'mKa']

        # for j in channels: #* Do not run the quantification for Fluor cannels
        #     try: channelImagePath = timeSubset[timeSubset['Channel'] == j].iloc[0,0]
        #     except IndexError:
        #         continue

        #     fluorophoreSet = root.createElement('fluoSet')
        #     fluorophoreSet.setAttribute('type', j)
        #     fileSet.appendChild(fluorophoreSet)


        #     fluorophoreImage = root.createElement('fluoImage')
        #     fluorophoreText = root.createTextNode(channelImagePath)
        #     fluorophoreImage.appendChild(fluorophoreText)
        #     fluorophoreSet.appendChild(fluorophoreImage)

xml_str = root.toprettyxml()
save_path_file = "Target_" + str(date_today) + ".xml"

with open(save_path_file, "w") as f:
    f.write(xml_str)
