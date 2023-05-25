#%%
import os
import pandas as pd
import numpy as np
from pandas.core import frame
from pandas.core.arrays.sparse import dtype
from scipy.io import loadmat
from scipy.io.matlab.mio import savemat

#%%

#%%
# Temp = "C:/Users/Peter/Dropbox (Grant Brown's Lab)/Peter Bartlett Data/Code/Other versions_and_scripts"
# microfluidics_results =  "C:/Users/Peter/Desktop/Subset_Images/DEC15RESULTS"

def switch_slash(path):
	new = path.replace(os.sep, '/')
	return (new)

microfluidics_results = str(input("microfluidics_results"))
microfluidics_results = switch_slash(microfluidics_results)

os.chdir(microfluidics_results)

#%%
run_Allmasks = True
subset = "ALL"

#%%
if run_Allmasks == True:
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
	Allmasks.sort_values(by = "Unique_frame", inplace= True)
	Allmasks.set_index(["Unique_frame"], inplace = True)
	Allmasks.to_csv('Allmasks.csv')
	Allmasks.to_parquet('Allmasks.parquet')
else:
	Allmasks = pd.read_csv("Allmasks.csv", index_col = 0)

if subset == "ALL":
	pass
else:
	subset = subset.split(", ")
	Allmasks = Allmasks.loc[Allmasks["Unique_pos"].isin(subset)]

#%%
#%%


# Allmasks = pd.read_csv("Allmasks.csv", index_col = 0)

# def mask_load(m):
# 	mask = loadmat(m)
# 	data = mask['data']
# 	return (data)

# for r in range(len(Allmasks)):
# 	path = Allmasks.iloc[r, 0]
# 	f = path[path.find("ask_")+4:path.find(".mat")] # this could be changed to target the frame column of Allmasks just after the "f'"
# 	frame_pxMAT = mask_load(path)
# 	n = np.max(frame_pxMAT)
# 	copy = frame_pxMAT
# 	for c in range(1, n+1):
# 		smaller = (frame_pxMAT== c).nonzero()
# 		a,b = smaller
# 		l = len(a)
# 		for pix in range(l):
# 			copy[(a[pix]-4):(a[pix]+5), (b[pix]-4):(b[pix]+5)] = c
# 	save_loc = os.path.dirname(path)
# 	os.chdir(save_loc)
# 	savemat(f"frame_pxMAT_{f}.mat", {"px_data": copy})##The name of the matrix is now called "px_data" instead of "data" . This was done to stop accident misreads in.
# 	print(f"frame_pxMAT_{f}.mat")

# try: a
# except NameError: pass
# else: del a

# try: b
# except NameError: pass
# else: del b
#%%

from joblib import Parallel, delayed

def mask_load(m):
	mask = loadmat(m)
	data = mask['data']
	return (data)
#%%
def mask_expand(r):
	path = Allmasks.iloc[r, 0]
	f = path[path.find("ask_")+4:path.find(".mat")] # this could be changed to target the frame column of Allmasks just after the "f'"
	frame_pxMAT = mask_load(path)
	n = np.max(frame_pxMAT)
	copy = frame_pxMAT
	for c in range(1, n+1):
		smaller = (frame_pxMAT== c).nonzero()
		a,b = smaller
		l = len(a)
		for pix in range(l):
			copy[(a[pix]-4):(a[pix]+5), (b[pix]-4):(b[pix]+5)] = c
	save_loc = os.path.dirname(path)
	os.chdir(save_loc)
	savemat(f"frame_pxMAT_{f}.mat", {"px_data": copy})##The name of the matrix is now called "px_data" instead of "data" . This was done to stop accident misreads in.


pn = os.cpu_count()
pn = int(input(f"Machine has {pn} cores. What is the max cores to use? [if all, press ESC]"))

print(f"function_starting with {pn} cores")

Parallel(n_jobs = pn, verbose = 100)(delayed(mask_expand)(r) for r in range(len(Allmasks)))

# # %%
# n_f = "00063"
# savemat(f"frame_pxMAT_{n_f}.mat", copy)

# #%%
# test_load = loadmat("frame_pxMAT_00001.mat")['data']
#%%
# mask_folder = "C:/Users/Peter/Desktop/DATA_TEST/20180210/position740200_CellXOutput/GFP"
# os.chdir(mask_folder)

# TempTest = str(f"{microfluidics_results}/1TemTest")

# os.chdir(TempTest)
os.chdir(microfluidics_results)

maskpaths = []
count = 0
for root, dirs, files, in os.walk(os.getcwd()):
	for name in files:
		if name.endswith(".mat") and name.startswith("frame_pxMAT"):
			maskpaths.append({'Path': os.path.join(root, name)})
			count = count + 1
			print(count, end="\r")
		else:
			pass

def f_Pos_mask(z):
	start = z.find('RESULTS') + 10 #Note: The shift forward is dependant upon how you write out position
	end = z.find('GFP_')-1  #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
	pos = z[start:end]
	if len(pos) <= 17:
		return(pos)
	else:
		pass

def mask_num(m):
	num1 = m.find('pxMAT')+7
	num2 = m.find('.mat')
	num_ = ("f"+ m[num1:num2])
	return(num_)

Allmasks_exp = pd.DataFrame(maskpaths)
# Allmasks_exp["Date"] = pd.Series(Allmasks_exp.iloc[:,0]).apply(f_maskdate)
# Allmasks_exp["Position_ID"] = pd.Series(Allmasks_exp.iloc[:,0]).apply(f_Pos_mask)
Allmasks_exp["MaskNum"] = pd.Series(Allmasks_exp.iloc[:,0]).apply(mask_num)
# Allmasks_exp["Hour/Run_t"] = pd.Series(Allmasks_exp.iloc[:,0]).apply(mask_hour)
# Allmasks_exp["Unique_pos_hour"] = Allmasks_exp['Date'] + "h" + Allmasks_exp["Hour/Run_t"] + Allmasks_exp['Position_ID']

Allmasks_exp["Unique_pos"] = pd.Series(Allmasks_exp.iloc[:,0]).apply(f_Pos_mask)
Allmasks_exp["Unique_frame"] = Allmasks_exp["Unique_pos"] + Allmasks_exp['MaskNum'] #By performing this way, can avoid one extra step per index

Allmasks_exp.set_index(["Unique_frame"], inplace = True)

def f_maskdate(x):
	expdate = (x[0:5])
	return(expdate)
Allmasks_exp["Date"] = pd.Series(Allmasks_exp["Unique_pos"]).apply(f_maskdate)
# Allmasks_exp = Allmasks_exp[Allmasks_exp["Date"] == 'd0224']

Allmasks_exp.to_csv('Allmasks_exp.csv')
Allmasks_exp.to_parquet('Allmasks_exp.csv')
#%%
# from PIL import Image

# def mask_load(m):
#     mask = loadmat(m)
#     data = mask['data']
#     return (data)

# Allmasks = pd.read_csv('Allmasks.csv')
# Allmasks_exp = pd.read_csv('Allmasks_exp.csv')
# #%%
# r = 344
# path = Allmasks.iloc[r, 1]
# print(path)
# frame_pxMAT = mask_load(path)
# frame_pxMAT[frame_pxMAT >=1 ] = 1

# # x = Image.fromarray(frame_pxMAT)
# # x
# reg = Image.fromarray(np.uint8(frame_pxMAT)*255)
# reg

# #%%
# def bud_mask_load(m):
# 	mask = loadmat(m)
# 	data = mask['px_data']
# 	return (data)


# path = Allmasks_exp.iloc[r, 1]
# print(path)
# frame_pxMAT = bud_mask_load(path)
# frame_pxMAT[frame_pxMAT >=1 ] = 1

# exp = Image.fromarray(np.uint8(frame_pxMAT)*255)
# exp




# from matplotlib import pyplot as plt
# plt.imshow(frame_pxMAT, interpolation='nearest')
# # plt.show()
