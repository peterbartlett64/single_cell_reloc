#%%
from PIL import Image
import pandas as pd
import numpy as np
from numpy import asarray
import os
from scipy.io import loadmat
import math
from icecream import ic
from joblib import Parallel, delayed

#%%
# Load an image using OpenCV
def read_and_modify_intensity(image_path:str, cell_mask_path:str, quantification_path:str, unique_frame:str,normalization_method:str):
	"""This function will take in the image path and cell_mask to properly visualize cell measurements because of normalization methods. Function was written to be flexible to normalization method so that is can be used regardless of whether it was the background or the cell itself that the brightest pixel values done on. This will be modifying the image so the percentile on which calculations are made is blinded for the computer and user.
		Args:
			image_path (path): path of the image file to be loaded and corrected
			cell_mask_path (path): path of the mask file to be designate which cell is which
			quantification_path (path): path of the quantification file for that protein/frame within it so that the cell factors can be extracted
			normalization_method (str): column designating what was used as denominator in brightest pixel value normalization. Pixel values for idividual cells will be scaled by the associated factor for the cell-time
	"""

	#* Read in the image to be manipulated
	# image = np.array(Image.open(image_path))
	image_raw = Image.open(image_path)
	image = asarray(image_raw)
	image = image.astype('float32')
	if image is None:
		print(f"Error: Image {image_path} not found or could not be loaded.")
	else:
		print(f'Loaded {image_path} correctly and continuing')

	#* Get the paths required as variables from the lookup table
	try:
		quant_file = pd.read_parquet(quantification_path)
	except:
		quant_file = pd.read_csv(quantification_path)

	cell_mask  = loadmat(cell_mask_path)
	# cell_mask = cell_mask['data'] #! This is for regular. Testing with expanded
	cell_mask = cell_mask["px_data"]
	# cell_factor = np.copy(cell_mask).astype('float32')
	cell_factor = np.empty(cell_mask.shape).astype('float32')

	#* make a new file name for saving
	# filename = os.path.splitext(os.path.basename(image_path))[0] + "_cellNormalized.tif" #* Had to change from this becasue of shared position names between days. Using the barcoded version below
	# filename = unique_frame + "_cellNormalized.tif" #! This is for regualar. Testing with expanded
	filename = unique_frame + "_Expanded_cellNormalized.tif"

	#* Make a copy of the quant_file so that it does not get overwritten
	quant_copy = quant_file.copy()
	quant_file.set_index('Cell_Barcode', inplace = True)
	cells = pd.DataFrame(quant_copy['Cell_Barcode'].drop_duplicates())

	cells['N'] = pd.Series(cells['Cell_Barcode']).apply(lambda x: x[x.find('c')+1:]).astype('uint16')
	for i in range(len(cells)): #! This for generating the mask, but it does not do the multiplication yet
		cell_n = cells.iloc[i,:]['N']
		cell_bar = cells.iloc[i,:]['Cell_Barcode']
		factor = (quant_file.loc[cell_bar, normalization_method])
		factor_float = 1/factor
		np.putmask(cell_factor,(cell_mask == cell_n) & (cell_factor<1.2), factor_float)

	image_mod_arr = image * cell_factor# * 255 #* This is the multiplication of the file. This could be sped up by using a tensor in a future version
	image_mod = Image.fromarray(image_mod_arr)
	image_mod.save(filename)
	# cv2.imwrite(filename, image_mod)

	#, The below is an in progrss version using numba. Might be better to just write for the GPU
	#< cell_factor = np.empty(cell_mask.shape).astype('float32')

	#< @numba.njit
	#< def f(cell_factor, ref, factor, factor_float):
	#< 	cell_factor = np.empty(ref.shape)
	#< 	for i in range(ref.shape[0]):
	#< 		for j in range(ref.shape[1]):
	#< 			if ref[i, j] == factor:
	#< 				cell_factor[i,j] = factor_float
	#< 			else:
	#< 				pass
	#< 	return cell_factor

	#. Need to save the new image
	return(True)
#%%
def reloc_mask_gen(cell_mask_path: str, reloc_df_small, unique_frame:str):
	"""This function will take in t2023-10-28e image path and cell_mask to properly visualize cell measurements because of normalization methods. Function was written to be flexible to normalization method so that is can be used regardless of whether it was the background or the cell itself that the brightest pixel values done on. This will be modifying the image so the percentile on which calculations are made is blinded for the computer and user.
		Args:
			cell_mask_path (path): path of the mask file to be designate which cell is which
			reloc_df_small (pandas DataFrame): dataframe with the information on which cells are currently localized
			unique_frame(str): the frame barcode
	"""

	#* Get the paths required as variables from the lookup table
	cell_mask  = loadmat(cell_mask_path)
	# cell_mask = cell_mask['data']
	cell_mask = cell_mask["px_data"]
	cell_states = np.copy(cell_mask) #* This one can be left as int

	#* make a new file name for saving
	filename = os.path.splitext(os.path.basename(cell_mask_path))[0] + "_state.tif" #. Check that this line will not give '.tiff_cellNormalized.tiff'

	#* Make a copy of the quant_file so that it does not get overwritten
	reloc_copy = reloc_df_small.copy().reset_index(drop = False)
	reloc_df_small.set_index('Cell_Barcode', inplace = True)
	cells = pd.DataFrame(reloc_copy['Cell_Barcode'].drop_duplicates())

	cells['N'] = pd.Series(cells['Cell_Barcode']).apply(lambda x: x[x.find('c')+1:]).astype('uint16')

	#* The below is required to deal with the 100 and 200 cell before this gets crowded. This wasn't needed because referencing one array to make changes on the other
	#< temp = cells['Cell_Barcode'][0]
	#< prefix = (lambda x : x[:x.find("c")+1])(temp)
	#< prefix = prefix + "0"

	#< do_first = [100, 200]
	#< for d in do_first:
	#< 	cell_bar = prefix + d
	#< 	try:
	#< 		cell_n = d
	#< 		state = (reloc_df_small.loc[cell_bar, 'Relocalized']) + 1
	#< 		np.putmask(cell_states, (cell_mask == cell_n), state)
	#< 	except:
	#< 		continue

	for i in range(len(cells)): #! This for generating the mask, but it does not do the multiplication yet
		cell_n = cells.iloc[i,:]['N']
		cell_bar = cells.iloc[i,:]['Cell_Barcode']
		state = (reloc_df_small.loc[cell_bar, 'Relocalized']) + 1 #! The 1 is added so that the background is 0, the unlocalized is 1 =>100 and relocalized is 2 => 200
		np.putmask(cell_states, (cell_mask == cell_n) &(cell_states<2), state) # Check that the pixels have not already be labeled as relocalized

	image_mod_arr = cell_states * 100 #* This is the multiplication of the file. This could be sped up by using a tensor in a future version
	image_mod = Image.fromarray(image_mod_arr)
	image_mod.save(filename)
	# cv2.imwrite(filename, image_mod)

	#. Need to save the new image
	return(True)
#%%
def correction_manager(p_i: int, normalization_method:str): #* This is the manager to call the paths for use
	"""
	This function takes in an index for the merged table for image, track_masks, and quantification paths.

	Args:
		p_i: index of pad to be loaded
		normalization_method (str): column designating what was used as denominator in brightest pixel value normalization. Pixel values for idividual cells will be scaled by the associated factor for the cell-time
	"""
	z = read_and_modify_intensity(image_path = merged_indices.iloc[p_i,:]["Path_image"], cell_mask_path = merged_indices.iloc[p_i, :]["Path_mask"], quantification_path= merged_indices.iloc[p_i,:]['Path_quant'], unique_frame = merged_indices.iloc[p_i,:]['Unique_frame'],  normalization_method = normalization_method)
	ic(z)
	return(z)

def move_mask_manager(frame:str): #* This is the manager to call the paths for use
	"""
	This function takes in an index for the merged table for image, track_masks, and quantification paths.

	Args:
		p_i: index of pad to be loaded
		normalization_method (str): column designating what was used as denominator in brightest pixel value normalization. Pixel values for idividual cells will be scaled by the associated factor for the cell-time
	"""
	z = reloc_mask_gen(cell_mask_path=mask_results_merge.loc[frame, "Path"][0], reloc_df_small = mask_results_merge.loc[frame,:], unique_frame=frame)
	return(z)


#%%
if __name__ == '__main__':
	# Global_Variables = gv.global_manager()
	Global_variables = {'analyze': 'E:/Microfluidics/Analyze',
	'microfluidics_results': 'E:/Microfluidics/RESULTS',
	'post_path': 'E:/ALL_FINAL', #. gv.slash_switch(input("Post quant path?")) , #Todo: This needs to be changed to a input call
	'subset': False,
	'subset_by': '',
	'subset_collection': '',
	'cpu_se': int(math.floor(os.cpu_count()*0.7)),
	'timepoint_gap': 7.5,
	'percentiles': [95, 99],
	'multiplex': True,
	'figures_root': 'D:/Figures_root',
	'image_mod_folder': 'D:/Manipulated'}
	#Todo: Update the global variables call to have a key for image_mod_folder

	os.chdir(Global_variables['image_mod_folder']) #* This was

	norm_default = "factor_median_OBJ_GFP"
	norm_confirm = input(f"Was the normalization by {norm_default}")
	if norm_confirm.lower() == 'y' or norm_confirm.lower() == 'yes':
		normalization_method = norm_default
	else:
		normalization_method = input("What was the normalization method?")

	#* Get the masks for combining
	imgIndex = pd.read_parquet("E:/Microfluidics/RESULTS/imgIndex.parquet").reset_index(drop = False)
	# mask_index = pd.read_parquet("E:/Microfluidics/RESULTS/Allmasks.parquet").reset_index(drop = False) #! This is the regualar mask. Testing out the expanded below
	mask_index = pd.read_parquet("E:\Microfluidics\RESULTS\Allmasks_exp.parquet").reset_index(drop = False)
	quant_index = pd.read_parquet("D:\ALL_FINAL\Quantification_index.parquet")[['Path', 'Frame']].rename(columns={'Path':'Path_quant'}).sort_values(by = "Frame") #* Read in the quant_index and keep only the relevant columns
	results_final = pd.read_parquet("E:\ALL_FINAL\Combined_by_perc\FLR1_selected.parquet", columns=['Cell_Barcode', 'ImageID', 'Relocalized', 'Loc_score', 'log_Loc_score']).reset_index(drop = False) #* This will take a milwhile since it is a 4G file


	mask_results_merge = pd.merge(results_final, mask_index, left_on='ImageID', right_on='Unique_frame',suffixes = ('', '_mask'))

	frames = mask_results_merge['Unique_frame']
	mask_results_merge.set_index('Unique_frame', inplace = True)


	# quant_index = pd.read_parquet("D:/MOST FINAL- Second/Quant_ALL_index.parquet")[['Path', 'PositionID']].rename(columns={'Path':'Path_quant'})
	merged_indices = pd.merge(imgIndex, mask_index, left_on = 'Unique_frame', right_on = "Unique_frame",  suffixes = ('_image', '_mask'))[['Unique_frame', 'Path_mask', 'Path_image', 'Channel', 'Unique_pos_image']]
	merged_indices = merged_indices.loc[(merged_indices['Channel'] == 'GFP')| (merged_indices['Channel'] == 'gre')].dropna().reset_index(drop = True) #* Keep only the paths where cells are GFP
	merged_indices = pd.merge(merged_indices, quant_index, left_on = 'Unique_frame', right_on = 'Frame')
	# merged_indices = pd.merge(merged_indices, quant_index, left_on = 'Unique_pos_image', right_on = 'PositionID')

	# merged_indices = pd.read_parquet('merged_indices,parquet') #* This has not been created yet but is very simple to do so
	copy = merged_indices.copy()
	merged_indices = merged_indices.loc[(merged_indices['Unique_pos_image'] == "d0222r2p240200")]# | (merged_indices['Unique_pos_image'] == "d0222r2p300200")] #. These are barcodes for FLR1. Don't forget to change if looking for a different protien
	# y = Parallel(n_jobs=Global_variables['cpu_se'], verbose= 100)(delayed(correction_manager)(p_i = p_i, normalization_method= normalization_method) for p_i in range(len(merged_indices)))

	for p_i in range(len(merged_indices)):
		ic(correction_manager(p_i = p_i, normalization_method= normalization_method))

	for f in frames:
		ic(move_mask_manager(frame = f))

# %%
