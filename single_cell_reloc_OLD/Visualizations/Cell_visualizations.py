#%%
from joblib import delayed, Parallel
from PIL import Image
import cv2
import pandas as pd
from scipy.io import loadmat
import numpy as np

#%%
def crop(img, matix,cell_n):
	a,b = (matix== cell_n).nonzero()
	x1 = np.min(a)
	x2 = np.max(a)
	y1 = np.min(b)
	y2 = np.max(b)
	img = img[x1:x2+1, y1:y2+1]
	return(img)

#%%
#, Mask loading function which will be moved to the gloabl functions group
def mask_load(m):
	mask = loadmat(m)
	data = mask['data']
	return (data)

#, Create a visual representation of the pixels that are higher or lower than the set loc_thresholds
#? This could be modified to accept the composite images and modifying just the GFP layer then reconsitute. This is incomplete but the composite branch is partially formed
#- barcode: can be in the form of either a frame barcode or a postitional barcode
#- composite: a boolean to designate whether individual or multiframe composites are being modified
#- myo_chan: designates which fluorescent channel the myo1 is being caught on. This would be flexible to demuliplexing tag
#- blackout: a boolean to designate whether the cells from the other strain should be blacked out
#- spec cell: accepts a cell barcode to crop for. If non-true, the function will generate a smaller image where just the one cell is visible. The size of this image will vary based on cell size which will also vary accross the timecourse
#? The above point may be worth correcting in a future form

def visual_images(barcode, myo_chan, composite = False, blackout = False, spec_cell = False) -> None:
	ReadInList = [f"Raw_factorUpper_{myo_chan}", f"Raw_factorLower_{myo_chan}", "Myo1Indentity"]

	temp_thresh_table = pd.read_csv(f"Movement_treat_course_{barcode}.csv", usecols= ReadInList) #* This can either be the composite
	row =  temp_thresh_table.iloc[0:1] #* leave a small dataframe instead of converting to a series
	threshold_higher = row["threshold_higher"]
	threshold_lower = row["threshold_lower"]

	if composite == False:
		image_path = row["Path"]
		image_frame = row["image_barcode"]

	if composite == True: #* right now the program is only written to accept a single image, and it would require previous and future work to do with the composite
		img = cv2.imread(image_path) #* keep the image as an image instead of an np.array()
		img = img[:,:,1] #* select the green layer. Ie. the GFP which is tagged on the proteins of interest
		image_position = row["Position_id"] #! Check to make sure that the this


	elif composite == False:
		img = cv2.imread(image_path) #* this would be for reading in the greyscale image

	else:
		return(f"{barcode}.Error exit") #* If composite is neither true or false. then the program should end instead of throwing a terminal error

	if blackout == True or spec_cell != False: #* In either case of blacking out or cell cropping, the image mask must be loaded
		seg_index = seg_index
		frame_pxMAT = mask_load(seg_index.loc[barcode,["Path"]].values[0])

		#< Loop for creation of croppped images
		if spec_cell != False:
			cell_n = int(spec_cell[spec_cell.find("c"):]) #* Convert the cell barcode to a cell number which can be found in the mask

			def crop(img_f, matrix_f,cell_n_f): #* This might not be needed in the new version of joblib. With cloudpickle, it might be available with global declaration
				a,b = (matrix_f == cell_n_f).nonzero() #* Ge the x and y locations where the mask is equal to the cell of interest
				x1 = np.min(a)
				x2 = np.max(a)
				y1 = np.min(b)
				y2 = np.max(b)
				img_r = img_f[x1:(x2+1), y1:(y2+1)]
				mask_r =matrix_f[x1:(x2+1), y1:(y2+1)]
				return(img_r, mask_r)
			img, frame_pxMAT = crop(img_f = img, matrix_f = frame_pxMAT, cell_n_f = cell_n) #* this the funtion that should be modified
			#>

		#< Loop for blacking out of the incorrect strain
		if blackout == True:
			blackout_table = temp_thresh_table.loc[temp_thresh_table["Myo1Identity"] != myo_chan].reset_index(drop = False)
			blackout_table["Cell_n"] = pd.Series(blackout_table["Cell_Barcode"]).apply(lambda x: int(x[x.find("c"):]))
			blackout = frame_pxMAT

			for b_out in blackout_table["Cell_n"]:
				a,b = (frame_pxMAT == b_out).nonzero() #* Ge the x and y locations where the mask is equal to the cell of interest
				x1 = np.min(a)
				x2 = np.max(a)
				y1 = np.min(b)
				y2 = np.max(b)
				blackout[x1:(x2+1), y1:(y2+ 1)] = 0

			blackout[np.nonzero(blackout)] = 1 #* Set all locations which have not been set to zero as 1, including the non-cell space. The only spots that are set to zero by this point should be the incorrect strain

			img = cv2.bitwise_and_(img, img, mask = blackout) #* Perform a mask mutliplication
			#? Add padding: https://stackoverflow.com/questions/43391205/add-padding-to-images-to-get-them-into-the-same-shape

			try:
				del blackout
			except:
				pass
			#>

	#* Select the spots where the pixels are grater than the threshold
	image_threshold_higher = (img > threshold_higher)
	image_threshold_lower =  (img < threshold_lower)

	if composite == False:
		#* Save the modified images in non-composite form
		image_threshold_higher.save(f"{image_frame}_image_threshold_higher.tif")
		image_threshold_lower.save(f"{image_frame}_image_threshold_lower.tif")
	elif composite == True:
		#* Save the modified composite images
		image_threshold_higher.save(f"{image_position}_composite_threshold_higher.tif")
		image_threshold_lower.save(f"{image_position}_composite_threshold_lower.tif")

# def create_thresh_table():
# 	pd.read_csv()

def visualization_manager(threshold_table, blackout = False): #, Manager for the visual cell creation, allows the results to called with the speciic params
	image_prompt_visual = input("What is the image that you would like to create a threshold set for?")
	image_prompt_visual = image_prompt_visual.split(", ")

	if blackout == True:
		seg_index_path = input("Where is the cell index?")
		seg_index_path = global_functions.slash_switch(seg_index_path) #, define the global functions group
		seg_index = pd.read_csv(seg_index_path).set_index(["Unique_frame"], drop = True)

	temp_thresh_table = threshold_table.loc[threshold_table[:].isin(image_prompt_visual)]
	x = Parallel(n_jobs= pr, verbose= 100)(delayed(visual_images)(temp_thresh_table, i) for i in temp_thresh_table)
	return(x) #* this will show the completion log

#! the calling should print xand save the log to a text file