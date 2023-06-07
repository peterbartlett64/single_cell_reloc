#%%
import os
import PIL
import numpy as np
from joblib import Parallel, delayed
import cv2
from matplotlib import pyplot as plt
import pandas as pd
# import single_cell_reloc_paraquet.global_functions.global_variables as gv
from PIL import Image

def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
	new = path.replace(os.sep, '/')
	return (new)

#%%
###############################< Code modified from Medium Aricle: "Image Shifting using NumPy from Scratch" by Sameer
#* https://medium.com/analytics-vidhya/image-shifting-using-numpy-from-scratch-8bd52663da52

def pad_vector(vector, how, depth, fill_value): #! Check to confirm that the fill value has not changed raw_val
	vect_shape = vector.shape[:2]
	type_use = vector.dtype #* Added for variable bit depth
	if (how == 'upper') or (how == 'top'):
		pp = np.full(shape=(depth, vect_shape[1]), fill_value=fill_value, dtype= type_use) #* Modified for variable bit depth
		pv = np.vstack(tup=(pp, vector))
	elif (how == 'lower') or (how == 'bottom'):
		pp = np.full(shape=(depth, vect_shape[1]), fill_value=fill_value, dtype= type_use)
		pv = np.vstack(tup=(vector, pp))
	elif (how == 'left'):
		pp = np.full(shape=(vect_shape[0], depth), fill_value=fill_value, dtype= type_use)
		pv = np.hstack(tup=(pp, vector))
	elif (how == 'right'):
		pp = np.full(shape=(vect_shape[0], depth), fill_value=fill_value, dtype= type_use)
		pv = np.hstack(tup=(vector, pp))
	else:
		return vector
	return pv

def read_this(image_file, gray_scale=True): #? This is fine
	image_src = cv2.imread(image_file, cv2.IMREAD_ANYDEPTH) #* Modified to read in any depth
	if not gray_scale:
		image_src = cv2.cvtColor(image_src, cv2.COLOR_BGR2RGB)
	return image_src

def shifter(vect, y, y_, fill_value):
	if (y > 0):
		image_trans = pad_vector(vector=vect, how='upper', depth=y_, fill_value =fill_value)
		image_trans = image_trans[:-y_,:]#.
	elif (y < 0):
		image_trans = pad_vector(vector=vect, how='lower', depth=y_, fill_value= fill_value)
		image_trans = image_trans[y_:,:]#.
	else:
		image_trans = vect
	return image_trans

def shift_image(image_src, at):
	x, y = at
	x_, y_ = abs(x), abs(y)

	fill_value = np.mean(image_src)
	if (x > 0):
		left_pad = pad_vector(vector=image_src, how='left', depth=x_, fill_value=fill_value)
		image_trans = shifter(vect=left_pad, y=y, y_=y_, fill_value = fill_value)
		image_trans = image_trans[:,:-x_]#.
	elif (x < 0):
		right_pad = pad_vector(vector=image_src, how='right', depth=x_, fill_value = fill_value)
		image_trans = shifter(vect=right_pad, y=y, y_=y_, fill_value = fill_value)
		image_trans = image_trans[:,x_:] #.
	else:
		image_trans = shifter(vect=image_src, y=y, y_=y_)

	return image_trans

def translate_this(image_file, at, with_plot=False, gray_scale=True):
	if len(at) != 2: return False

	image_src = read_this(image_file=image_file, gray_scale=gray_scale)

	if not gray_scale:
		r_image, g_image, b_image = image_src[:, :, 0], image_src[:, :, 1], image_src[:, :, 2]
		r_trans = shift_image(image_src=r_image, at=at)
		g_trans = shift_image(image_src=g_image, at=at)
		b_trans = shift_image(image_src=b_image, at=at)
		image_trans = np.dstack(tup=(r_trans, g_trans, b_trans))
	else:
		image_trans = shift_image(image_src=image_src, at=at)


	#. Make copy
	# dot = image_file.find(".")
	# name = image_file[:dot]
	# image_ext = image_file[dot:]
	# temp_name = name + "_test" + image_ext

	# cv2.imwrite(temp_name, image_trans)
	cv2.imwrite(image_file, image_trans) #. This was added to the code by Peter. Usually you would want a backup, but in this case, there are multiple stored copies of the orignal images for safety. It is easier for downstream to change the iamges in the actual "Analyze" folder

	if with_plot:
		cmap_val = None if not gray_scale else 'gray'
		fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 20))

		ax1.axis("off")
		ax1.title.set_text('Original')

		ax2.axis("off")
		ax2.title.set_text("Translated")

		ax1.imshow(image_src, cmap=cmap_val)
		ax2.imshow(image_trans, cmap=cmap_val)
		return True

	return image_trans
################################################################################>
#%%
def miniImgIndex():
	def f_Frame_int(z):
		start = z.find('_time')+5
		end = z.find('_time')+9
		return(int(z[start:end])) #* This is a modification of normal string returned version
	def get_file_name(x):
		os.path.join(x)
		return(name)
	def f_skip_test(x):
		if "test" in x:
			return(None)
		else:
			return(x)

	count = 0
	mini_imgIndex = []
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name.endswith(".tif"):
				mini_imgIndex.append({'File_name': os.path.join(name)})
				count = count + 1
				print(count, end="\r")
			else :
				pass
		break

	del dirs
	del root

	mini_imgIndex = pd.DataFrame(mini_imgIndex)
	mini_imgIndex['Frame'] = pd.Series(mini_imgIndex.iloc[:,0]).apply(f_Frame_int) #* Skip where the file was a test to avoid recursive testing
	mini_imgIndex['File_name'] = pd.Series(mini_imgIndex.iloc[:,0]).apply(f_skip_test)
	mini_imgIndex.dropna(inplace=True)
	# mini_imgIndex['File_name'] = pd.Series(mini_imgIndex[:,0]).apply(lam) #. Change this line
	# mini_imgIndex['File_name'] = pd.Series(mini_imgIndex[:,0]).apply(lambda x: os.path.basename(x).split('/')[-1])

	return(mini_imgIndex)

#%%

# i = "D:/Microfluidics/References/Reference Shifts/d0215r1/Sub_position070300_time0031 - Copy.tif"
# vector_temp = [-455, -4]
# translate_this(image_file=i, at =vector_temp, with_plot=True)
#%%
if __name__ == "__main__":
	folder = input("Where is the directory to be done?")
	# vector = input('What are the vector aspects {x, y}').split(sep= ", ")
	# vector[0] = int(vector[0]) * -1
	# vector[1] = int(vector[1]) * -1
	# vector = vector_temp
	vector = [-183, -17]
	stop = int(input("What is the max frame to be shift"))
	folder = slash_switch(folder)
	os.chdir(folder)
	try:
		# past_folder
		# if past_folder == folder:
		# 	print(f"The '{past_folder}' has already been run!")
		pd.read_csv("already_shifted.txt")
		print(f"The '{folder}' has already been run!")
	except:
		pass

	parallel_state = input("Should the files be shifted in parallel? [y/n]")
	if parallel_state == "y" or parallel_state == "yes" or parallel_state =="Y" or parallel_state == "Yes" or parallel_state == "YES":
		parallel_state = True
	else:
		parallel_state = False
	use_cores_len = os.cpu_count()
	mini_imgIndex = miniImgIndex() #* Run search in current folder
	mini_imgIndex_subset = mini_imgIndex.loc[mini_imgIndex['Frame'] == stop]

	if parallel_state:
		Parallel(n_jobs=use_cores_len, verbose = 100, prefer='threads')(delayed(translate_this)(image_file = i, at = vector, with_plot=False, gray_scale=True) for i in mini_imgIndex_subset['File_name']) #* This should prefer threads as it largely IO limited.
	else:
		for i in mini_imgIndex_subset['File_name']:
			translate_this(image_file= i, at = vector, with_plot=False, gray_scale=True )
	already_shifted = pd.Series(vector)
	already_shifted.to_csv("already_shifted.txt")#* this will save stuipid mistakes if the folder has already been run
	# # past_folder = folder