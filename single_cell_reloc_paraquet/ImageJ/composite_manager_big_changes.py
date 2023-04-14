#%%
import pandas as pd
import os
import imagej
# from joblib import Parallel, delayed
import single_cell_reloc.global_functions.global_variables as glv
import logging
from scyjava import jimport
from joblib import Parallel, delayed
import scyjava
import random
#%%
#! Correct this to prompt user for the local instalation of imageJ

#// def prep_er() -> None:
#// 	global ij
#// 	ij = imagej.init('C:/Users/pcnba/OneDrive/Desktop/Fiji.app', mode='headless') #* Load the local instalation of imageJ.

#// 	global imgaIndex
#// 	imgIndex = pd.read_parquet('imgIndex.parquet')
#// 	return(None)
def ij_initiate()->str:
	try:
		ij
		text = 'IJ is already loaded'
	except NameError:
		exec("global ij\nij = imagej.init('C:/Users/Nikon/Desktop/Fiji.app')")
		#!ij = imagej.init('C:/Users/pcnba/OneDrive/Desktop/Fiji.app')")
		#* This seems roundabout but it gets around the automatic python syntax checking which would fail if global is called unwrapped
		text = 'IJ has been initiated'
	return(text)

def denoise_macro_define(a_Ka, d_Ka, s_Ka, a_KO, d_KO, s_KO, directory, pos_core) ->str:
	composite_script_mKa = f"""
	img_KO = run("Image Sequence...", "dir={directory} filter=(mKa_{pos_core}) sort");
	run("PureDenoise", "parameters = '3 1' estimation = 'Manual {a_Ka} {d_Ka} {s_Ka}' ")
	rename("Output_img_Ka")
	"""
	composite_script_mKO = f"""
	img_KO = run("Image Sequence...", "dir={directory} filter=(mKO_{pos_core}) sort");
	run("PureDenoise ", "parameters = '3 1' estimation = 'Manual {a_KO} {d_KO} {s_KO}' ")
	rename("Output_img_KO")
	"""
	return(composite_script_mKa, composite_script_mKO)

#// ij.py.run_plugin("PureDenoise...", args={''})
#// def denoise_imageJ(position, channel):
#// 	FolderOpener_py = jimport('ij.plugin.FolderOpener')
#// 	BIG_denoise_py = jimport('ij.plugin.PureDenoise')
#// 	mKa_set = FolderOpener_py.open(f"D:/Testing_myo_14_11_2022/d0222r1p940200/Testing_comp/", " filter=({channel}_*)");
#// 	mKa_set.setTitle(f"{channel}_{position}")

def series_macro_define_run(directory, position_full, channel, analyze)->str:
	pos_core = revert_pos_barc(position_full)
	directory = back_slash(directory)
	series_script = f'''img = File.openSequence("{directory}/", "filter=({channel}_position{pos_core}_time)");
	saveAs("Tiff", "{analyze}/denoise/{channel}/{position_full}.tif"); close();''' #// It was img.close()
	try:
		ij.py.run_macro(series_script)
		return(["PASS",directory, pos_core, channel])
	except:
		return(["FAIL",directory, pos_core, channel])

def series_macro_define_run_parallel(position_full, channel, analyze)->str:
	directory = os.path.dirname(imgIndex.loc[imgIndex["Unique_pos"] == position_full]["Path"][0])
	pos_core = revert_pos_barc(position_full)
	series_script = f'''img = File.openSequence("{directory}/", "filter=({channel}_position{pos_core}_time)");
	saveAs("Tiff", "{analyze}/denoise/{channel}/{position_full}.tif"); close();'''
	# ij = imagej.init()
	try:
		ij.py.run_macro(series_script)
		# ij.dispose()
		return(["PASS",directory, pos_core, channel])
	except:
		# ij.dispose()
		return(["FAIL",directory, pos_core, channel])

def denoise_composite_manger():
	#// Ka_denoise_params = input("FOR mKa: define: Alpha, delta detector offset, standard deviation of gausian noise")
	#// params_Ka  = str.split(Ka_denoise_params, ", ")
	#// if len(params_Ka) != 3:
	#// 	Ka_denoise_params = input("define: Alpha, delta detector offset, standard deviation of gausian noise")
	#// else:
	#// 	pass

	#!These could be defined in the Global but it is getting crowded
	a_Ka = 1.86 #params_Ka[0]
	d_Ka = 131 #params_Ka[1]
	s_Ka = 12 #params_Ka[2]

	#// KO_denoise_params = input(" FOR mKO: define: Alpha, delta detector offset, standard deviation of gausian noise")
	#// params_KO  = str.split(KO_denoise_params)
	#// if len(params_KO) != 3:
	#// 	KO_denoise_params = input("Alpha, delta detector offset, standard deviation of gausian noise")
	#// else:
	#// 	pass

	a_KO = 15 #params_KO[0]
	d_KO = 266 #params_KO[1]
	s_KO = 26 #params_KO[2]

	os.chdir(Global_Variables["microfluidics_results"]) #* When this is called, it will have been defined
	# imgIndex = pd.read_parquet('imgIndex.parquet')
	imgIndex = pd.read_csv('imgIndex.csv')
	positions = imgIndex['Unique_pos'].unique()

	# if ij.WindowManager.getIDList() is None:
	# 	ij.py.run_macro('newImage("dummy", "8-bit", 1,1,1);')

	for p in positions:
		try:
			original_core = revert_pos_barc(position = p)
			temp_directory = os.path.dirname(imgIndex.loc[imgIndex["Unique_pos"] == p]["Path"][0])
			# temp_directory = os.path.dirname(temp_path)
			composite_script_mKa, composite_script_mKO = denoise_macro_define(a_Ka, d_Ka, s_Ka, a_KO, d_KO, s_KO, temp_directory, original_core)
			ij.py.run_macro(composite_script_mKa) #* confirm if this is entered properly
			ij.py.run_macro(composite_script_mKO)
			# composite_imageJ(original_core=original_core, directory=temp_directory, mKa_image_name= 'Output_img_Ka', mKO_image_name='Output_img_KO') #* I might be able to grab the position name pos_core from the input file path
			GPU_composite_imageJ(original_core=original_core, directory=temp_directory, mKa_image_name= 'Output_img_Ka', mKO_image_name='Output_img_KO')

			ij.window().clear()
		except:
			ij.window().clear()

			pass #! THis should be replaced with a logging function

def just_composite_manager():
	imgIndex = pd.read_csv('imgIndex.csv')
	positions = imgIndex['Unique_pos'].unique()

	for p in positions:
		original_core = revert_pos_barc(position = p)
		temp_directory = os.path.dirname(imgIndex.loc[imgIndex["Unique_pos"] == p]["Path"][0])
		GPU_composite_imageJ(orignal_core = original_core, directory=temp_directory, mKa_image_name = p, mKO_image_name=p)

def time_series_manger():
	os.chdir(Global_Variables["analyze"]) #* When this is called, it will have been defined
	# imgIndex = pd.read_parquet('imgIndex.parquet')
	# global imgIndex
	# imgIndex = pd.read_csv('imgIndex.csv')

	positions = imgIndex['Unique_pos'].unique()
	channels = imgIndex['Channel'].unique()

	nc = os.cpu_count()
	length_pos = len(positions) * len(channels)

	if length_pos <= nc:
		use_cores_len = nc
		Parallel_func = False
	else:
		Parallel_func = True

	Parallel_func = False #* This is temporary

	#, Create the folders to store in. Skip if exists already, allowing for batch to be split up without error
	try:
		denoise_folder = os.path.join(Global_Variables['analyze'], "denoise")
		os.mkdir(denoise_folder)
	except FileExistsError:
		pass

	if Parallel_func == False:
		for tc in channels:
			try:
				folder = os.path.join(denoise_folder, tc)
				os.mkdir(folder)
			except FileExistsError:
				pass
			for p in positions:
				temp_directory = os.path.dirname(imgIndex.loc[imgIndex["Unique_pos"] == p]["Path"].unique()[0])
				try:
					res = series_macro_define_run(directory=temp_directory,position_full=p, channel = tc, analyze=Global_Variables["analyze"])
					print(res)
				except:
					pass #! THis should be replaced with a logging function
	if Parallel_func == True:
		res = Parallel(n_jobs=use_cores_len, verbose = 100, prefer='threads')(delayed(series_macro_define_run_parallel)(position_full=p, channel = tc, analyze=Global_Variables["analyze"]) for p in positions for tc in channels) #* This should prefer threads as it largely IO limited.
		print(res)

def revert_pos_barc(position):
	original_core = position[position.find("p")+1:]
	return(original_core)

def composite_imageJ(original_core, directory, mKa_image_name, mKO_image_name) -> None:
	macro_composite_save = f"""
	res = imageCalculator("Average create stack", "{mKa_image_name}","{mKO_image_name}");
	rename("myoDenAdd")
	selectWindow(myoDenAdd)
	StackWriter.save(imp, "{directory}", "format=tiff name=Adddivide{original_core}_time");
	"""
	ij.py.run_macro(macro_composite_save) #*Run the defined macro above

	#?different way
	# ij.py.run_plugin("Calculator Plus", "i1=mKa i2=mKO.tif operation=[Add: i2 = (i1+i2) x k1 + k2] k1=0.5 k2=0 create";) #* This would be the requirement for FIJI, but it does not have the right compiler
	return(None)

def GPU_composite_imageJ(original_core, directory, mKa_image_name, mKO_image_name) -> None:
	macro_composite_save = f'''
	// add images
	image_mKa = open("{directory}/mKa/{mKa_image_name}.tif
	//image_mKa = "{mKa_image_name}";
	Ext.CLIJ2_push(image_mKa);
	image2 = open("{directory}/mKa/{mKO_image_name}.tif
	//image2 = "{mKO_image_name}";
	Ext.CLIJ2_push(image2);Fimg
	image3 = "Sum_images";
	Ext.CLIJ2_addImages(image_mKa, image2, image3);
	Ext.CLIJ2_pull(image3);

	// multiply image and scalar
	Ext.CLIJ2_push(image3);
	image4 = "Composite_scaled";
	scalar = 0.5;
	Ext.CLIJ2_multiplyImageAndScalar(image3, image4, scalar);
	image5 = Ext.CLIJ2_pull(image4);
	StackWriter.save(image5, "{directory}", "format=tiff name={original_core}_time start = 1");
	'''
	ij.py.run_macro(macro_composite_save) #*Run the defined macro above
	return(None)

#! Temporary functions
def position_range(): # Input in form "dMMDDrNpNN-NN", use the gloabl version of imgIndex
	range_string = input("Provide comma deliminated postion ranges in form dMMDDrNpNN-NN. (3-digit codes will be accepted)")
	position_list = imgIndex['Unique_pos'].unique()
	options = []
	range_string = range_string.split(", ")
	for x in range_string:
		delim = x.find("-")
		first_part = x[:x.find("p")]
		middle = int(x[x.find("p")+1:delim])
		end_part = int(x[delim+1:])
		for n in range(middle,end_part): #* This nested for loop shouldn't be too big of an issue
			if n < 100:
				options.append(first_part + "p" + str(n).zfill(2) + "0") #* must add a trailing zero
			else:
				options.append(first_part + "p" + str(n).zfill(3) + "0")
	accept = []
	for o in options:
		for p_t in position_list:
			if p_t.startswith(o):
				accept.append(p_t)
	return(accept)

#TODO: Move to global functions
def back_slash(s):
	new = s.replace('\\','/')
	return (new)
#%%
if __name__ == "__main__":
	# Global_Variables = glv.global_manager()
	Global_Variables = {
		"analyze": "F:/Microfluidics/Missing_Analyze2",
		"microfluidics_results": "F:/Microfluidics/Missing_RESULTS2",
		# "post_path": "C:/Users/pcnba/OneDrive/Documents/Testing",
		"subset": False,
		'subset_by': 'range',
		'subset_collection': '',
		"cpu_se": 16,
		"timepoint_space": 7.5,
		"percentiles": [95, 99],
		"multiplex": True}
	os.chdir(Global_Variables['microfluidics_results'])
	imgIndex = pd.read_csv('imgIndex.csv')
	sample_size = int(input("Sample Size?"))
	pos_list = imgIndex['Unique_pos'].unique()
	len_pos = len(pos_list)
	rand_ints = random.sample(range(0, len_pos), sample_size)
	samp_pos = pos_list[rand_ints]
	imgIndex = imgIndex.loc[(imgIndex['Channel']!='mTf')& (imgIndex['Channel']!='mTF') & (imgIndex['Channel']!='Sub') & (imgIndex['Channel']!='out')]
	# accepted = position_range() #! This is temporarily here
	Global_Variables['subset_collection'] = list(samp_pos)  #+  list(accepted)
	imgIndex = imgIndex.loc[imgIndex['Unique_pos'].isin(Global_Variables['subset_collection'])]
	print(ij_initiate())
	time_series_manger()
	print("Functions completed")

