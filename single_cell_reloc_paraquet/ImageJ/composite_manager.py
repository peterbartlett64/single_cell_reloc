#%%
#, This will only run with the Microfluidics_Pipe_highest environment
#%%
import pandas as pd
import os
import imagej
from joblib import Parallel, delayed
import single_cell_reloc.global_functions.global_variables as glv
import logging
from scyjava import jimport

#%%
#! Correct this to prompt user for the local instalation of imageJ

#// def prep_er() -> None:
#// 	global ij
#// 	ij = imagej.init('C:/Users/pcnba/OneDrive/Desktop/Fiji.app', mode='headless') #* Load the local instalation of imageJ.

#// 	global imgaIndex
#// 	imgIndex = pd.read_parquet('imgIndex.parquet')
#// 	return(None)

def macro_define(a_Ka, d_Ka, s_Ka, a_KO, d_KO, s_KO, directory, core):
	composite_script_mKa = f"""
	img_KO = IJ.FolderOpener({directory})
	rename("Input_img_Ka")
	selectWindow(Input_img_Ka)
	run("PureDenoise ", "parameters = '3 1' estimation = 'Manual {a_Ka} {d_Ka} {s_Ka}' ")
	rename("Output_img_Ka")
	"""
	composite_script_mKO = f"""
	img_KO = run("Image Sequence...", "dir={directory} filter=(mKO_{core}) sort");

	IJ.FolderOpener({directory})
	rename("Input_img_KO")
	selectWindow(Input_img_KO)
	run("PureDenoise ", "parameters = '3 1' estimation = 'Manual {a_KO} {d_KO} {s_KO}' ")
	rename("Output_img_KO")
	"""
	return(composite_script_mKa, composite_script_mKO)

#// args = {"":}
#// result = ij.py.run_script("Groovy", composite_script, args)
#// stats = result.getOUtput("stats")
#// macro_composite_save = f"""
#// """

#%%
def composite_manger():
	global ij
	ij = imagej.init('C:/Users/pcnba/OneDrive/Desktop/Fiji.app', mode='headless')
	# ij = imagej.init('C:/Users/pcnba/OneDrive/Desktop/ImageJ/ImageJ.exe', mode='headless')

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
			temp_path = imgIndex.loc[imgIndex["Unique_pos"] == p]["Path"][0]
			temp_directory = os.path.dirname(temp_path)
			composite_script_mKa, composite_script_mKO = macro_define(a_Ka, d_Ka, s_Ka, a_KO, d_KO, s_KO, temp_directory, original_core)
			ij.py.run_macro(composite_script_mKa) #* confirm if this is entered properly
			ij.py.run_macro(composite_script_mKO)
			composite_imageJ(original_core=original_core, directory=temp_directory) #* I might be able to grab the position name core from the input file path

			ij.window().clear()
		except:
			ij.window().clear()

			pass #! THis should be replaced with a logging function

#// ij.py.run_plugin("PureDenoise...", args={''})
#// def denoise_imageJ(position, channel):
#// 	FolderOpener_py = jimport('ij.plugin.FolderOpener')
#// 	BIG_denoise_py = jimport('ij.plugin.PureDenoise')
#// 	mKa_set = FolderOpener_py.open(f"D:/Testing_myo_14_11_2022/d0222r1p940200/Testing_comp/", " filter=({channel}_*)");
#// 	mKa_set.setTitle(f"{channel}_{position}")

def revert_pos_barc(position):
	# date = position[position.find("d")+1:position.find("r")]
	pos = position[position.find("p")+1:]
	original_core = pos + "time_"
	return(original_core)

def composite_imageJ(original_core, directory) -> None:
	macro_composite_save = f"""
	res = imageCalculator("Average create stack", "Output_img_Ka","Output_img_KO");
	rename("myoDenAdd")
	selectWindow(myoDenAdd)
	StackWriter.save(imp, "{directory}", "format=tiff name={original_core}_time");
	"""
	ij.py.run_macro(macro_composite_save) #*Run the defined macro above

	#?different way
	# ij.py.run_plugin("Calculator Plus", "i1=mKa i2=mKO.tif operation=[Add: i2 = (i1+i2) x k1 + k2] k1=0.5 k2=0 create";) #* This would be the requirement for FIJI, but it does not have the right compiler
	# ij.py.run_plugin('imageCalculator("Average create stack", "mKO.tif","mKO.tif");')
	return(None)


if __name__ == "__main__":
	# Global_Variables = glv.global_manager()
	Global_Variables = {
		"analyze": "C:/Users/pcnba/OneDrive/Documents/Testing",
		"microfluidics_results": "C:/Users/pcnba/OneDrive/Documents/Testing",
		"post_path": "C:/Users/pcnba/OneDrive/Documents/Testing",
		"subset": False,
		'subset_by': '',
		'subset_collection': '',
		"cpu_se": 16,
		"timepoint_space": 7.5,
		"percentiles": [95, 99],
		"multiplex": True}
	composite_manger()