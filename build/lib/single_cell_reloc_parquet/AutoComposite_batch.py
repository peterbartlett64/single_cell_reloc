
#%%
import imagej
import os

os.path.

#%%
ij = imagej.init('C:/Users/pcnba/OneDrive/Desktop/Fiji.app')
channels = ['GFP', 'mKa', 'mKO']
for c in channels:
	temp_path = os.path.join(Global_Variables["analyze"], "/{c}") #*This is defined in the gloabl functions
	os.mkdir(temp_path)

def path_combine_macro_er(postion_reg, channel, directory):
	macro_combine = f'File.openSequence("{directory}", " filter=({channel}_position{postion_reg})"); \nrun("Save");\nsaveAs("Tiff", "{Global_Variables["analyze"]}/{channel}_{postion_reg}.tif");'
	return(macro_combine)

def GPU_composite():
	marco_composite = f'F'
	return(marco_composite)

test_position = "440200"
test_path = "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code"
macro = path_combine_macro_er(test_position, test_path)
print(macro)

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
	channels = ['GFP', 'mKa', 'mKO']
	for c in channels:
		temp_path = os.path.join(Global_Variables["analyze"], "/{c}") #*This is defined in the gloabl functions
		os.mkdir(temp_path)
	



# def random_sample_training(Channel):
# 	new_path = os.path.join(Img_root, f"/Training/{Channel}")
# 	os.chdir(new_path)

# 	return()