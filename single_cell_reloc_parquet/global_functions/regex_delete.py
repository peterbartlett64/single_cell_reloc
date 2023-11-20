#%%
import os
import re
import pandas as pd

def switch_slash(path):
	new = path.replace(os.sep, '/')
	return (new)

print('CONTINUING THIS FILE WILL DESTROY TRACKING FOR MATCH CASES. IT WILL ALSO REMOVE ALL SEGMENTATION FILES UP TARGET POINT')
continue_state = input('Continue?')
if continue_state != 'y':
	directory_root = 'purposely crashing backup'
	exit()
else:
	directory_root = input("Root directory?")
directory_root = switch_slash(directory_root)
os.chdir(directory_root)
# #%%
# def delete_files_matching_pattern(directory, pattern):
# 	# Iterate over all files in the directory

# 	for root, dirs, files in os.walk(os.getcwd()):
# 		for filename in files:
# 			filepath = os.path.join(root, filename)

# 			# Check if the file matches the pattern
# 			if re.search(pattern, filename):
# 				# Delete the file
# 				os.remove(filepath)
# 				print(f"*Would have* Deleted file: {filepath}")
# 			else:
# 				print(f"Would not have deleted {filepath}")

# def delete_files_NOT_matching_pattern(directory, pattern):
# 	# Iterate over all files in the directory

# 	for root, dirs, files in os.walk(os.getcwd()):
# 		for filename in files:
# 			filepath = os.path.join(root, filename)

# 			# Check if the file matches the pattern
# 			if re.search(pattern, filename):
# 				print(f"Would not have deleted {filepath}")
# 			else:
# 				# os.remove(filepath)
# 				print(f"*Would have* Deleted file: {filepath}")

# # Directory to search in


# # Regular expression pattern

# pattern_number = r'^[A-Za-z]+_000([3-9]\d|\d{3,}).[A-Za-z]+$'  # Example pattern: starts with "example" and ends with ".txt"
# pattern_type = r'^(track)'

# for root, dirs, files in os.walk(os.getcwd()):
# 	for d in dirs:
# 		# delete_files_matching_pattern(directory = d, pattern = pattern_type)
# 		delete_files_NOT_matching_pattern(directory = d, pattern = pattern_number)
#%%
targ_number = int(input("What is the target number?"))
target_daterun = input("What is the target daterun? [d{mmdd}r{n}]")

file_list = []
count = 0
for root, dirs, files in os.walk(os.getcwd()):
	for name in files:
		file_list.append({'Path': os.path.join(root, name)})
		count = count + 1
		print(count, end="\r")
file_list = pd.DataFrame(file_list)
#%%

def get_number(x):
	try:
		post_label = x[x.find("GFP_mKO, mKa")+12:]
		number = post_label[post_label.find("_0")+1:post_label.find(".")]
		number = int(number)
	except:
		number = None
	return(number)

def get_type(x):
	post_label = x[x.find("GFP_mKO, mKa")+12:]
	if 'track' in post_label or 'TracX' in post_label:
		return("T")
	if 'CellX' in post_label:
		return("S")
	elif 'pxMAT' in post_label:
		return("T")
	elif 'Track' in post_label:
		return("T")
	else:
		return("S")

def get_daterun(z):
	start = z.find('d02') #Note: The shift forward is dependant upon how you write out position
	end = z.find('p')  #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
	pos = z[start:end]
	return(pos)

file_list['daterun'] = pd.Series(file_list['Path']).apply(get_daterun)
file_list['number'] = pd.Series(file_list['Path']).apply(get_number)
file_list['type'] = pd.Series(file_list['Path']).apply(get_type)

file_list = file_list.loc[file_list['daterun'] == target_daterun]
file_list = file_list.loc[(file_list["number"] <= targ_number) | (file_list['type'] == "T")]

for fp in file_list['Path']:
	os.remove(fp)
