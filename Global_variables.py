def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
	new = path.replace(os.sep, '/')
	return (new)

def get_paths():
	
# post_path =  "C:/Users/Peter/Desktop/Subset_Images/JAN2/2022-01-03"
microfluidic_results = str(input("Microfluidics results folder"))
# microfluidic_results = "D:/Microfluidics/RESULTS"

post_path= str(input("Post_quant results folder"))
# post_path= "D:/Microfluidics/RESULTS/2022-03-09(APR3) - Copy"

pn = os.cpu_count()
pn = int(input(f"Sytem has {pn} cores. How many would you like to use"))

# info_microfluidics_results = "F:/Microfluidics/RESULTS"
# microfluidics_results= "C:/Users/Peter/Desktop/Subset_Images/2021-12-03"
percentile = int(input("What percentile would you like to run on"))

#Change the paths to the corrected version becuase of the windows input
microfluidic_results = slash_switch(microfluidic_results)
post_path = slash_switch(post_path)
os.chdir(post_path) # This would normally be the path set in the quantification scipt