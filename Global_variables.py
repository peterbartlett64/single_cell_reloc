
#%%
from os import path, chdir, cpu_count
#%%
def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
	new = path.replace(os.sep, '/')
	return (new)

def global_vars():
	# global analyze
	# global microfluidic_results
	# global post_path
	# global cpu_se
	# global percentiles

	a = 0
	while a == 0:
		analyze = input("Where are the images stored?")
		analyze = slash_switch(analyze)
		if path.isdir(analyze) == True:
			a = 1
		elif analyze == '':
			a = 1
		else:
			print("Non permissible")
			pass
	m = 0
	while m == 0:
		microfluidics_results =  input("Microfluidics results folder")
		microfluidics_results = slash_switch(microfluidics_results)
		if path.isdir(microfluidics_results) == True:
			m = 1
		elif microfluidics_results == '':
			m = 1
		else:
			print("Non permissible")
			pass

	p = 0
	while p == 0:
		post_path = input("Post_quant results folder")
		post_path = slash_switch(post_path)
		if path.isdir(post_path) == True:
			p = 1
		elif post_path == '':
			p = 1 
		else:
			print("Non permissible")
			pass
	
	cpu_number = cpu_count()
	
	cp = 0
	while cp == 0:
		try:
			cpu_se  = int(input(f"Sytem has {cpu_number} cores. How many would you like to use"))
		except ValueError:
			continue
		if cpu_se <= cpu_number:
			cp = 1
		elif analyze == '':
			cp = 1
		else:
			print(f"That is an invalid number of cpu cores. Please select a number <= {cpu_number}")
	
	perc = 0
	while perc == 0:
		percentiles_input = input("Pipeline set to auto-run on 95th and 99th pecentile. Would you like to run on EXTRA intensity? If yes, enter the numerical intensity. If no, press ENTER")
		if percentiles_input == "":
			percentiles = [95, 99]
			perc = 1
		else:
			percentiles_input = str.split(percentiles_input, ",")
			percentiles = [int(i) for i in percentiles_input]

			if percentiles == [95, 99]:
				perc = 1
			else:
				if type(percentiles) is str or int:
					percentiles = [percentiles, 95, 99] #* The system shoudl still run on 95h and 99th percentiles
					perc = 1
				elif type(percentiles) is list:
					percentiles = [percentiles] + [95, 99]
					perc = 1
				else:
					print("Non permissible")
					pass

	global Global_variables
	Global_variables = {"analyze": analyze,
						"microfluidics_results": microfluidics_results,
						"post_path": post_path,
						"cpu_se": cpu_se,
						"percentiles": percentiles}

	return(None)
	# return(f"percentiles are {percentiles}.")
	# return(f"Analyze at {analyze},\nmicrofluidics_results at {microfluidics_results},\npost_path at {post_path},\npercentiles are {percentiles}.\nRunning with {cpu_se} cpu cores out of {cpu_number}")

def global_continue(): #* Confirm that the paths given are correct
	print(f"Analyze at {Global_variables['analyze']},\nmicrofluidics_results at {Global_variables['microfluidics_results']},\npost_path at {Global_variables['post_path']},\npercentiles are {Global_variables['percentiles']}.\nRunning with {'cpu_se'} cpu cores out of {Global_variables['cpu_se']}")
	cont = input("Would you like to continue?")
	if cont == "Yes" or "yes" or 1:
		cont = 1
	return(cont)
	
	#//
	#! None
	#! Analyze at ,
	#! microfluidics_results at ,
	#! post_path at ,
	#! percentiles are [[12], 95, 99].
	#! Running with cpu_se cpu cores out of 12
	#//
	
if __name__ == "__main__": #* Allow the program to be run individually to chang the global variables
	print(global_vars())
	cont = 0
	while cont != 1:
		cont  = global_continue()
	del cont

# %%
