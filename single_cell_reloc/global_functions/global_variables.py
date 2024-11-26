#%%
from os import path, cpu_count
import json
import os

#%%
def slash_switch(path): ## This function is currently unused but could be usefull in the future for the cwd setting
	new = path.replace(os.sep, '/')
	return (new)

def Experimental_info():
	while True:
		try:
			timepoint_space = float(input("What is the time between images?")) #* If the number is converted to a float, then it should be acceptable downstream
		except ValueError:
			continue
		conf_ts = input(f"To confirm, are there {timepoint_space} minutes between frames?")

		if conf_ts.lower() == "yes" or conf_ts.lower() == "y":
			break
		else:
			continue

	cell_tracking = True

	while True:
		user_input = input('Segment with fluorescent channels? True/False')
		if user_input.lower() == 'true':
			fluorescent_seg = True
			break
		elif user_input.lower() == 'false':
			fluorescent_seg = False
			break
		else:
			continue

	Experimental_info = {
		"timepoint_space" : timepoint_space,
		"cell_tracking" : cell_tracking,
		"fluorescent_seg" : fluorescent_seg
	}

	return(Experimental_info)

def global_vars():
	# global analyze
	# global microfluidic_results
	# global post_path
	# global cpu_se
	# global percentiles

	a = 0
	loop = 1
	while a == 0:
		if loop >= 3:
			try_again = input("There have been 3 failed atempts to input variable. Would you like to try again? [y/n]")
			if try_again.lower() == 'y' or try_again.lower() == 'yes':
				pass
			else:
				return('Failed to input all global variables')
		analyze = input("Where are the images stored?")
		analyze = slash_switch(analyze)
		if path.isdir(analyze) == True:
			a = 1
		elif analyze == '':
			a = 1
		else:
			print("Non permissible. The directory given does not exist\n")
			pass
		loop += 1

	m = 0
	loop = 1
	while m == 0:
		if loop >= 3:
			try_again = input("There have been 3 failed atempts to input variable. Would you like to try again? [y/n]")
			if try_again.lower() == 'y' or try_again.lower() == 'yes':
				pass
			else:
				return('Failed to input all global variables')
		microfluidics_results =  input("Microfluidics results folder")
		microfluidics_results = slash_switch(microfluidics_results)
		if path.isdir(microfluidics_results) == True:
			m = 1
		elif microfluidics_results == '':
			m = 1
		else:
			print("Non permissible\n")
			pass
		loop += 1

	p = 0
	loop = 1
	while p == 0:
		if loop >= 3:
			try_again = input("There have been 3 failed atempts to input variable. Would you like to try again? [y/n]")
			if try_again.lower() == 'y' or try_again.lower() == 'yes':
				pass
			else:
				return('Failed to input all global variables')
		post_path = input("Post_quant results folder")
		post_path = slash_switch(post_path)
		if path.isdir(post_path) == True:
			p = 1
		elif post_path == '':
			p = 1
		else:
			print("Non permissible\n")
			pass
		loop += 1

	cpu_number = cpu_count()

	cp = 0
	loop = 1
	while cp == 0:
		if loop >= 3:
			try_again = input("There have been 3 failed atempts to input variable. Would you like to try again? [y/n]")
			if try_again.lower() == 'y' or try_again.lower() == 'yes':
				pass
			else:
				return('Failed to input all global variables')
		try:
			cpu_se  = int(input(f"Sytem has {cpu_number} cores. How many would you like to use?"))
		except ValueError:
			continue
		if cpu_se <= cpu_number:
			cp = 1
		elif analyze == '':
			cp = 1
		else:
			print(f"That is an invalid number of cpu cores. Please select a number <= {cpu_number}\n")
		loop += 1

	perc = 0
	loop = 1
	while perc == 0:
		if loop >= 3:
			try_again = input("There have been 3 failed atempts to input variable. Would you like to try again? [y/n]")
			if try_again.lower() == 'y' or try_again.lower() == 'yes':
				pass
			else:
				return('Failed to input all global variables')
		percentiles_input = input("Pipeline set to auto-run on 95th and 99th pecentile. Would you like to run on EXTRA intensity? If yes, enter the numerical intensity. If no, press ENTER")
		if percentiles_input == "":
			percentiles = [95, 99]
			perc = 1
		else:
			percentiles_input = str.split(percentiles_input, ",")
			percentiles_input = [int(i) for i in percentiles_input]

			if percentiles_input == [95, 99]:
				percentiles = [95,99]
				perc = 1
			else:
				if type(percentiles_input) is list:
					percentiles = percentiles_input + [95, 99]

					perc = 1 #* Check that all values are indeed percentiles
					for p in percentiles:
						if 0 < p <= 100:
							pass
						else:
							print("Non permissible values given")
							perc = 0
							continue
				else:
					print("Non permissible values given")
		loop += 1

	sub = 0
	loop = 1
	while sub == 0:
		if loop >= 3:
			try_again = input(f"There have been {loop} failed atempts to input variable. Would you like to try again? [y/n]")
			if try_again.lower() == 'n' or try_again.lower() == 'no':
				return('Failed to input all global variables')
			else:
				pass

		else:
			pass

		subset = input("Is this experiment subseted? [y/n]")
		if subset.lower() == "yes" or subset.lower() == "y":
			subset = True
			subset_by = input("What should the analysis be subsetted by?[Date, Run, OR position_barcode (d<mmdd>r<run>p<xxxxxxx>)]")
			subset_collection = input(f"Enter the list of {subset_by}s that you would like to analyze. [Full structure is d<mmdd>|date|r<run>|run|p<xxxxxxx>|position| ]")
			try:
				subset_collection.split(", ") #TODO: Complete the confirmation that the entries are permissible
			except:
				continue
			sub =  1
		elif subset.lower() == "no" or subset.lower() == "n":
			subset = False
			subset_by = ''
			subset_collection = ''
			sub = 1
		else:
			loop += 1

	gap = 0
	while gap == 0:
		timepoint_gap = int(input("How many minutes between frames? [Enter integer]"))
		response = input(f"{timepoint_gap} minutes? [y/n]")
		if response.lower == "yes" or response.lower == "y":
			gap = 1
		else:
			pass
		del response

	mult = 0
	loop = 1
	while mult == 0:
		if loop >= 3:
			try_again = input("There have been 3 failed atempts to input variable. Would you like to try again? [y/n]")
			if try_again.lower() == 'y' or try_again.lower() == 'yes':
				pass
			else:
				return('Failed to input all global variables')

		multiplex = input("Does this experiment have muliplexed samples?")
		if multiplex.lower() == 'yes':
			multiplex = True
			mult = 1
		if multiplex.lower() == 'no':
			multiplex = False
			mult = 1 #* This will exit
		loop += 1

	global Global_variables
	Global_variables = {"analyze": analyze,
						"microfluidics_results": microfluidics_results,
						"post_path": post_path,
						"subset": subset,
						'subset_by': subset_by,
						'subset_collection': subset_collection,
						"cpu_se": cpu_se,
						"timepoint_space": timepoint_gap,
						"percentiles": percentiles,
						"multiplex": multiplex}

	os.chdir(microfluidics_results)
	with open("Global_variables.json", "w") as write_file:
		json.dump(Global_variables, write_file, indent=4) #*Write the global variables to a JSON file

	return(Global_variables)
	# return(f"percentiles are {percentiles}.")
	# return(f"Analyze at {analyze},\nmicrofluidics_results at {microfluidics_results},\npost_path at {post_path},\npercentiles are {percentiles}.\nRunning with {cpu_se} cpu cores out of {cpu_number}")

def global_continue(Global_variables): #* Confirm that the paths given are correct
	if Global_variables['multiplex'] == True:
		mult_a = "ARE"
	else:
		mult_a = "ARE NOT"

	if Global_variables["subset"] == True:
		wording = 'WILL'

		print(f"Analyze at {Global_variables['analyze']};\nmicrofluidics_results at {Global_variables['microfluidics_results']};\npost_path at {Global_variables['post_path']}; \nThere {wording} be subsetting by{Global_variables['subset_by']} and will include {Global_variables['subset_collection']}\nRunning with {'cpu_se'} cpu cores out of {Global_variables['cpu_se']};\nThere are {Global_variables['timepoint_gap']} minutes between frames;\npercentiles are {Global_variables['percentiles']};\nSamples {mult_a} multiplexed")

	else:
		wording = 'WILL NOT'

		print(f"Analyze at {Global_variables['analyze']};\nmicrofluidics_results at {Global_variables['microfluidics_results']};\npost_path at {Global_variables['post_path']}; \nThere {wording} be subsetting by{Global_variables['subset_by']} and will include {Global_variables['subset_collection']};\nRunning with {'cpu_se'} cpu cores out of {Global_variables['cpu_se']};\nThere are {Global_variables['timepoint_gap']} minutes between frames;\npercentiles are {Global_variables['percentiles']};\nSamples {mult_a} multiplexed")
	cont = input("Are all these values correct? [y/n]")
	if cont.lower == "yes" or cont.lower == "y":
		cont = 1
	return(cont)

def global_manager():
	global_vars()
	cont = 0
	while cont != 1:
		cont = global_continue()

	#// I don't know what this is from
	#! None
	#! Analyze at ,
	#! microfluidics_results at ,
	#! post_path at ,
	#! percentiles are [[12], 95, 99].
	#! Running with cpu_se cpu cores out of 12
	#//

#? There is a version that makes global varaibles and there is a version that stores the globals in a dictionary
if __name__ == "__main__": #* Allow the program to be run individually to change the global variables
	global_manager()
	print(Global_variables)
	# cont = 0
	# while cont != 1:
	# 	cont = global_continue()
	# del cont
