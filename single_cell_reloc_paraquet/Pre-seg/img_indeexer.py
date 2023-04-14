#%%
import pandas as pd
import os
import single_cell_reloc.global_functions.global_variables as glv
glv.slash_switch
# def instances():
# 	global subdat
# 	if Experimental_info["data_subset"] == True:
# 		subdat = str(input("Enter position barcodes to be analyzed (comma seperated and in the form dMMDDpPosition)"))

# 		instances = subdat.split(", ")
# 		print(f"Analysis will be performed for positions which correspond to the following pos_barcodes: {instances}")
# 	else:
# 		pass

def imgIndex_er(image_path, microfluidics_results):
	os.chdir(image_path)

	def f_lazer(z):
		start = z.find('_position')-3
		end = z.find('_position')
		return(z[start:end])

	def f_Position_ID(z):
		start = z.find('_position')+9 #Note: The shift forward is dependant upon how you write out position
		end = z.find('_time')   #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
		return(z[start:end])

	def f_Frame(z):
		start = z.find('_time')+5
		end = z.find('_time')+9
		return(z[start:end])

	def f_expdate(x):
		dstart = x.find('_hr')-4
		dend = x.find('_hr')
		expdate = x[dstart:dend]
		expdate = expdate.translate({ord('-'): None})
		expdate = "d" + expdate
		return(expdate)

	def f_runnumb(d):
		rsta = d.find("hr")+2
		ren = d.find("hr")+6
		return(d[rsta:ren])

	# def f_non_den(p):
	# 	substring = "denoise"
	# 	if substring in p:
	# 		return(None)
	# 	else:
	# 		return(p)

	count = 0
	imgIndex = []
	for root, dirs, files in os.walk(os.getcwd()):
		for name in files:
			if name.endswith(".tif"):
				imgIndex.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else :
				pass

	imgIndex = pd.DataFrame(imgIndex)
	# imgIndex["Path"] = pd.Series(imgIndex.iloc[:,0]).apply(f_non_den)
	imgIndex['Date'] = pd.Series(imgIndex.iloc[:,0]).apply(f_expdate)
	imgIndex['PositionID'] = pd.Series(imgIndex.iloc[:,0]).apply(f_Position_ID)
	imgIndex['Frame'] = pd.Series(imgIndex.iloc[:,0]).apply(f_Frame)
	imgIndex['Channel'] = pd.Series(imgIndex.iloc[:,0]).apply(f_lazer)
	imgIndex['Hour/Run_t'] =pd.Series(imgIndex.iloc[:,0]).apply(f_runnumb)
	#imgIndex["temp_date_hour"] = imgIndex["Date"] + imgIndex["Hour/Run_t"]
	imgIndex["Hour/Run_t"] = imgIndex["Hour/Run_t"].apply(pd.to_numeric)

	z = imgIndex[["Date", "Hour/Run_t"]]
	z = z.drop_duplicates().reset_index()

	z["Run"] = 0
	z["for"] = 1
	z = z.sort_values(by = ["Date", "Hour/Run_t"], ascending = True)

	for q in range(len(z)):
		r = 1
		if z.iloc[q,:]['for'] == 1:
			if z.iloc[q,1] == z.iloc[(q-1),1]:
				z.iloc[q,3] = (z.iloc[q-1,3] + 1)
				z.iloc[q,4] = 0
			else:
				z.iloc[q,3] = r
				z.iloc[q,4] = 0
		else:
			pass

	z = z[["index", "Date", "Run"]] #only grab the columns which are usefull, and not the ones created just for the logical loops

	t = len(z) - 1
	imgIndex ["Run"] = 0
	for g in range(len(z)):
		if g == 0:
			i1 = 0
			i2 = int(z.loc[g, ["index"]]) - 1
			imgIndex.loc[i1:i2, ["Run"]]  = z.iloc[g, 2]
		if g == t:
			i1 = int(z.loc[g, ["index"]])
			i2 = len(imgIndex) - 1
			imgIndex.loc[i1:i2, ["Run"]]  = z.iloc[g, 2]
		else:
			g2 = g + 1
			i1 = int(z.loc[g, ["index"]])
			i2 = int(z.loc[g2, ["index"]]) -1
			imgIndex.loc[i1:i2, ["Run"]]  = z.iloc[g, 2]
	del z

	del count
	del dirs
	del files
	del name

	imgIndex["Run"] = imgIndex["Run"].astype('str')
	imgIndex["Hour/Run_t"] = imgIndex["Hour/Run_t"].astype("str")
	imgIndex["Unique_pos"] = imgIndex['Date'] +"r" + imgIndex["Run"] + "p" + imgIndex['PositionID']
	imgIndex["Unique_pos_hour"] = imgIndex['Date'] + "h" + imgIndex["Hour/Run_t"] + "p" + imgIndex['PositionID']
	imgIndex["Unique_frame"] = imgIndex["Unique_pos"] + "f" + imgIndex['Frame'] #By performing this way, can avoid one extra step per index
	imgIndex.set_index(["Unique_frame"], inplace = True)
	imgIndex

	os.chdir(microfluidics_results)
	series = imgIndex["Unique_pos"]
	series = pd.unique(series)
	l = len(series)

	# with open('Positions.txt', 'a+') as posit:
	#     for p in series:
	#         posit.write(f'{p} ')
	# posit.close()
	imgIndex.to_csv("imgIndex.csv")
	# imgIndex.to_paraquet("imgIndex.paraquet")
	return(imgIndex)

if __name__ == "__main__":
	analyze = glv.slash_switch(input("Where is anlyze/are the images stored?"))
	microfluidics_results = glv.slash_switch(input("Where is the data output?"))
	imgIndex_er(analyze, microfluidics_results)

# %%
