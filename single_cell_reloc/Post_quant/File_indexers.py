import pandas as pd
import os
import datetime

def Quantification_index_er():
	#! commented out because only Quant_ALL files are present
	Quantification_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith(".csv") and name.startswith("Quantification_d"):
				Quantification_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		# break #This makes the program run non-recursively and not decend into daughter folders

	Quantification_index = pd.DataFrame(Quantification_index)

	year = str(datetime.datetime.now().year)
	decade = year[:3] # This assumes that the analysis is done within the same decade as starting
	del year

	def f_Position_ID(z):
		start = z.find('tion_')+5 #Note: The shift forward is dependant upon how you write out position
		end = z.find("_"+ decade)-5   #Assume that this pipeline will only be used this century! Make sure that the 'n' is present in the array to confirm no place 0 has been lost
		return(z[start:end])

	def f_Frame(z):
		start = z.find('tion_')+5
		end = z.find('_' + "20")
		return(z[start:end])

	def f_expdate(x):
		dend = x.find('r')
		expdate = x[0:dend]
		return(expdate)

	def Mod_epoch(f):
		t = os.path.getmtime(f)
		return(t)
		# return(time.ctime(t))

	def f_non_sync(p):
		substring = ".sync"
		if substring in p:
			return(None)
		else:
			return(p)

	def epoch_convert(timestamp):
		date_time = datetime.datetime.fromtimestamp(timestamp)
		d = date_time.strftime("%m/%d/%Y")
		return(d)

	# Quantification_index["Pad"] =
	Quantification_index["PositionID"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Position_ID)
	Quantification_index["Date"] = pd.Series(Quantification_index["PositionID"]).apply(f_expdate)
	Quantification_index["Frame"] = pd.Series(Quantification_index.iloc[:,0]).apply(f_Frame)

	#* NEWLY ADDED
	Quantification_index["Mod_epoch"] = pd.Series(Quantification_index["Path"]).apply(Mod_epoch)
	Quantification_index["Mod_date"] = pd.Series(Quantification_index["Mod_epoch"]).apply(epoch_convert)

	Quantification_index["Max_epoch_frame"] = Quantification_index.groupby("Frame")["Mod_epoch"].transform('max') #? This isn't sylish but works
	Quantification_index = Quantification_index.loc[Quantification_index["Mod_epoch"] == Quantification_index["Max_epoch_frame"]]
	Quantification_index.drop(columns=["Max_epoch_frame"])


	if Experimental_info["data_subset"] == True:
		Quantification_index = Quantification_index[Quantification_index["PositionID"].isin(instances)]
		print(f"Post-quant will be performed for positions which correspond to the following pos_barcodes: {instances}")
	else:
		pass
	Quantification_index.sort_values(by = "PositionID", inplace = True)
	Quantification_index.to_csv("Quantification_index.csv")
	return(Quantification_index)

def combine_pos(pos):
	Quant_frame_comb = pd.DataFrame([])
	subset = Quantification_index.loc[Quantification_index["PositionID"] == pos]
	subset_s = subset.sort_values(by = "Frame")
	for qf in range(len(subset)):
		try:
			q = pd.read_csv(subset_s.iloc[qf,0])
		except:
			continue
		Quant_frame_comb = pd.concat([Quant_frame_comb, q])
	Quant_frame_comb.to_csv(f"Quant_{pos}_ALL.csv")
	return(f"{pos} frame merging is complete")

def Quant_ALL_index_er():
	year = str(datetime.datetime.now().year) #this assumes the analysis is completed within the same decade of starting. A reasonable expectation, but something to check if running aroudn the new year
	decade = year[:3]
	del year

	def f_Position_ID_qALLi(z):
		start = z.find('ant_')+4
		end = z.find("_ALL")
		return(z[start:end])


	Quant_ALL_index  = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith("_ALL.csv") and name.startswith("Quant"): # fix naming
				Quant_ALL_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		break #This makes the program run non-recursively and not decend into daughter folders

	Quant_ALL_index = pd.DataFrame(Quant_ALL_index)
	Quant_ALL_index["PositionID"] = pd.Series(Quant_ALL_index.iloc[:,0]).apply(f_Position_ID_qALLi)
	Quant_ALL_index.sort_values(by = "PositionID", inplace = True)
	Quant_ALL_index.to_csv("Quant_ALL_index.csv")
	return(Quant_ALL_index)

def Quant_prim_index_er():
	def f_Position_ID_qALLi(z):
		start = z.find('ant_')+4
		end = z.find("_prim")
		return(z[start:end])

	# def f_col(p):
	# 	Col = int(p[p.find("p")+1:-4])
	# 	return(Col)

	# def f_expdate(x):
	# 	dend = x.find('r')
	# 	expdate = x[0:dend]
	# 	return(expdate)

	# def f_run(x):
	# 	s = x.find("r")+1
	# 	e = s+1
	# 	r = x[s:e]
	# 	return(r)
	os.chdir(post_path)
	Quant_prim_index = []
	count = 0
	for root, dirs, files, in os.walk(os.getcwd()):
		for name in files:
			if name.endswith("mary.csv") and name.startswith("Quant"): # fix naming
				Quant_prim_index.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
		break #This makes the program run non recursively and not decend into daughter folders

	Quant_prim_index = pd.DataFrame(Quant_prim_index)
	Quant_prim_index["PositionID"] = pd.Series(Quant_prim_index.iloc[:,0]).apply(f_Position_ID_qALLi)
	# Quant_prim_index["Position"] = pd.Series(Quant_prim_index["Path"]).apply(f_Position_ID_qALLi)
	# Quant_prim_index["Run"] = pd.Series(Quant_prim_index["Position"]).apply(f_run)
	# Quant_prim_index["Col"] = pd.Series(Quant_prim_index["Position"]).apply(f_col)
	# Quant_prim_index["Date"] = pd.Series(Quant_prim_index["Position"]).apply(f_expdate)

	Quant_prim_index.sort_values(by = "PositionID", inplace = True)
	Quant_prim_index.to_csv("Quant_prim_index.csv")# , index = False)
	return(Quant_prim_index)





def read_condition_informaiton():
	Condition_information = pd.read_excel("MicrofluidicsMap_wCol.xlsx", sheet_name='ProteinLocations', dtype={'Date' : str})
	Condition_information = Condition_information.drop(columns=["Notes", "Predicted localization Change"]).dropna()
	Condition_information.drop(Condition_information.columns[Condition_information.columns.str.contains('Unnamed',case = False)],axis = 1, inplace= True)

	time_per_frame = Global_variables["timepoint_space"]
	year = str(datetime.datetime.now().year)
	decade = year[:3]
	del year

	def f_convert_date(d):
		d = str(d)
		s = d.find(decade)+5
		return("d" + d[s:])
	Condition_information["Date"] = pd.Series(Condition_information["Date"]).apply(f_convert_date)
	return(Condition_information)

