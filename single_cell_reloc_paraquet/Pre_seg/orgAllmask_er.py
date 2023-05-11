import os
import pandas as pd
import time
import single_cell_reloc_parquet.global_functions.global_variables as gv

os.chdir(microfluidics_results)

def orgAllmask_er():
	orgmaskpaths = []
	count = 0
	fluors = "GFP_mKO_mKa"
	exclude = list(['No_fluor', 'No_flour', 'mKO_mKa']) #TODO: Change in final version to not expilicity subset like this
	for root, dirs, files, in os.walk(os.getcwd(), topdown= True):
		[dirs.remove(d) for d in list(dirs) if d in exclude]
		for name in files:
			if name.endswith(".mat") and name.startswith("mask"):
				orgmaskpaths.append({'Path': os.path.join(root, name)})
				count = count + 1
				print(count, end="\r")
			else:
				pass
	print(count)
	del count

	def f_Pos_mask(z):
		start = z.find('p') #Note: The shift forward is dependant upon how you write out position
		end = z.find('GFP')-1  #Make sure that the 'n' is present in the array to confirm no place 0 has been lost
		pos = z[start:end]
		# if len(pos) <= 8:
		#     return(pos)
		# else:
		#     pass
		return(pos)


	def mask_inf(m):
		num1 = m.find('ask_')+5
		num2 = m.find('.mat')
		num_ = ("f"+ m[num1:num2])
		return(num_)

	def mask_hour(m):
		hourS = m.find("_")+2
		hourE = m.find("_p")
		return(m[hourS:hourE])

	def mask_barcode(m):
		start =m.find('d02')
		end = m.find('GFP') -1
		return (m[start:end])

	def f_maskdate(x):
		dstart = x.find('d0')
		dend = dstart + 5
		expdate = (x[dstart:dend])
		return(expdate)

	orgAllmasks = pd.DataFrame(orgmaskpaths)
	orgAllmasks["Date"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_maskdate)
	orgAllmasks["Position_ID"] = pd.Series(orgAllmasks.iloc[:,0]).apply(f_Pos_mask)
	orgAllmasks["MaskNum"] = pd.Series(orgAllmasks.iloc[:,0]).apply(mask_inf)
	# Allmasks["Hour/Run_t"] = pd.Series(Allmasks.iloc[:,0]).apply(mask_hour)
	# Allmasks["Unique_pos_hour"] = Allmasks['Date'] + "h" + Allmasks["Hour/Run_t"] + Allmasks['Position_ID']

	orgAllmasks["Unique_pos"] = pd.Series(orgAllmasks.iloc[:,0]).apply(mask_barcode)
	orgAllmasks["Unique_frame"] = orgAllmasks["Unique_pos"] + orgAllmasks['MaskNum'] #By performing this way, can avoid one extra step per index

	del orgmaskpaths

	orgAllmasks.dropna(inplace=True) # temp drop because there are some extra segmentations
	orgAllmasks.reset_index(inplace=True, drop = True)
	return(orgAllmasks)

if __name__ == "__main__":
	Global_variables  = gv.global_variables()
	orgAllmasks = orgAllmasker()