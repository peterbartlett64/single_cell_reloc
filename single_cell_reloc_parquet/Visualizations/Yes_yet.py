working = test.copy() #. full_data.copy()
working = working.sort_values(by = "ImageID")
df_all_rem = working.loc[working["Frames_post_treatment"] < 0]
#>
working = working.loc[working["Frames_post_treatment"] >= 0]

# def complex_yet(x):
# 	ind = x["Relocalized"].idxmax()
# 	does = x.loc[ind]
# 	x["yet"] = 0
# 	x.loc[:ind, "yet"] = 0
# 	x.loc[ind:, "yet"] = does
# 	return

def reloc_yet(x):
	ind = x.idxmax() #* Get the first occurrence of max value. This is 1 when relocalised, so the first time in the 'Relocalized' series where true
	does = x.loc[ind] #* Get the value of max value. ie. If max value is 0, then it never relocalizes, whereas if 1 then it does at some point
	x.loc[:ind] = 0 #* Set all values less than the max value as 0
	x.loc[ind:] = does #* Set all values after the max value as 1
	return(x)

def workaround(ind):
	return(working.loc[ind, "Unique_Frame"])

def does_workaround(ind):
	return(working.loc[ind,"Relocalized"])

working["Yet"] = working.groupby(by = "Cell_Barcode")["Relocalized"].transform(reloc_yet) # This repesents wether there has been relocaliztion yet
working["ind"] = working.groupby(by = "Cell_Barcode")["Relocalized"].transform('idxmax')
working["Does"] = pd.Series(working["ind"]).apply(does_workaround) #this will work for now but should make it come out of one of the other functions
working["When"] = pd.Series(working["ind"]).apply(workaround)
working.drop(columns='ind', inplace = True)

#. These need to be fixed
# working["Yes_yet"] = working.groupby('Unique_Frame')['Yet'].apply(lambda x: x[x == 1].count())
# working['No_yet'] = working.groupby('Unique_Frame')['Yet'].apply(lambda x: x[x == 0].count())

#* This will work for now
info_yes = working.set_index('Frame').Yet.eq(1).sum(level=0).astype(int).reset_index()
info_yes.set_index("Frame", inplace = True)
info_yes.rename(columns={"Yet": "Yes_yet"}, inplace = True)


info_no = working.set_index('Frame').Yet.eq(0).sum(level=0).astype(int).reset_index()
info_no.set_index('Frame', inplace = True)
info_no.rename(columns={"Yet": "No_yet"}, inplace = True)

newest = pd.merge(working, info_yes, left_on = "Frame", right_index = True)
newest = pd.merge(newest, info_no, left_on = "Frame", right_index = True)