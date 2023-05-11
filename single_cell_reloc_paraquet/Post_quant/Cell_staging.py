#%%
#, This is the updated version which relies on the output of the TracX composite
import pandas as pd
import os
from joblib import Parallel, delayed
from math import sqrt
import random
import openpyxl
import scipy.stats as stats
import numpy as np
import time

#%%
#, Temporary read in for testing
cell_cycle_df = pd.read_excel("C:/Users/pcnba/Desktop/Temp_data/2023111175439_CellCycleResults_test_composite_d0216r1p030300.xlsx", sheet_name='2023111175439_CellCycleResults_', usecols = ["track", "parent", "daughter", "tree", "generation", "age", "g2_st", "g2_e", "g2_dur", "div", "seg_st_track", "seg_st_daughter"]).set_index("track")

#%%
#? Might be able to incorparate a distance element in a future version. Distance_g1 = frame - <previous g2_e (with obvious issues capturing and output when it is is the first division)> ; Distance_g2 = frame - g2_st
#? Distance g1 G2_start -13========> 0 (G2_e) =========> 4
#? Distance g2 G2_start 0  ========> 13(g2_e) =========> 0

def stage_from_frame(cell, frame) ->int: #? I should try implementing this as a parrellel pandas
	sub = cell_cycle_df.loc[cell,:][["g2_st", "g2_e"]] #* Subset based on cell position_barcode and grab the relavant columns
	for i in range(len(sub)):
		#* Get required stage times to find where the current frame fits in
		g = sub.iloc[i:(i+2), :]
		g2_st1 = g.iloc[0]["g2_st"]
		g2_e1 = g.iloc[0]["g2_e"]
		try:
			g2_st2 = g.iloc[1]["g2_st"]
			g2_e2 = g.iloc[1]["g2_e"]

			#< Start of the multirow conditions within the loop. 1 ==> G1/S, 2 ==> G2/M
			if frame < g2_st1:
				return("G1")
			elif g2_st1 <= frame < g2_e1:
				return("G2")
			elif g2_e1 < frame <= g2_st2:
				return("G1")
			elif g2_st2 <= frame < g2_e2: #* This could be left for the next line and is equivalent to the final conditon
				return("G2")
			else:
				continue #* If the cell_cycle_df was able to pull an n+1 row, and it was stil not the final division, then there must be a row which follows, so the loop should continue
			#> End of the multirow conditions
		except IndexError: #* This is mostly to capture cells which only divided once, and to label the post final division as G1
			g2_st = g.iloc[0]["g2_st"]
			g2_e = g.iloc[0]["g2_e"]
			if frame < g2_st:
				return("G1") #* This is for the pre first division
			if g2_st <= frame < g2_e: #* There is no need to do the initial comparison because this will the last division
				return("G2")
			elif frame > g2_e:
				return("G1")
			else:
				continue
	return(None) #* The only way the function can get past this point is if there is not cell stage information at that frame for that cell

#. This should be implemented at the Quant_ALL stage of post quant. At this point, it is still based on the postion, which is the lowest level requirement, and still contains all the required information. It is not worth carrying all the needed rows through and not having access to the cell stage information
def cell_stage_assocation(position_barcode):
	try:
		cell_cycle_df = pd.read_parquet(f"{position_barcode}/", usecols=["track", "parent", "daughter", "tree", "generation", "age", "g2_st", "g2_e", "g2_dur", "div", "seg_st_track", "seg_st_daughter"]) #* Removed columns which do not have information based on input to TracX
		cell_cycle_df["Cell_Barcode"] = position_barcode + pd.Series("track").apply(lambda x: str(x).zfill(4)) #*Create the cell barcodes based on the tracks and the the positional position_barcode from function
		cell_cycle_df.set_index("Cell_Barcode", inplace= True) #* Set the table to be a lookup table to avoid having to the do a complex merge which would results in repetition
		localization_df = pd.read_parquet(f"Movement_treat_course_{position_barcode}.parquet", usecols= ["Cell_Barcode", "Relocalized", "Frame_x"]) #! Decide if this is the correct input
		localization_df["Cell_stage"] = localization_df.apply(lambda x: stage_from_frame(x["Cell_Barcode"], x["Frame_x"])) #* This is done like a lookup table
		localization_df.to_parquet("localization_cell_stage.parquet")
	except Exception as e:
		return(e, position_barcode) #* This will be returned to the user in the terminal

def cell_stage_manager(all_pos, pr):
	results, position_barcode = Parallel(n_jobs= pr, verbose= 50, prefer="treads")(delayed(cell_stage_assocation)(position) for position in all_pos) #* this is likely going to be io bound so set to prefer threads over the loky backend
	results = results.append(results) # This might work or it may fail
	return(position_barcode) # return the bacode upon completion

def stats_cell(Protein): #* This will have to be reffered to later on. It takes in a protein as input, which simplifies statistics
	Protein_df = pd.read_parquet(f"{Protein}.parquet")
	Protein_df.astype({"Cell_stage": 'category', "Relocalized": "category"})
	temp_for_stats = Protein_df.dropna(subset=['Cell_stage']) #* Remove cells for which we do not know the cell stage

	F_reloc, p_reloc = stats.f_oneway(temp_for_stats["Relocalized"], temp_for_stats["Cell_stage"]) #* Perform a one way ANOVA comparison of two categorical variables
	return(Protein,F_reloc,p_reloc) #*Based on this output, the program can be run as a vector and the results can be compiled into a dataframe

