import numpy as np
import pandas as pd
from PIL import Image
from scipy.io import loadmat

lookup_table = pd.read_paraquet("test.paraquet", index = True)
def Quant_addon(pos):
	pos_look = lookup_table.loc[pos].set_index("Frame")
	for m in frames:
		pos_frame_path = pos_look.loc[m]["Path"]
		mask = loadmat(pos_frame_path)
		#* this runs into problem because the location and sizes of the masks differ



def Add_myo():

	return(myo-KO, myo_Ka)