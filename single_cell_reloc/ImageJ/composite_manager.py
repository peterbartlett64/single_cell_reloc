import pandas as pd
import os
import imagej


def composite_manger(position):
	os.chdir(Global_Variables["microfluidics_results"])
	imgIndex =   pd.read_csv('imgIndex.csv')
