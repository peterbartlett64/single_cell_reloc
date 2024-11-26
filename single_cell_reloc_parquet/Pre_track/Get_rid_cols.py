#%%
#! This does not seem to work well as the program will try to read in the .mat version too which is harder to remove columns from

import os
import pandas as pd
import glob

for name in sorted(glob.glob('cells_*[0-9].txt')):
	temp = pd.read_csv(name, sep = '\t', usecols= ['cell.frame', 'cell.index', 'cell.center.x', 'cell.center.y', 'cell.majoraxis', 'cell.minoraxis', 'cell.orientation', 'cell.area', 'cell.volume', 'cell.perimeter', 'cell.eccentricity', 'cell.fractionOfGoodMembranePixels', 'cell.mem.area', 'cell.mem.volume', 'cell.nuc.radius', 'cell.nuc.area'])
	temp.to_csv(name, sep = '\t')
