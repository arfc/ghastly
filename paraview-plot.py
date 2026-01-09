import os
from paraview.simple import *
import subprocess
import glob

#this is meant to be added to paraview as a python script and run within the
#paraview GUI.  Commented out lines will need to be adjusted by the user
#before running on their own system.

#csvpath = '~/pathto/directorywith/pebcoord/csvfiles'
#number of char in csv filenames.  For example files padded to 15 char have
# 19 characters (15 + the 4 in '.csv')
#sort_int = -19
cpath = os.path.expanduser(csvpath)
unsorted_csvnames = glob.glob(os.path.join(cpath, '*.csv'))
csv_fnames = sorted(unsorted_csvnames, key=lambda x:x[sort_int:])
csvreader = CSVReader(FileName=csv_fnames)
csvreader.HaveHeaders=0
points = TableToPoints(Input=csvreader)
points.XColumn = 'Field 1' 
points.YColumn = 'Field 2' 
points.ZColumn = 'Field 3'
points.UpdatePipeline()

