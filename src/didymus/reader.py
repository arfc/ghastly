#import here
import numpy as np
from didymus import pebble as peb


#defining code-specific readers here:
class OpenmcReader():
	"""
	Class for reading Openmc pebble placement
	"""
	
	def __init__(self,coord_array):
		self.coord_array = coord_array
	
	def pebs(self,peb_rad,mat_ids,uniq_ids):
		peb_array = []
		for i, coord, in enumerate(self.coord_array):
			peb_array.append(peb.Pebble(coord,peb_rad,mat_ids[i],uniq_ids[i]))
			
		return peb_array
