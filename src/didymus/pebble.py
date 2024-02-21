#imports

class Pebble:
	'''
	Contains data on a single pebble
	'''
	
	def __init__(self,coords,rad,mat_id,uniq_id,recirc=False):
		self.coords = coords
		self.rad = rad
		self.mat_id = mat_id
		self.uniq_id = uniq_id
		self.recirc = recirc
	
