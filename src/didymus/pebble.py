# imports

class Pebble:
    '''
    Class for a single pebble object.  Contains information
    necessary to determine the pebble location within the core
    and track its material composition and unique id number.
    '''

    def __init__(self, coords, rad, mat_id, uniq_id, recirc=False):
        '''
        Initializes a single instance of a Pebble object.

        Parameters
		----------

		coords : array
            Array containing the x, y, and z coordinates of the pebble.
        rad : float
            Radius of the pebble, with the same units as core measurements.
        mat_id : int or str
            User-defined integer or string for the material in the pebble,
            which multiple pebble objects can share.
        uniq_id : int
            User-defined integer for identifying a specific pebble.
            Each uniq_id should be distinct.
        recirc : bool
            Boolean for determining if a given pebble should recirculate,
            for use in multi-pass cycles.  Defaults to False.

        '''
        self.coords = coords
        self.rad = rad
        self.mat_id = mat_id
        self.uniq_id = uniq_id
        self.recirc = recirc
