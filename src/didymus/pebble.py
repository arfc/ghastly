# imports

class Pebble:
    '''
    Class for a single pebble object. Contains information
    necessary to determine the pebble location within the core
    and track its material composition and unique id number.
    '''

    def __init__(self, coords, radius, mat_id, pebble_id, recirc=False):
        '''
        Initializes a single instance of a Pebble object.

        Parameters
        ----------

        coords : array
            Array containing the x, y, and z coordinates of the pebble.
        radius : float
            Radius of the pebble, with the same units as core measurements.
        mat_id : int or str
            User-defined integer or string for the material in the pebble,
            which multiple pebble objects can share.
        pebble_id : int
            User-defined integer for identifying a specific pebble.
            Each uniq_id should be distinct.
        recirc : bool
            Boolean for determining if a given pebble should recirculate,
            for use in multi-pass cycles.  Defaults to False.

        '''
        self.coords = coords
        self.radius = radius
        self.mat_id = mat_id
        self.pebble_id = pebble_id
        self.recirc = recirc
