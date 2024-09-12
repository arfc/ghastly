# imports

class Pebble:
    '''
    Class representing a single pebble. Contains information
    necessary to determine the pebble location within the core
    and track its material composition and unique id number.
    '''

    def __init__(self, coords, r, reg_id, pass_num, l_type, pebble_id):
        '''
        Initializes a single instance of a Pebble object.

        Parameters
        ----------
        coords : array
            Array containing the x, y, and z coordinates of the pebble.
        radius : float
            Radius of the pebble, with the same units as core measurements.
        reg_id : str
            Region identifier corresponding to the region the pebble is located in
            based upon its center coordinates.
        pass_num : int
            Number of passes the pebble has completed.  For example, a fresh
            pebble's pass_num would be 0.
        l_type : int
            The 'type' number of the pebble, as defined by lammps.  Each
            l_type corresponds to a specific material.
        pebble_id : int
            User-defined integer for identifying a specific pebble.
            Each pebble_id should be distinct.

        '''
        self.coords = coords
        self.r = r
        self.reg_id = reg_id
        self.pass_num = pass_num
        self.l_type = l_type
        self.pebble_id = pebble_id
