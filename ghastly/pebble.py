# imports

class Pebble:
    '''
    Class representing a single pebble. Contains information
    necessary to determine the pebble location within the core
    and track its material composition and unique id number.
    '''

    def __init__(self, uid, coords, velocity, zone, layer, 
                 pass_num, recirc, l_type, history = []):
        '''
        Initializes a single instance of a Pebble object.

        Parameters
        ----------
        coords : array
            Array containing the x, y, and z coordinates of the pebble.
        radius : float
            Radius of the pebble, with the same units as core measurements.
        reg_id : str
            Region identifier corresponding to the region the pebble is located
            in based upon its center coordinates.
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
        self.uid = int(uid)
        self.coords = coords
        self.velocity = velocity
        self.zone = int(zone)
        self.layer = int(layer)
        self.pass_num = int(pass_num)
        self.recirc = int(recirc)
        self.l_type = int(l_type)
        self.history = history

