# imports

class Pebble:
    '''
    Class representing a single pebble. Contains information
    necessary to determine the pebble location within the core
    and track its material composition and unique id number.
    '''

    def __init__(self, uid, coords, velocity, zone, layer, 
                 pass_num, recirc, history = []):
        '''
        Initializes a single instance of a Pebble object.

        Parameters
        ----------
        uid : int
            Unique integer ID of a pebble
        coords : array
            Array containing the x, y, and z coordinates of the pebble.
        velocity : array
            Array containing the x, y, and z coordinates of the pebble.        
        zone : int
            Integer value corresponding to the label of the radial zone the
            pebble began transit in
        layer : int
            Integer value corresponding to the label of the axial layer the
            pebble began transit in
        pass_num : int
            Number of passes the pebble has completed.  A fresh pebble's 
            pass_num would be 0.
        recirc : int
            Number of times the pebble has recirculated since transit began.
        history : list
            A list containing the radial zone ID the pebble was located in for
            the majority of previous passes.  The first element is the first
            pass, the second is the second pass, and so on.

        '''
        self.uid = int(uid)
        self.coords = coords
        self.velocity = velocity
        self.zone = int(zone)
        self.layer = int(layer)
        self.pass_num = int(pass_num)
        self.recirc = int(recirc)
        self.history = history

