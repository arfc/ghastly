import numpy as np

class Core:
    '''
    Class for a single Core object, which will be referenced
    in moving Pebble objects.
    '''
    def __init__(self, pebble_radius, origin=np.zeros(3), down_flow = True):
        '''
        Initializes a single instance of a Core object.  As
        this does not specify the shape of the core, the Core
        parent class should not be used directly.
        
        Parameters
        ----------
        pebble_radius : float
            Radius of pebbles in core.
        origin : numpy array
            A numpy array consisting of 3 elements: the x, y,
            and z coordinates of the core's origin.  Default
            is centered at (0.0,0.0,0.0).  Units should be meters.
        down_flow : bool
            Whether axial flow in the core is upward or downward.
            True means flow is downward, False means it is upward.
        
        '''
        self.pebble_radius = pebble_radius
        self.origin = origin
        self.down_flow = down_flow
    
class CylCore(Core):
    '''
    Class for a cylindrically shaped core section, with its axis parallel
    to the z-axis.
    '''
    def __init__(self, radius, height, 
                 regions *args, **kwargs):
        '''
        Initializes a single instance of a CylCore object.  All dimensions
        should be in meters.
        
        Parameters
        ----------
        radius : float
            Radius of the core.
        height : float
            Height of the core.
        regions : list
            List containing the region_id of each element within the given
            CylCore object.
        '''
        super().__init__(*args, **kwargs)
        self.core_radius = core_radius
        self.core_height = core_height
        self.regions = regions

class ConeCore(Core):
    '''
    Class for a right truncated cone, such as those found in discharge chutes.
    The central axis is parallel to the z-axis.  
    '''
    def __init__(self, upper_radius, lower_radius, height, regions):
        '''
        Intializes a single instance of a ConeCore object.  All distances are
        in meters.
        Parameters
        ----------
        upper_radius : float
            Radius of the top of the cone.  Does not necessarily correspond to 
            the largest radius.
        lower_radius : float
            Radius of the bottom of the cone.  Does not necessarily correspond
            to the smallest radius.
        height : float
            Height of the cone, along central axis.
        regions : list
            List containing the region_id of each element within the given
            ConeCore object.
        '''
        super().__init__(*args, **kwargs)
        self.upper_radius = upper_radius
        self.lower_radius = lower_radius
        self.height = height
        self.regions = regions


        
        
        
