class Core:
    '''
    Parent class for Core objects.
    '''
    def __init__(self, down_flow = True):
        '''
        Initializes a single instance of a Core object.  As
        this does not specify the shape of the core, the Core
        parent class should not be used directly.
        
        Parameters
        ----------
        down_flow : bool
            Whether axial flow in the core is upward or downward.
            True means flow is downward, False means it is upward.
        
        '''
        self.down_flow = down_flow
    
class CylCore(Core):
    '''
    Class for a cylindrically shaped core section, with its axis parallel
    to the z-axis.
    '''
    def __init__(self, x_c, y_c, r, z_max, z_min, regions, *args, **kwargs):
        '''
        Initializes a single instance of a CylCore object.  All dimensions
        should be in meters.
        
        Parameters
        ----------
        x_c : float
            Coordinate of the cylinder's center on the x-axis.
        y_c : float
            Coordinate of the cylinder's center on the y-axis.
        r : float
            Radius of the cylinder.
        z_max : float
            Z-coordinate of the cylinder's top.
        z_min : float
            Z-coordinate of the cylinder's bottom.
        regions : list
            List containing the region_id of each element within the given
            CylCore object.
        '''
        super().__init__(*args, **kwargs)
        self.x_c = x_c
        self.y_c = y_c
        self.r = r
        self.z_max = z_max
        self.z_min = z_min
        self.height = abs(z_min) + abs(z_max)
        self.regions = regions

class ConeCore(Core):
    '''
    Class for a right truncated cone, such as those found in discharge chutes.
    The central axis is parallel to the z-axis.  
    '''
    def __init__(self, x_c, y_c, r_upper, r_lower, z_max, z_min, 
                 regions, *args, **kwargs):
        '''
        Intializes a single instance of a ConeCore object.  All distances are
        in meters.
        Parameters
        ----------
        x_c : float
            Coordinate of the cone's center on the x_axis
        y_c : float
            Coordinate of the cone's center on the y_axis
        r_upper: float
            Radius of the top of the cone.  Does not necessarily correspond to 
            the largest radius.
        r_lower : float
            Radius of the bottom of the cone.  Does not necessarily correspond
            to the smallest radius.
        z_max : float
            Z_coordinate of the cone's top
        z_min : float
            Z_coordinate of the cone's bottom
        regions : list
            List containing the region_id of each element within the given
            ConeCore object.
        '''
        super().__init__(*args, **kwargs)
        self.x_c = x_c
        self.y_c = y_c
        self.r_upper = r_upper
        self.r_lower = r_lower
        self.z_max = z_max
        self.z_min = z_min
        self.height = abs(z_min)+abs(z_max)
        self.regions = regions
