import numpy as np


class Core:
    '''
    Parent class for Core objects.  Because the Core class doesn't define
    a specific geometry, the Core class should not be used directly.
    '''

    def __init__(self, x_c, y_c, z_max, z_min):
        '''
        Initializes a single instance of a Core object.

        Parameters
        ----------
        x_c : float
            Coordinate of the cone's center on the x_axis.
        y_c : float
            Coordinate of the cone's center on the y_axis.
        z_max : float
            Z-coordinate of the cylinder's top.
        z_min : float
            Z-coordinate of the cylinder's bottom.

        '''
        self.x_c = x_c
        self.y_c = y_c
        self.z_max = z_max
        self.z_min = z_min
        self.h = abs(z_min) + abs(z_max)


class CylCore(Core):
    '''
    Class for a cylindrically shaped core section, with its axis parallel
    to the z-axis.
    '''

    def __init__(self, r, regions, *args, **kwargs):
        '''
        Initializes a single instance of a CylCore object.  All dimensions
        should be in meters.

        Parameters
        ----------
        r : float
            Radius of the cylinder [m].
        regions : list
            List containing the region_id of each element within the given
            CylCore object.
        '''
        super().__init__(*args, **kwargs)
        self.r = r
        self.volume = np.pi * (r**2) * self.h
        self.regions = regions


class AnnularCore(Core):
    '''
    Class for an annularly shaped core region, with its axis parallel to the
    z-axis.
    '''

    def __init__(self, r_outer, r_inner, regions, *args, **kwargs):
        '''
        Initializes an AnnularCore object.  Note that unlike an AnnularRegion,
        this assumes it is a full annular shell, not a sector.  Dimensions are
        in meters.

        Parameters
        ----------
        r_outer : float
            Outer radius of the annulus [m].
        r_inner : float
            Inner radius of the annulus [m].
        regions : list
            List containing the reg_id of each region within the core element.
        '''
        super().__init__(*args, **kwargs)
        self.r_outer = r_outer
        self.r_inner = r_inner
        self.volume = (np.pi * (r_outer**2 - r_inner**2) * self.h)
        self.regions = regions


class ConeCore(Core):
    '''
    Class for a right truncated cone, such as those found in discharge chutes.
    The central axis is parallel to the z-axis.
    '''

    def __init__(self, r_upper, r_lower, regions, *args, **kwargs):
        '''
        Intializes a single instance of a ConeCore object.  All distances are
        in meters.

        Parameters
        ----------
        r_upper: float
            Radius of the top of the cone [m].  Does not necessarily correspond
            to the largest radius.
        r_lower : float
            Radius of the bottom of the cone [m].  Does not necessarily
            correspond to the smallest radius.
        regions : list
            List containing the region_id of each element within the given
            ConeCore object.
        '''
        super().__init__(*args, **kwargs)
        self.r_upper = r_upper
        self.r_lower = r_lower
        self.volume = ((1 / 3) * np.pi * self.h
                       * (r_upper**2 + r_lower**2 + r_upper * r_lower))
        self.regions = regions


class AnnConeCore(Core):
    '''
    Class for an annular right truncated cone, with its centerline parallel
    to the z-axis
    '''

    def __init__(self, r_out_up, r_in_up, r_out_low, r_in_low,
                 *regions, *args, **kwargs):
        '''
        Initializes an AnnConeCore object.  All distances should be in meters.

        Parameters
        ----------
        r_out_up : float
            The outer radius at the upper part of the annular cone [m].
        r_in_up : float
            The inner radius at the upper part of the annular cone [m].
        r_out_low : float
            The outer radius at the lower part of the annular cone [m].
        r_in_low : float
            The inner radius at the lower part of the annular cone [m].
        regions : list
            List containing the reg_id of each element within the core element.
        '''
        super().__init__(*args, **kwargs)
        self.r_out_up = r_out_up
        self.r_in_up = r_in_up
        self.r_out_low = r_out_low
        self.r_in_low = r_in_low
        self.volume = ((1 / 3) * np.pi * self.h * (
            (r_out_up**2 + r_out_low**2 + r_out_up * r_out_low)
            - (r_in_up**2 + r_in_low**2 + r_in_up * r_in_low)))
