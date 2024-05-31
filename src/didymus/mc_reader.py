# import here
import numpy as np
from didymus.pebble import Pebble


# defining code-specific readers here:
class OpenmcReader():
    """
    Class for reading Openmc pebble coordinates from the
    OpenMC pack_spheres() function into and arrary of
    distinct Pebble objects.
    """

    def __init__(self, coord_array):
        '''
        Initializes the OpenmcReader.

        Parameters
        ----------
        coord_array : numpy.ndarray
            A 2-D numpy array with shape (N, 3) of cartesian coordinates 
            for the centroid of a sphere, obtained from
            :function:`openmc.pack_spheres()`, where N is the number
            of spheres.

        '''
        self.coord_array = coord_array

    def generate_pebbles(self, pebble_radius, mat_ids, pebble_ids):
        '''
        Creates an array of didymus Pebble objects, using
        the central coordinates in coord_array.

        Parameters
        ----------
        pebble_radius : float
            Radius of a single pebble, with units matching those
            used to create center coordinates.  Assumes all
            pebbles are the same size.
        mat_ids : array of int or str
            Array containing the mat_ids associated with each pebble
            center.  While multiple pebbles can share the same mat_id,
            the length of the mat_ids array should match coord_array.
        pebble_ids : array of int
            Array containing the uniq_ids associated with each pebble
            center.  Each pebble should have a distinct pebble_id.

        Returns
        -------
        pebble_array : List of :class:`didymus.Pebble` objects
            A list of unique :class:`didymus.Pebble` objects for each
            set of centroid coordinates provided in coord_array.

        '''
        peb_array = []
        for i, coord, in enumerate(self.coord_array):
            peb_array.append(
                Pebble(
                    coord,
                    peb_rad,
                    mat_ids[i],
                    uniq_ids[i]))

        return peb_array
