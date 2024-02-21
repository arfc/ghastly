# import here
import numpy as np
from didymus import pebble as peb


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
        coord_array : array
            Array of center coordinates, in x, y, and z, from
            the pack_spheres() function of OpenMC.

        '''
        self.coord_array = coord_array

    def pebs(self, peb_rad, mat_ids, uniq_ids):
        '''
        Creates an array of didymus Pebble objects, using
        the central coordinates in coord_array.

        Parameters
        ----------
        peb_rad : float
            Radius of a single pebble, with units matching those
            used to create center coordinates.  Assumes all
            pebbles are the same size.
        mat_ids : array of int or str
            Array containing the mat_ids associated with each pebble
            center.  While multiple pebbles can share the same mat_id,
            the length of the mat_ids array should match coord_array.
        uniq_ids : array of int
            Array containing the uniq_ids associated with each pebble
            center.  Each pebble should have a distinct uniq_id.

        Returns
        -------
        peb_array : array of Pebble objects
            Array containing a unique Pebble object for each
            coordinate provided in coord_array.

        '''
        peb_array = []
        for i, coord, in enumerate(self.coord_array):
            peb_array.append(
                peb.Pebble(
                    coord,
                    peb_rad,
                    mat_ids[i],
                    uniq_ids[i]))

        return peb_array
