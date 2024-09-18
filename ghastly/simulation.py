import numpy as np
import openmc
from ghastly import core

class Sim:
    '''
    Class for containing simulation-wide parameters and methods.
    '''

    def __init__(self, r_pebble, t_final,pf, down_flow=True):
        '''
        Initializes the Sim class.

        Parameters
        ----------
        r_pebble : float
            Radius of the pebbles in the simulation [m].
        t_final : float
            Total reactor-time that is being simulated [s].
        pf : float
            Target packing fraction in core region.  Determines total number
            of pebbles, but packing fraction may not match this value in
            all areas due to settling.
        down_flow : bool
            The direction of flow in the reactor core.  True means it is
            flowing downward, False means it is flowing upward.  Default is
            down.
        '''
        self.r_pebble = r_pebble
        self.pebble_volume = (4/3)*np.pi*(r_pebble**3)
        self.t_final = t_final
        self.pf = pf
        self.down_flow = down_flow

    def run_sim(self):
        '''
        Run a ghastly simulation
        '''

        pass

    def pack_cyl(self, element):
        '''
        packs a cylindrical core element
        '''
        #make the openmc region based on the core
        

        #pack region to pf

        #return the coords you get
        
    
    def pack_core(self, core_elements):
        '''
        initial pack for core.  openmc packs cylindrical/annular regions,
        then passes to LAMMPS to fill the rest needed.
        '''
        rough_pack = []
        for element in core_elements:
            if type(element) == core.CylCore:
                pass
                #rough_pack = self.pack_cyl(element)






