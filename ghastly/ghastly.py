import numpy as np
import openmc
from jinja2 import Environment, FileSystemLoader
from ghastly import read_input
from ghastly import simulation
from ghastly import core
from ghastly import region

def fill_core(input_file):
    '''
    Fill core with pebbles using parameters from input file, and a combination
    of OpenMC packing and LAMMPS gran package's pour fix
    '''

    inputblock = read_input.InputBlock(input_file)
    simblock = inputblock.create_obj()

    #should this also pack the outtake?  does pf/number of pebbles consider
    #only the main core when used in literature?

    #for now - require pf, fill number needed to get core_main+core_outtake
    #to that pf.  you can report pf if you want but lammps will also log it.
    rough_pack = []
    core_volume = 0
    pack_zones = simblock.core_main | simblock.core_outtake
    for element in pack_zones.values():
        core_volume += element.volume
        if type(element) == ghastly.core.CylCore:
            coords = self.pack_cyl(simblock, element)
            rough_pack += coords
        else:
            pass

    n_pebbles = int((self.pf*core_volume)/self.pebble_volume)

    pebbles_left = n_pebbles - len(rough_pack)

    x_b, y_b, z_b = self.find_box_bounds(simblock)

    self.fake_dump_file(rough_pack, "rough-pack.txt", "ff ff ff",
                            x_b, y_b, z_b)



def pack_cyl(simblock, element):
        '''
        packs a cylindrical core element
        '''
        sides = openmc.ZCylinder(r=element.r)
        top = openmc.ZPlane(z0=element.z_max)
        bottom = openmc.ZPlane(z0=element.z_min)
        region_bounds = -sides & -top & +bottom
        
        coords = openmc.model.pack_spheres(simblock.r_pebble, 
                                           region = region_bounds,
                                           pf = simblock.pf,
                                           contraction_rate = simblock.k_rate)

        return coords

def find_box_bounds(simblock):
        '''
        given lists of core component elements, determine the appropriate size
        of the lammps bounding box
        '''

        core_list = (simblock.core_intake | 
                     simblock.core_main | 
                     simblock.core_outtake)
        x_list = []
        y_list = []
        z_list = []
        for element in core_list.values:
            z_list += [element.z_min, element.z_max]
            if type(element) == ghastly.core.CylCore:
                x_list += [-element.r, element.r]
                y_list += [-element.r, element.r]
            elif type(element) == ghastly.core.ConeCore:
                x_list += [-element.r_upper, element.r_upper, 
                           -element.r_lower, element.r_lower]
                y_list += [-element.r_upper, element.r_upper,
                           -element.r_lower, element.r_lower]
        x_b = {"low": (min(x_list) - 0.05*min(x_list)), 
               "up": (max(x_list) + 0.05*max(x_list))}

        y_b = {"low": (min(y_list) - 0.05*min(y_list)), 
               "up": (max(y_list) + 0.05*max(y_list))}

        z_b = {"low": (min(z_list) - 0.05*min(z_list)), 
               "up": (max(z_list) + 0.05*max(z_list))}

        return x_b, y_b, z_b

def fake_dump_file(coords, dump_filename, bound_conds,
                       x_b, y_b, z_b):
        '''
        given coord array, create a "fake" dump file that can be imported into
        LAMMPS
        '''

        peb_list = [{"id":i, "x":v[0], "y":v[1], "z":v[2]} 
                    for i, v in enumerate(coords)]

        dump_template = environment.get_template("dump_template.txt")
        dump_text = dump_template.render(n_rough_atoms = len(coords),
                                         bound_conds = bound_conds,
                                         x_b = x_b,
                                         y_b = y_b,
                                         z_b = z_b,
                                         peb_list = peb_list)

        with open(dump_filename, mode='w') as f:
            f.write(dump_text)

