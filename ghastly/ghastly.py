import numpy as np
import openmc
from jinja2 import Environment, PackageLoader
import ghastly
from ghastly import read_input
from ghastly import region

env = Environment(loader=PackageLoader('ghastly','templates'))

def fill_core(input_file):
    '''
    Fill core with pebbles using parameters from input file, and a combination
    of OpenMC packing and LAMMPS gran package's pour fix
    '''

    input_block = read_input.InputBlock(input_file)
    sim_block = input_block.create_obj()

    #should this also pack the outtake?  does pf/number of pebbles consider
    #only the main core when used in literature?

    #for now - require pf, fill number needed to get core_main+core_outtake
    #to that pf.  you can report pf if you want but lammps will also log it.
    rough_pack = []
    core_volume = 0
    pack_zones = sim_block.core_main | sim_block.core_outtake
    for element in pack_zones.values():
        core_volume += element.volume
        if type(element) == ghastly.core.CylCore:
            coords = pack_cyl(sim_block, element)
            rough_pack += coords
        else:
            pass

    n_pebbles = int((sim_block.pf*core_volume)/sim_block.pebble_volume)

    pebbles_left = n_pebbles - len(rough_pack)

    x_b, y_b, z_b = find_box_bounds(sim_block, pour = True)

    fake_dump_file(rough_pack, "rough-pack.txt", "ff ff ff",
                            x_b, y_b, z_b)

    #next: you have the input file read into input block, so now the next step
    #is feeding everything into jinja templates.

    #first, make the variable block file

    variables = input_block.lammps_var
    variables["r_pebble"] = sim_block.r_pebble
    variables["seed"] = sim_block.seed

    variables_template = env.get_template("variable_template.txt")
    variable_text = variables_template.render(variables = variables)
    variable_filename = "pour_variables.txt"
    with open(variable_filename, mode='w') as f:
        f.write(variable_text)

    #now regions:
    reg_files = []
    reg_names = []
    for element_name, element in pack_zones.items():
        reg_names.append(str(element_name))
        if type(element) == ghastly.core.CylCore:
            reg_template = env.get_template("cylcore_template.txt")
            reg_text = reg_template.render(region_name = element_name,
                                           x_c = element.x_c,
                                           y_c = element.y_c,
                                           r = element.r,
                                           z_min = element.z_min,
                                           z_max = element.z_max,
                                           open_bottom = element.open_bottom)
            reg_filename = str(element_name)+"_region.txt"
            reg_files.append(reg_filename)
            with open(reg_filename, mode='w') as f:
                f.write(reg_text)
        elif type(element) == ghastly.core.ConeCore:
            reg_template = env.get_template("conecore_template.txt")
            reg_text = reg_template.render(region_name = element_name,
                                           x_c = element.x_c,
                                           y_c = element.y_c,
                                           r_lower = element.r_lower,
                                           r_upper = element.r_upper,
                                           z_min = element.z_min,
                                           z_max = element.z_max,
                                           open_bottom = element.open_bottom)
            reg_filename = str(element_name)+"_region.txt"
            reg_files.append(reg_filenames)
            with open(reg_filename, mode='w') as f:
                f.write(reg_text)

        else:
            raise TypeError(str(element_name)+" is not a CylCore or ConeCore")

    #now the main file:
    
    match sim_block.down_flow:
        case True:
            flow_vector = "0 0 -1"
        case False:
            flow_vector = "0 0 1"
        case _:
            raise TypeError("down_flow should be true or false")
    
    #you need to determine the dimensions of the pour region based on whatever
    #core element in the core_main section is on top.  you will also
    #need to use that region's x_c and y_c (arguably you can average the box
    #bounds for that, but if you get any rounding it might throw things off,
    #so it's safer to copy the x_c/y_c of the region you are pouring into

    #for pour zone height - you have a bigger gap on top if you are pouring,
    #so the height of the pour region is just inside z_b.up to just above
    #the z_max of the z-most core_main section.

    #Godspeed.

    main_template = env.get_template("pour_main.txt")
    main_text = main_template.render(variable_filename = variable_filename,
                                     x_b = x_b,
                                     y_b = y_b,
                                     z_b = z_b,
                                     region_files = reg_files,
                                     n_regions = len(reg_files),
                                     region_names = reg_names,
                                     flow_vector = flow_vector,

                                     pebbles_left = pebbles_left)



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

        return list(coords)

def find_box_bounds(sim_block, pour = False):
        '''
        given lists of core component elements, determine the appropriate size
        of the lammps bounding box
        '''

        core_list = (sim_block.core_intake | 
                     sim_block.core_main | 
                     sim_block.core_outtake)
        x_list = []
        y_list = []
        z_list = []
        for element in core_list.values():
            z_list += [element.z_min, element.z_max]
            if type(element) == ghastly.core.CylCore:
                x_list += [-element.r, element.r]
                y_list += [-element.r, element.r]
            elif type(element) == ghastly.core.ConeCore:
                x_list += [-element.r_upper, element.r_upper, 
                           -element.r_lower, element.r_lower]
                y_list += [-element.r_upper, element.r_upper,
                           -element.r_lower, element.r_lower]

        #error!  you need the magnitude of the disance bt x_min and x_max to
        # always increase, but this will not always do that
        x_b = {"low": (min(x_list) - 0.05*min(x_list)), 
               "up": (max(x_list) + 0.05*max(x_list))}

        y_b = {"low": (min(y_list) - 0.05*min(y_list)), 
               "up": (max(y_list) + 0.05*max(y_list))}
        if pour = True:
            z_b = {"low": (min(z_list) - 0.05*min(z_list)), 
                    "up": (max(z_list) + 0.2*max(z_list))}
        else:
            z_b = {"low": (min(z_list) - 0.05*min(z_list)), 
                    "up": (max(z_list) + 0.05*max(z_list))}

        return x_b, y_b, z_b

def fake_dump_file(coords, dump_filename, bound_conds,
                       x_b, y_b, z_b):
        '''
        given coord array, create a "fake" dump file that can be imported into
        LAMMPS
        '''

        #environment = Environment(loader=PackageLoader('./ghastly/templates'))


        peb_list = [{"id":i, "x":v[0], "y":v[1], "z":v[2]} 
                    for i, v in enumerate(coords)]

        dump_template = env.get_template("dump_template.txt")
        dump_text = dump_template.render(n_rough_pack = len(coords),
                                         bound_conds = bound_conds,
                                         x_b = x_b,
                                         y_b = y_b,
                                         z_b = z_b,
                                         peb_list = peb_list)

        with open(dump_filename, mode='w') as f:
            f.write(dump_text)

