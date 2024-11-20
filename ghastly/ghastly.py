import numpy as np
import openmc
import ghastly
from jinja2 import Environment, PackageLoader
from ghastly import read_input
from ghastly import region
from lammps import lammps

env = Environment(loader=PackageLoader('ghastly', 'templates'))


def fill_core(input_file, rough_pf):
    '''
    This function uses OpenMC's model.pack_spheres function to "rough" pack
    the full core region defined by the Ghastly input file, then uses a series
    of Jinja templates to create a LAMMPS input file that will run a physics
    simulation, modeling not only the settling of the roughly packed pebbles
    due to gravtiy, but also handle the addtion of any extra pebbles needed to
    reach the target number of pebbles needed.

    Parameters
    ----------
    input_file : str
        Name of ghastly input file, which should be in JSON format, following
        the layout and parameter naming convention given in the example inputs.
    rough_pf : float
        Packing fraction used in OpenMC packing function.  Can match the target
        packing fraction provided in the input file, but users should take care 
        supplying a rough_pf that is greater than the input's target pf.
    '''

    input_block = read_input.InputBlock(input_file)
    sim_block = input_block.create_obj()

    rough_pack = []
    pack_zones = sim_block.core_main | sim_block.core_outtake
    core_volume = sum([i.volume for i in pack_zones.values()])
    for element in pack_zones.values():
        if type(element) == ghastly.core.CylCore:
            coords = pack_cyl(sim_block, element, rough_pf)
            rough_pack += coords
        else:
            pass

    n_pebbles = int((sim_block.pf*core_volume)/sim_block.pebble_volume)
    print("Equivalent number of pebbles is "+str(n_pebbles))

    pebbles_left = n_pebbles - len(rough_pack)

    if pebbles_left < 0:
        raise ValueError('''Negative pebbles left - core overfilled.  Reduce
                         rough_pf and try again.''')

    x_b, y_b, z_b = find_box_bounds(sim_block, pour=True)

    fake_dump_file(rough_pack, "rough-pack.txt", "ff ff ff",
                   x_b, y_b, z_b)

    # next: you have the input file read into input block, so now the next step
    # is feeding everything into jinja templates.

    # first, make the variable block file

    variables = input_block.lammps_var
    variables["r_pebble"] = sim_block.r_pebble
    variables["seed"] = sim_block.seed

    variables_template = env.get_template("variable_template.txt")
    variable_text = variables_template.render(variables=variables)
    variable_filename = "pour_variables.txt"
    with open(variable_filename, mode='w') as f:
        f.write(variable_text)

    # now regions:
    reg_files = []
    reg_names = []
    for element_name, element in pack_zones.items():
        reg_names.append(str(element_name))
        if type(element) == ghastly.core.CylCore:
            reg_template = env.get_template("cylcore_template.txt")
            reg_text = reg_template.render(region_name=element_name,
                                           x_c=element.x_c,
                                           y_c=element.y_c,
                                           r=element.r,
                                           z_min=element.z_min,
                                           z_max=element.z_max,
                                           open_bottom=element.open_bottom)
            reg_filename = str(element_name)+"_region.txt"
            reg_files.append(reg_filename)
            with open(reg_filename, mode='w') as f:
                f.write(reg_text)
        elif type(element) == ghastly.core.ConeCore:
            reg_template = env.get_template("conecore_template.txt")
            reg_text = reg_template.render(region_name=element_name,
                                           x_c=element.x_c,
                                           y_c=element.y_c,
                                           r_lower=element.r_lower,
                                           r_upper=element.r_upper,
                                           z_min=element.z_min,
                                           z_max=element.z_max,
                                           open_bottom=element.open_bottom)
            reg_filename = str(element_name)+"_region.txt"
            reg_files.append(reg_filename)
            with open(reg_filename, mode='w') as f:
                f.write(reg_text)

        else:
            raise TypeError(str(element_name)+" is not a CylCore or ConeCore")

    # now the main file:

    match sim_block.down_flow:
        case True:
            flow_vector = "0 0 -1"
        case False:
            flow_vector = "0 0 1"
        case _:
            raise TypeError("down_flow should be True or False.")
    main_core_z_max = max([(key, element.z_max)
                           for key, element in sim_block.core_main.items()])
    main_top = sim_block.core_main[main_core_z_max[0]]

    x_c_pour = main_top.x_c
    y_c_pour = main_top.y_c
    z_max_pour = 0.95*z_b["up"]
    z_min_pour = 1.05*main_core_z_max[1]
    if type(main_top) == ghastly.core.CylCore:
        r_pour = 0.75*main_top.r
    elif type(main_top) == ghastly.core.ConeCore:
        r_pour = 0.75*main_top.r_upper

    main_template = env.get_template("pour_main.txt")
    main_text = main_template.render(variable_filename=variable_filename,
                                     x_b=x_b,
                                     y_b=y_b,
                                     z_b=z_b,
                                     region_files=reg_files,
                                     n_regions=len(reg_files),
                                     region_names=reg_names,
                                     flow_vector=flow_vector,
                                     x_c_pour=x_c_pour,
                                     y_c_pour=y_c_pour,
                                     r_pour=r_pour,
                                     z_min_pour=z_min_pour,
                                     z_max_pour=z_max_pour,
                                     pebbles_left=pebbles_left)

    main_filename = "pour_main_input.txt"
    with open(main_filename, mode='w') as f:
        f.write(main_text)

    lmp = lammps()
    lmp.file(main_filename)


def pack_cyl(sim_block, element, rough_pf):
    '''
    Given a cylindrical core object, this function creates a corresponding
    OpenMC region, and the the pack_spheres function to generate a list of
    non-overlapping pebble centroids that fit within that region.

    Parameters
    ----------
    sim_block : ghastly Sim object
        Sim class object that contains simulation-specific simulation,
        read from input_file.
    element : ghastly CylCore object
        CylCore object for a portion of the core model.
    rough_pf : float
        Packing fraction that OpenMC will pack this cylinder to.  Not the
        same as the packing fraction in the input.

    '''
    sides = openmc.ZCylinder(x0=element.x_c,
                             y0=element.y_c,
                             r=element.r)
    top = openmc.ZPlane(z0=element.z_max)
    bottom = openmc.ZPlane(z0=element.z_min)
    region_bounds = -sides & -top & +bottom

    coords = openmc.model.pack_spheres(sim_block.r_pebble,
                                       region=region_bounds,
                                       pf=rough_pf,
                                       contraction_rate=sim_block.k_rate)

    return list(coords)


def find_box_bounds(sim_block, pour=False):
    '''
    Using the core_zones contained in the sim_block, this function
    uses the dimensions of all core elements in the model to determine the
    appropriate size of the bounding box used in LAMMPS simulations.
    Note that when used to determine the bounding box for a simulation
    using the pour LAMMPS fix, the bounding box is extended in the
    positive z direction, to allow for a temporary insertion region to
    be added.

    Parameters
    ----------
    sim_block : Sim object
        Sim object created from parameters in the input file for ghastly.
    pour : bool
        Default False.  Whether or not the LAMMPS simulation the bounds
        are being generated for is using the pour fix in LAMMPS

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
            x_list += [(element.x_c - element.r),
                       (element.x_c + element.r)]
            y_list += [(element.y_c - element.r),
                       (element.y_c + element.r)]
        elif type(element) == ghastly.core.ConeCore:
            x_list += [(element.x_c - element.r_upper),
                       (element.x_c + element.r_upper),
                       (element.x_c - element.r_lower),
                       (element.x_c + element.r_lower)]
            y_list += [(element.y_c - element.r_upper),
                       (element.y_c + element.r_upper),
                       (element.y_c - element.r_lower),
                       (element.y_c + element.r_lower)]

    match pour:
        case True:
            f = 1.05
            f_zup = 1.2
        case _:
            f = 1.05
            f_zup = 1.05
    x_b = {"low": ((1-f)*element.x_c + f*min(x_list)),
           "up": ((1-f)*element.x_c + f*max(x_list))}
    y_b = {"low": ((1-f)*element.y_c + f*min(y_list)),
           "up": ((1-f)*element.y_c + f*max(y_list))}
    z_b = {"low": ((1-f)*0.5*(max(z_list)+min(z_list)) + f*min(z_list)),
           "up": ((1-f)*0.5*(max(z_list)+min(z_list)) + f_zup*max(z_list))}

    return x_b, y_b, z_b


def fake_dump_file(coords, dump_filename, bound_conds,
                   x_b, y_b, z_b):
    '''
    Using the coordinate array and simulation boundary conditions and
    dimentions, this function uses a jinja template
    to create a LAMMPS dumpfile that be read into LAMMPS.

    Parameters
    ----------
    coords : list
        List of pebble centroid coordinates
    dump_filename : str
        Desired name of LAMMPS dumpfile created.
    bound_conds : str
        String providing the exact surface boundary conditions to be used
        in the LAMMPS bounding box.  See LAMMPS documentation for more
        information on bounding box conditions.
    x_b : dict
        Dictionary with key: value pairs giving the upper and lower
        bounds in the x-direction for the bounding box.
    y_b : dict
        As x_b, but for the y direction.
    z_b : dict
        as x_b but for the z direction.

    '''

    peb_list = [{"id": i, "x": v[0], "y": v[1], "z": v[2]}
                for i, v in enumerate(coords)]

    dump_template = env.get_template("dump_template.txt")
    dump_text = dump_template.render(n_rough_pack=len(coords),
                                     bound_conds=bound_conds,
                                     x_b=x_b,
                                     y_b=y_b,
                                     z_b=z_b,
                                     peb_list=peb_list)

    with open(dump_filename, mode='w') as f:
        f.write(dump_text)
