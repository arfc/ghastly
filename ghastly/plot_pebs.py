import numpy as np
import openmc
import ghastly
from ghastly import read_input
from ghastly import pebble
from ghastly import velocity
from jinja2 import Environment, PackageLoader
from os import path
from paraview.simple import *
import subprocess
import glob
import string
import matplotlib as mpl

rng = np.random.default_rng()

def plot_bed_xs(filepath, velpath, coordpath, recircpath, inputfile, 
                vel_plot=False, color_crit=[], start = 0, end = None,
                ):
    '''
    filepath : str
        file path leading to ghastly output /coords directory.  use ~.  
        example: '~/sims/model1/coords'
    inputfile: ghastly input file, for things like r_peb, reactor geometry, etc
    '''
    fpath = path.expanduser(filepath)
    inp_path = path.expanduser(inputfile)
    inp_block = read_input.InputBlock(inp_path)
    sim_block = inp_block.create_obj()

    unsorted_names = glob.glob(path.join(fpath, "*.bin"))
    filenames = sorted(unsorted_names, key=lambda x:x[-19:])
    
    peb_mats = create_peb_mats(color_crit)
    all_peb_lists = []
    for i, file in enumerate(filenames[start:end]):
        if i%1 == 0:
            subprocess.run(['binary2txt', file])
            txtfile = glob.glob(path.join(fpath, "*.txt"))[0]
            data = np.loadtxt(txtfile, delimiter=' ', skiprows=9)
            pebbles = []
            for d in data:
                pebbles.append(pebble.Pebble(uid=d[0],
                                            coords=100*d[-3:],
                                            velocity=0.0,
                                            zone=max(d[2], 1),
                                            layer=d[3],
                                            pass_num=d[4],
                                            recirc=d[5],
                                            l_type=d[1]))
            peb_list = create_plot_pebs(pebbles, peb_mats,
                                    sim_block, color_crit)
            all_peb_lists.append(peb_list)

            #cleanup
            subprocess.run(['rm', txtfile])
    return sim_block, peb_mats, all_peb_lists

def create_peb_mats(color_crit):
    '''
    create openmc materials based on color criterion
    '''
    id_counter = 10
    peb_mats = []
    for crit in color_crit:
        mat_crit_list = []
        for i in crit[1]:
            name = crit[0]+str(i)
            temp_mat = openmc.Material(name=name, material_id=id_counter)
            temp_mat.set_density('g/cm3', 1.0)
            temp_mat.add_element('U', 1.00)
            mat_crit_list.append(temp_mat)
            id_counter += 10
        peb_mats.append(mat_crit_list)
    return peb_mats


def create_plot_pebs(pebbles, peb_mats, sim_block, color_crit):
    '''
    create the geometries and mats and suchlike needed to plot pebbles
    '''
    
    peb_or = 100*sim_block.r_pebble
    peb_out = openmc.Sphere(r=peb_or)
    half_peb = openmc.XPlane()
    left_reg = -peb_out & -half_peb
    right_reg = -peb_out & +half_peb
    
    peb_dict = {}
    check_keys = []
    for i, v1 in enumerate(color_crit[0][1]):
        for j, v2 in enumerate(color_crit[1][1]):
            #check_keys.append((v1, v2))
            key = (v1, v2)
            peb_dict[key] = {}
            peb_dict[key]['coords'] = []

            left_cell = openmc.Cell(region = left_reg, 
                                    fill = peb_mats[0][i])
            right_cell = openmc.Cell(region = right_reg, 
                                     fill = peb_mats[1][j])
            peb_univ = openmc.Universe(cells=[left_cell, right_cell])
            peb_dict[key]['univ'] = peb_univ

    crit1 = color_crit[0][0].casefold()
    crit2 = color_crit[1][0].casefold()
    for peb in pebbles:
        match crit1:
            case 'zone':
                v1 = peb.zone
            case 'recirc':
                v1 = peb.recirc
            case 'pass_num':
                v1 = peb.pass_num
            case 'layer':
                v1 = peb.layer
        match crit2:
            case 'zone':
                v2 = peb.zone
            case 'recirc':
                v2 = peb.recirc
            case 'pass_num':
                v2 = peb.pass_num
            case 'layer':
                v2 = peb.layer
        peb_dict[(v1, v2)]['coords'].append(peb.coords)
    
    peb_list = []
    for crit_key, group in peb_dict.items():
        peb_list += [openmc.model.TRISO(peb_or, group['univ'], coord) 
                     for coord in group['coords']]

    return peb_list


def plot_vel_xs(directory, velpath, coordpath, recircpath, inputfile,
                sort_int=-19, n_skip=1, delimiter=' ', skiprows=9,
                n_recirc=2500000, n_dump=372366, dt=2.6855e-07,
                recirc_hz=0.014, start = 0, end = None, colormap='managua'):
    '''
    only prints a frame for every 10 data points - with 10 points per second,
    that's 1 frame per second
    filepath : str
        file path leading to ghastly output /coords directory.  use ~.  
        example: '~/sims/model1/coords'
    inputfile: ghastly input file, for things like r_peb, reactor geometry, etc
    '''
    inp_path = path.expanduser(inputfile)
    inp_block = read_input.InputBlock(inp_path)
    sim_block = inp_block.create_obj()
    peb_r = 100*sim_block.r_pebble

    mod = openmc.Material(name='moderator')
    mod.set_density('g/cm3', 1.0)
    mod.add_element('C', 1.0)

    cool = openmc.Material(name='coolant')
    cool.set_density('g/cm3',1.0)
    cool.add_element('C', 1.0)

    for i_ele in sim_block.core_intake.values():
        if type(i_ele) == ghastly.core.CylCore:
            icyl_r = 100*i_ele.r
            icyl_h = 100*i_ele.h
        if type(i_ele) == ghastly.core.ConeCore:
            icone_h = 100*i_ele.h
            icone_ru = 100*i_ele.r_upper
            icone_rl = 100*i_ele.r_lower
    for m_ele in sim_block.core_main.values():
        if type(m_ele) == ghastly.core.CylCore:
            mcyl_r = 100*m_ele.r
            mcyl_h = 100*m_ele.h
        if type(m_ele) == ghastly.core.ConeCore:
            mcone_h = 100*m_ele.h
            mcone_ru = 100*m_ele.r_upper
            mcone_rl = 100*m_ele.r_lower
            mcone_zl = 100*m_ele.z_min
    for o_ele in sim_block.core_outtake.values():
        if type(o_ele) == ghastly.core.CylCore:
            ocyl_r = 100*o_ele.r
            ocyl_h = 100*o_ele.h

    #intake
    icyl_side = openmc.ZCylinder(r = icyl_r)
    icyl_top = openmc.ZPlane(z0 = (mcone_zl + mcone_h + mcyl_h + icone_h + icyl_h))
    icyl_bot = openmc.ZPlane(z0 = (mcone_zl + mcone_h + mcyl_h + icone_h))
    icone_H = -(icone_rl*icone_h)/(icone_ru-icone_rl)
    icone_z0 = mcone_zl + mcone_h + mcyl_h + icone_H
    icone_r2 = (icone_rl/icone_H)**2
    icone_side = openmc.ZCone(z0 = icone_z0, r2 = icone_r2)

    #main
    mcyl_side = openmc.ZCylinder(r = mcyl_r)
    mcyl_top = openmc.ZPlane(z0 = (mcone_zl + mcone_h + mcyl_h))
    mcyl_bot = openmc.ZPlane(z0 = (mcone_zl + mcone_h))
    mcone_H = -(mcone_ru*mcone_h)/(mcone_rl-mcone_ru)
    mcone_z0 = mcone_h - mcone_H
    mcone_r2 = (mcone_ru/mcone_H)**2
    mcone_side = openmc.ZCone(z0 = mcone_z0, r2 = mcone_r2)
    mcone_bot = openmc.ZPlane(z0 = mcone_zl)

    #outtake
    ocyl_side = openmc.ZCylinder(r = ocyl_r)
    ocyl_bot = openmc.ZPlane(z0 = (mcone_zl - ocyl_h))


    icyl_reg = -icyl_side & -icyl_top & +icyl_bot
    icone_reg = -icone_side & -icyl_bot & +mcyl_top
    mcyl_reg = -mcyl_side & -mcyl_top & +mcyl_bot
    mcone_reg = -mcone_side & -mcyl_bot & +mcone_bot
    ocyl_reg = -ocyl_side & -mcone_bot & +ocyl_bot
    active_core_reg = icyl_reg | icone_reg | mcyl_reg | mcone_reg | ocyl_reg
    active_core = openmc.Cell(region=active_core_reg)

    refl_side = openmc.ZCylinder(r = (1.5*mcyl_r))
    refl_top = openmc.ZPlane(z0 = (mcone_zl + mcone_h + mcyl_h + 
                               icone_h + icyl_h + 0.25*mcyl_h))
    refl_bot = openmc.ZPlane(z0 = (mcone_zl - ocyl_h - 0.25*mcyl_h))

    refl_reg = -refl_side & -refl_top & + refl_bot & ~active_core_reg
    refl = openmc.Cell(region = refl_reg, fill = mod)
    
    vel_all_time = velocity.vel_pebbles(directory, 
                                        velpath, 
                                        coordpath, 
                                        recircpath, 
                                        inputfile,
                                        sort_int=sort_int, 
                                        n_skip=n_skip,
                                        delimiter=delimiter, 
                                        skiprows=skiprows,
                                        n_recirc=n_recirc, 
                                        n_dump=n_dump, 
                                        dt=dt,
                                        recirc_hz=recirc_hz)
    n_files = len(vel_all_time.keys())
    vmax = max([max(list(vel.values())) for vel in vel_all_time.values()])
    vmin = min([min(list(vel.values())) for vel in vel_all_time.values()])

    cpath = path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(path.join(cpath, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])

    for tstep in range(n_files):
        if tstep%n_skip==0:
            data = read_input.read_lammps_bin(c_fnames[tstep], cpath)
            pebbles = {}
            for d in data:
                if int(d[0]) in vel_all_time[tstep] and abs(d[-2]) <= peb_r:
                    uid = int(d[0])
                    coord = 100*d[-3:]
                    v_z = vel_all_time[tstep][uid]
                    zone = max(d[2], 1)
                    layer = d[3]
                    pass_n = d[4]
                    recirc = d[5]
                    l_type = d[1]
                    pebbles[uid] = pebble.Pebble(uid=uid,
                                                 coords=coord,
                                                 velocity=v_z,
                                                 zone=zone,
                                                 layer=layer,
                                                 pass_num=pass_n,
                                                 recirc=recirc,
                                                 l_type=l_type)
            
            peb_mats, peb_list, peb_colors = create_vel_pebs(pebbles, 
                                                             peb_r,
                                                             colormap,
                                                             vmin,
                                                             vmax,
                                                             tstep)

            bbox_reg = icyl_reg | mcyl_reg | ocyl_reg
            bbox = openmc.Cell(region = bbox_reg)
            lower_left_core, upper_right_core = bbox.region.bounding_box
            shape_core = (3, 3, 3)
            pitch_core = (upper_right_core - lower_left_core)/shape_core
            core_lattice = openmc.model.create_triso_lattice(peb_list, 
                                                             lower_left_core, 
                                                             pitch_core, 
                                                             shape_core, 
                                                             cool)

            active_core.fill=core_lattice


            universe = openmc.Universe(cells=[active_core, refl])

            geometry = openmc.Geometry(universe)
            geometry.export_to_xml()

            materials = peb_mats + [cool, mod]
            openmc.Materials(materials).export_to_xml()

            settings = openmc.Settings()
            settings.run_mode = 'plot'
            settings.export_to_xml()

            xz_mats=openmc.Plot()
            xz_mats.basis='xz'
            xz_mats.origin=(0, 0, 50 )
            xz_mats.width = (260, 400)
            xz_mats.pixels = (520, 800)
            xz_mats.color_by='material'
            xz_mats.colors = {mod: 'darkgrey', 
                              cool: 'gainsboro'} | peb_colors
            pad = 15-len(str(tstep+start))
            fname = pad*'0'+str(tstep+start)+'_vel_xs.png'
            xz_mats.filename = fname
            plots = openmc.Plots([xz_mats])
            plots.export_to_xml()
            openmc.run()
    
    plt.scatter ([-1, 1], [-1, 1], c = [vmin, vmax], cmap = 'managua')
    plt.savefig('vel_xs_colorbar.png')
    plot.close()
    print(vmin, vmax)



def create_vel_pebs(pebbles, peb_r, colormap, vmin, vmax, tstep):
    '''
    create vel plot pebbles mats and cells
    '''
    
    id_counter = 10 + 1000*tstep
    abc = string.ascii_letters
    peb_out = openmc.Sphere(r=peb_r)
    peb_reg = -peb_out
    v_range = vmax - vmin
    cmap = mpl.colormaps[colormap].resampled(len(pebbles)*1000) 
    
    peb_mats = []
    peb_colors = {}
    peb_univs = {}
    for uid, peb in pebbles.items():
        alpha = (abc[rng.integers(low=0, high=len(abc))]+
                 abc[rng.integers(low=0, high=len(abc))]+
                 abc[rng.integers(low=0, high=len(abc))])
        numeric = (str(rng.integers(low=100, high=1000))+
                   str(rng.integers(low=100, high=1000))+
                   str(rng.integers(low=100, high=1000)))
        name = alpha+numeric
        peb_mat = openmc.Material(name=name, material_id=id_counter)
        peb_mat.set_density('g/cm3', 1.0)
        peb_mat.add_element('C', 1.00)
        peb_mats.append(peb_mat)

        f = 1 - (vmax - peb.velocity)/v_range
        c1 = int(255*cmap(f)[0])
        c2 = int(255*cmap(f)[1])
        c3 = int(255*cmap(f)[2])
        color = tuple([c1, c2, c3])
        peb_colors[peb_mat] = color


        id_counter += 10

        peb_cell = openmc.Cell(region = peb_reg, 
                               fill = peb_mat)
        peb_univ = openmc.Universe(cells=[peb_cell])
        peb_univs[uid] = peb_univ
    
    peb_list = []
    for uid in pebbles.keys():
        peb_list += [openmc.model.TRISO(peb_r, 
                                        peb_univs[uid], 
                                        pebbles[uid].coords)]
    

    return peb_mats, peb_list, peb_colors


def create_coord_csv(directory, coordpath, inputfile,
                     sort_int=-19, n_skip=1, delimiter=' ', skiprows=9,
                     n_recirc=2500000, n_dump=372366, dt=2.6855e-07,
                     recirc_hz=0.014, colormap='magma'):
    '''
    '''

    inp_path = path.expanduser(inputfile)
    inp_block = read_input.InputBlock(inp_path)
    sim_block = inp_block.create_obj()
    peb_r = 100*sim_block.r_pebble

    cpath = path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(path.join(cpath, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])
    #n_files = len(c_fnames)
    n_files = 3

    for tstep in range(n_files):
        if tstep%n_skip==0:
            data = read_input.read_lammps_bin(c_fnames[tstep], cpath)
            #pebbles = [['uid', 'coord', 'zone', 'layer', 'pass_n', 'recirc']]
            pebbles = []
            for d in data:
                if abs(100*d[-2]) <= peb_r:
                    uid = int(d[0])
                    x = 100*d[-3]
                    y = 100*d[-2]
                    z = 100*d[-1]
                    #zone = int(max(d[2], 1))
                    #layer = int(d[3])
                    #pass_n = int(d[4])
                    #recirc = int(d[5])
                    #l_type = int(d[1])
                    pebbles.append([peb_r, x, y, z])
            pad = 15 - len(str(tstep))
            fname = pad*'0'+str(tstep)+'_bed.csv'
            np.savetxt(fname, pebbles, delimiter = ',')

def para_bed_plot(csvpath):
    '''
    '''

    cpath = path.expanduser(csv_path)

    unsorted_csvnames = glob.glob(path.join(cpath, '*.csv'))
    csv_fnames = sorted(unsorted_csvnames, key=lambda x:x[sort_int:])

    #csvreader = OpenDataFile(csv_fnames)

    #points = tableToPoints()

    #pebbles = Sphere(Input=points)

    










            











    


