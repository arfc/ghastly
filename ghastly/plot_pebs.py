import numpy as np
import openmc
from ghastly import read_input
from ghastly import pebble
from jinja2 import Environment, PackageLoader
from os import path
import subprocess
import glob

def plot_bed_xs(filepath, inputfile, color_crit, start = None, end = None):
    '''
    only prints a frame for every 10 data points - with 10 points per second,
    that's 1 frame per second
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


            











    


