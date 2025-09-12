import numpy as np
import openmc
from ghastly import read_input
from ghastly import pebble
from jinja2 import Environment, PackageLoader
from os import path
import glob
import matplotlib.pyplot as plt

def vel_profiler(directory, velpath, inputfile, vel_widths, sort_int=-19,
                 skip_factor=1, delimiter = ' ', skiprows = 9):
    '''
    generate a velocity profile plot given a path to a top level velocity
    output directory, a list of the velocity subdirectories within it 
    (no slashes), the associated ghastly inputfile, and the widths of the
    v_regs, as a function of pebble diameter
    (so a region of width 2*d_peb is 2)
    '''
    
    sim_block, vel_fnames, vel_dir = _setup_vel(directory,
                                                velpath,
                                                inputfile,
                                                sort_int)

    peb_d = 2*100*sim_block.r_pebble
    r_regs = len(velpath)
    n_files = len(vel_fnames[velpath[0]])
    
    vel_all_time = {}
    for i in range(n_files):
        if i%skip_factor==0:
            vel_map = {}
            for j, vpath in enumerate(velpath):
                bin_dir = path.join(vel_dir, vpath)
                bin_fname = vel_fnames[vpath][i]
                data = read_input.read_lammps_bin(bin_fname,bin_dir)
            
                raw_data = {}
                for d in data:
                    row = [int(d[0]),100*d[1], 100*d[2:]]
                    q = int(row[1]//peb_d)
                    if q in vel_map:
                        pass
                    else:
                        vel_map[q] = r_regs*[float("nan")]

                    if q in raw_data:
                        raw_data[q].append(row[2])
                    else:
                        raw_data[q] = [row[2]]

                for bin_k, bin_ in raw_data.items():
                    n_pebs = len(bin_v)
                    tot_vel = sum(bin_v)[2]
                    avg = tot_vel/n_pebs
                    vel_map[bin_k][j] = avg

        
            ticks = []
            for j, _ in enumerate(vel_widths):
                ticks.append(sum(vel_widths[:j]))
            r = peb_d*np.array(ticks+[ticks[-1]+vel_widths[-1]])
        
            sorted_vel = dict(sorted(vel_map.items(), reverse=True))

            y = peb_d*(np.array(list(sorted_vel.keys())+[-23.0]))
            z_vel = list(sorted_vel.values())
        
            plt.pcolormesh(r, y, z_vel, vmin=-100, vmax=50, shading='flat')
            plt.set_cmap('magma')
            plt.gca().set_aspect('equal')
            plt.xticks([0,40,80,120])
            plt.xlabel('Radius [cm]')
            plt.ylabel('Height[cm]')
            plt.ylim(-138.0,200.0)
            plt.title("Radial and Axial Velocity Profile")
            plt.colorbar()
            pad = 15 - len(str(i))
            fname = "vel_"+pad*"0"+str(i)+".png"
            plt.savefig(fname, bbox_inches='tight', dpi = 600)
            plt.close()

def _vel_setup(directory, velpath, inputfile, sort_int):
    '''
    initialize/read in what vel functions need to run (sans the vel_fnames)
    '''
    dpath = path.expanduser(directory)
    
    vel_fnames = {}
    for vpath in velpath:
        unsorted = glob.glob(path.join(dpath, vpath, "*.bin"))
        sorted_fnames = sorted(unsorted, key=lambda x:x[sort_int:])
        vel_fnames[vpath] = sorted_fnames


    inp_path = path.expanduser(inputfile)
    inp_block = read_input.InputBlock(inp_path)
    sim_block = inp_block.create_obj()

    return sim_block, vel_fnames, dpath

def track_avalanches(directory, velpath, inputfile, sort_int=-19):
    '''
    use velocity data to look at the KE of the pebs in the main core
    and try to find avalanches
    '''


    
        

        







