import numpy as np
import openmc
from ghastly import read_input
from ghastly import pebble
from jinja2 import Environment, PackageLoader
from os import path
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
from cycler import cycler

default_cmap = mpl.colormaps['magma']
default_cycler = (cycler(color=default_cmap(np.linspace(0,1,6))))

plt.rc('axes', prop_cycle=default_cycler)

def vel_profiler(directory, velpath, recircpath,
                 inputfile, vel_widths, sort_int=-19,
                 n_skip=1, delimiter = ' ', skiprows = 9,
                 n_dump = 372366, n_recirc = 2500000,
                 dt=2.6855e-07, recirc_hz=0.014):
    '''
    generate a velocity profile plot given a path to a top level velocity
    output directory, a list of the velocity subdirectories within it 
    (no slashes), the associated ghastly inputfile, and the widths of the
    v_regs, as a function of pebble diameter
    (so a region of width 2*d_peb is 2)
    '''
    
    sim_block, vel_fnames, vel_dir = _vel_setup(directory,
                                                velpath,
                                                inputfile,
                                                sort_int)

    peb_d = 2*100*sim_block.r_pebble
    r_regs = len(velpath)
    n_files = len(vel_fnames[velpath[0]])

    rpath = path.expanduser(recircpath)
    unsorted_r_fnames = glob.glob(path.join(rpath, "*.bin"))
    r_fnames = sorted(unsorted_r_fnames, key=lambda x:x[sort_int:])
    P = [float(read_input.read_lammps_bin(rfile, 
                                    rpath, 
                                    skiprows=3, 
                                    max_rows=1)) for rfile in r_fnames]
    
    vel_all_time = {}
    for i in range(n_files):
        i_r = int((i*(n_dump))//n_recirc)
        t_scale = P[i_r]/(n_recirc*dt*recirc_hz)
        if i%n_skip==0:
            vel_map = {}
            for j, vpath in enumerate(velpath):
                bin_dir = path.join(vel_dir, vpath)
                bin_fname = vel_fnames[vpath][i]
                data = read_input.read_lammps_bin(bin_fname,bin_dir)
            
                raw_data = {}
                for d in data:
                    v_z = 100*d[-1]/t_scale
                    z = 100*d[1]
                    q = int(z//peb_d)
                    if q in vel_map:
                        pass
                    else:
                        vel_map[q] = r_regs*[float("nan")]

                    if q in raw_data:
                        raw_data[q].append(v_z)
                    else:
                        raw_data[q] = [v_z]

                for bin_k, bin_v in raw_data.items():
                    n_pebs = len(bin_v)
                    tot_vel = sum(bin_v)
                    avg = tot_vel/n_pebs
                    vel_map[bin_k][j] = avg

        
            ticks = []
            for j, _ in enumerate(vel_widths):
                ticks.append(sum(vel_widths[:j]))
            r = peb_d*np.array(ticks+[ticks[-1]+vel_widths[-1]])
        
            sorted_vel = dict(sorted(vel_map.items(), reverse=True))

            y = peb_d*(np.array(list(sorted_vel.keys())+[-23.0]))
            z_vel = list(sorted_vel.values())
        
            plt.pcolormesh(r, y, z_vel, vmin=-100, vmax=100, shading='flat', 
                           norm='symlog')
            plt.set_cmap('managua')
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

def vel_line_plotter(directory, velpath, coordpath, recircpath, inputfile, 
                     sort_int=-19, n_skip=1,
                     delimiter=' ', skiprows=9, 
                     z_bin_coeff=5, r_bin_coeff=1,
                     n_recirc=2500000, n_dump=372366, 
                     dt=2.6855e-07, recirc_hz=0.014):
    '''
    plot a courser set of velocities as line plots (possibly also 2d on 3d
    bar graphs??)
    '''

    sim_block, vel_fnames, vel_dir = _vel_setup(directory,
                                                velpath,
                                                inputfile,
                                                sort_int)

    peb_d = 2*100*sim_block.r_pebble
    n_files = min([len(vel_fnames[velpath[n]]) for n in range(len(velpath))])
    #n_files = 100
    r_width = r_bin_coeff*peb_d
    z_width = z_bin_coeff*peb_d
    
    rpath = path.expanduser(recircpath)
    unsorted_r_fnames = glob.glob(path.join(rpath, "*.bin"))
    r_fnames = sorted(unsorted_r_fnames, key=lambda x:x[sort_int:])
    P = [float(read_input.read_lammps_bin(rfile, 
                                    rpath, 
                                    skiprows=3, 
                                    max_rows=1)) for rfile in r_fnames]
    
    vel_all_time = {}
    for i in range(n_files):
        i_r = int((i*(n_dump))//n_recirc)
        t_scale = P[i_r]/(n_recirc*dt*recirc_hz)
        if i%n_skip==0:
            vel_map = {}
            sort_data = []
            raw_data = {}
            for vpath in velpath:
                bin_dir = path.join(vel_dir, vpath)
                bin_fname = vel_fnames[vpath][i]
                data = read_input.read_lammps_bin(bin_fname,bin_dir)
            
                for d in data:
                    uid = int(d[0])
                    z = 100*d[1]
                    v_z = 100*d[-1]/t_scale
                    q = z_width*(z//z_width + 0.5)
                    if q in raw_data:
                        raw_data[q][uid] = v_z
                    else:
                        raw_data[q] = {uid:v_z}

                sort_data.append(dict(sorted(raw_data.items(), reverse=True)))

            for combo_dict in tuple(sort_data):
                for z_bin, v_dict in combo_dict.items():
                    vel_map[z_bin] = {}
                    for uid, v_z in v_dict.items():
                        vel_map[z_bin][uid] = v_z
                
            vel_all_time[i] = vel_map

    cpath = path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(path.join(cpath, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])

    recirc_bin = {}
    for tstep, vmap in vel_all_time.items():
        recirc_bin[tstep] = {}
        data = read_input.read_lammps_bin(c_fnames[tstep], cpath)
        for zi, v_dict in vmap.items():
            recirc_bin[tstep][zi] = {}
            temp_r = {}
            for d in data:
                uid = int(d[0])
                if uid in v_dict:
                    r = ((100*d[-3])**2 + (100*d[-2])**2)**0.5
                    r_bin = r_width*(r//r_width + 0.5)

                    if r_bin in temp_r:
                        temp_r[r_bin].append(v_dict[uid])
                    else:
                        temp_r[r_bin] = [v_dict[uid]]

            temp_r = dict(sorted(temp_r.items(), reverse=True)) 
            for ri, vels in temp_r.items():
                v_z_avg = sum(vels)/len(vels)
                recirc_bin[tstep][zi][ri] = v_z_avg
    
    #know from experience the pebbles only take about .3 seconds to hit
    #the bed, so if the whole thing is 0.67 seconds, splitting it into 0.1
    # second *bins* is too broad - more bins, slightly fewer points per bin
    #n_dump not being a factor of n_recirc means it staggers over
    #successive cycles
    n_05s = int(0.05/dt)
    n_guess = n_recirc//n_05s
    factor_guess = np.linspace(n_guess, 2*n_guess, n_guess+1)
    factors = np.delete(factor_guess, np.where(n_recirc%factor_guess != 0))
    check = [abs(f/n_guess-1) for f in factors]
    n_bins = factors[np.where(check == min(check))]
    bin_w = n_recirc/n_bins
    print(n_bins, bin_w)
    

    cycle_bins = {}
    for tstep, v_bin in recirc_bin.items():
        i_recirc = (tstep*n_dump)//n_recirc
        i_bin = int((tstep*n_dump - i_recirc*n_recirc)//bin_w)
        if i_bin in cycle_bins:
            cycle_bins[i_bin].append(v_bin)
        else:
            cycle_bins[i_bin] = [v_bin]
    
    recirc_cycle = {}
    for i_bin, c_bin in cycle_bins.items():
        recirc_cycle[i_bin] = {}
        for all_cycles in tuple(c_bin):
            for z_bin, r_bins in all_cycles.items():
                if z_bin in recirc_cycle[i_bin]:
                    pass
                else:
                    recirc_cycle[i_bin][z_bin] = {}
                for r_bin, vel in r_bins.items():
                    if r_bin in recirc_cycle[i_bin][z_bin]:
                        recirc_cycle[i_bin][z_bin][r_bin].append(vel)
                    else:
                        recirc_cycle[i_bin][z_bin][r_bin] = [vel]

    recirc_avg ={}
    for i_bin, z_bins in recirc_cycle.items():
        recirc_avg[i_bin] = {}
        for z_bin, r_bins in z_bins.items():
            recirc_avg[i_bin][z_bin] = {}
            for r_bin, vel in r_bins.items():
                avg_v = sum(vel)/len(vel)
                recirc_avg[i_bin][z_bin][r_bin] = avg_v 

    for i_bin, z_bins in recirc_avg.items():
        zs = list(z_bins.keys())
        r_fig = (max(list(z_bins[zs[0]].keys())))/peb_d
        z_fig = (max(zs) - min(zs))/peb_d
        fig, axs = plt.subplots(len(zs), 1, 
                                sharex = True, 
                                figsize=(0.2*r_fig, 0.2*z_fig),
                                tight_layout=True)
        for i, zi in enumerate(zs):
            rs = list(z_bins[zi].keys())
            v_zs = list(z_bins[zi].values())
            axs[i].plot(rs, v_zs)
            axs[i].set_ylim(-0.02, 0.002)
        
        fig.suptitle('Velocity Profile During Recirc Event')
        fig.supxlabel('R [cm]')
        fig.supylabel('Velocity [cm/s]')
        pad = 4 - len(str(i_bin))
        fname = pad*"0"+str(i_bin)+"cycle_vz.png"
        plt.savefig(fname, dpi = 450)
        plt.close()


def vel_pebbles(directory, velpath, coordpath, recircpath, inputfile,
                sort_int=-19, n_skip=1,
                delimiter=' ', skiprows=9,
                n_recirc=2500000, n_dump=372366, 
                dt=2.6855e-07, recirc_hz=0.014):


    '''
    get per-pebble velocity data, scaled to real time, to use with
    pebble bed xs plotting. (note to self, the pebble class has a
    velocity attribute you use there, so you need to pass the right
    velocity data as a function of uid
    '''

    sim_block, vel_fnames, vel_dir = _vel_setup(directory,
                                                velpath,
                                                inputfile,
                                                sort_int)

    peb_d = 2*100*sim_block.r_pebble
    #n_files = min([len(vel_fnames[velpath[n]]) for n in range(len(velpath))])
    n_files = 100
    r_width = r_bin_coeff*peb_d
    z_width = z_bin_coeff*peb_d
    
    rpath = path.expanduser(recircpath)
    unsorted_r_fnames = glob.glob(path.join(rpath, "*.bin"))
    r_fnames = sorted(unsorted_r_fnames, key=lambda x:x[sort_int:])
    P = [float(read_input.read_lammps_bin(rfile, 
                                    rpath, 
                                    skiprows=3, 
                                    max_rows=1)) for rfile in r_fnames]
    vel_all_time = {}
        for i in range(n_files):
            i_r = int((i*(n_dump))//n_recirc)
            t_scale = P[i_r]/(n_recirc*dt*recirc_hz)
            if i%n_skip==0:
                vel_map = {}
                sort_data = []
                raw_data = {}
                for vpath in velpath:
                    bin_dir = path.join(vel_dir, vpath)
                    bin_fname = vel_fnames[vpath][i]
                    data = read_input.read_lammps_bin(bin_fname,bin_dir)

                for d in data:
                    uid = int(d[0])
                    v_z = 100*d[-1]/t_scale
                    vel_map[uid] = v_z
                vel_all_time[i] = vel_map
    
    return vel_all_time

    

def __NOT_IMPLEMENTED_track_avalanches(directory, velpath, 
                                       coorddir, inputfile,
                                       sort_int=-19, n_skip = 1,
                                       delimiter = ' ', skiprows = 9):
    '''

    NOT IMPLEMENTED
    use velocity data to look at the KE of the pebs in the main core
    and try to find avalanches, assumes all pebbles identical, skips 1/2m
    factor in KE
    '''

    sim_block, vel_fnames, vel_dir = _vel_setup(directory,
                                                velpath,
                                                inputfile,
                                                sort_int)

    peb_d = 2*100*sim_block.r_pebble
    n_files = len(vel_fnames[velpath[0]])

    core_zmax = 100*max([core.z_max for core in sim_block.core_main.values()])
    core_zmin = 100*max([core.z_min for core in sim_block.core_main.values()])
    
    kinetic_energy = {}
    impulses = {}
    for i in range(n_files):
        if i%n_skip == 0:
            ke_map = {}
            for vpath in velpath:
                bin_dir = path.join(vel_dir, vpath)
                bin_fname = vel_fnames[vpath][i]
                data = read_input.read_lammps_bin(bin_fname,bin_dir)

                for d in data:
                    if 100*d[1] <= core_zmax and 100*d[1] >= core_zmin:
                        ke_map[d[0]] = sum((100*d[:2])**2)
            kinetic_energy[i] = ke_map

            avg_ke = sum(list(ke_map.values()))/len(ke_map.values())
            impulses[i] = {}
            for uid, ke in ke_map.items():
                if ke/avg_ke >=2.45:
                    impulses[i][uid] = {}
                    impulses[i][uid]['ke'] = ke

    for impulse in impulses.values():
        print(len(impulse))
    print()
    print()

    cdir =  path.expanduser(coorddir)
    unsorted_c_fnames = glob.glob(path.join(cdir, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])

    impulse_coords = {}
    for i, c_f in enumerate(c_fnames):
        data = read_input.read_lammps_bin(c_f, cdir)
        xs = []
        ys = []
        zs = []
        for d in data:
            if d[0] in impulses[i]:
                coord = 100*d[-3:]
                impulses[i][d[0]]['coord'] = coord
                xs.append(coord[0])
                ys.append(coord[1])
                zs.append(coord[2])

        #note to self for next time: try having a reference map of coords
        #for each time step to compare to.  Then, when you notice what might 
        #be an impulse based on KE, you can add the uid/coord to the pebble at
        #the current step, i, but also at the i+/-n ones.  maybe only put a
        #list of uiuds in for each timestep (check if uid already in there
        #first), then when it comes time to plot, pull coords out based on uid

        #maybe also try to do some sort of thing where you break ke up into
        #regions, and compare against local avg_ke?  could
        #also try adding the criterion that plotted pebbles need to be within 
        # a little more than a pebble diameter of another pebble?
        
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(xs, ys, zs, marker='o')
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('Y [cm]')
        ax.set_zlabel('Z [cm]')
        ax.set_zlim(54, 200)
        ax.set_xlim(-120,120)
        ax.set_ylim(-120,120)
        pad = 15 - len(str(i))
        filename = "impulse_"+pad*"0"+str(i)+".png"
        plt.savefig(filename)
        plt.close(fig)



    









    
        

        







