import numpy as np
import openmc
import ghastly
from ghastly import read_input
from ghastly import pebble
from ghastly import velocity
from jinja2 import Environment, PackageLoader
from os import path
import sys
from paraview.simple import *
import subprocess
import itertools
import glob
import string
import csv
import ast
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

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



def create_coord_csv(directory, coordpath, recircpath, inputfile,
                     sort_int=-19, n_skip=1, delimiter=' ', skiprows=9,
                     reactor_step = 24, n_recirc=2500000, n_dump=372366, 
                     dt=2.6855e-07, recirc_hz=0.014, colormap='magma'):
    '''
    '''

    inp_path = path.expanduser(inputfile)
    inp_block = read_input.InputBlock(inp_path)
    sim_block = inp_block.create_obj()
    peb_r = 100*sim_block.r_pebble

    rpath = path.expanduser(recircpath)
    unsorted_r_fnames = glob.glob(path.join(rpath, "*.bin"))
    r_fnames = sorted(unsorted_r_fnames, key=lambda x:x[sort_int:])
    P = [float(read_input.read_lammps_bin(rfile, 
                                    rpath, 
                                    skiprows=3, 
                                    max_rows=1)) for rfile in r_fnames]
    

    cpath = path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(path.join(cpath, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])
    n_files = len(c_fnames)
    #n_files = 10

    for tstep in range(n_files):
        i_r = int((tstep*(n_dump))//n_recirc)
        t_scale = P[i_r]/(n_recirc*dt*recirc_hz)
        sim_time = tstep*n_dump*dt
        reactor_time = sim_time*t_scale
        if tstep%n_skip==0:
            data = read_input.read_lammps_bin(c_fnames[tstep], cpath)
            pebbles = []
            for d in data:
                uid = int(d[0])
                x = 100*d[-3]
                y = 100*d[-2]
                z = 100*d[-1]
                zone = int(max(d[2], 1))
                layer = int(d[3])
                pass_n = int(d[4])
                recirc = int(d[5])
                pebbles.append([x, y, z, peb_r,
                          zone, layer, pass_n, recirc])


        pad = 15 - len(str(tstep))
        fname = pad*'0'+str(tstep)+'_bed.csv'
        with open(fname, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(pebbles)

def interp_coords(directory, coordpath, recircpath,
                     sort_int=-19, delimiter=' ', skiprows=9,
                     n_recirc=2500000, n_dump=372366, 
                     dt=2.6855e-07, recirc_hz=0.014):
    '''
    '''
    rpath = path.expanduser(recircpath)
    unsorted_r_fnames = glob.glob(path.join(rpath, "*.bin"))
    r_fnames = sorted(unsorted_r_fnames, key=lambda x:x[sort_int:])
    P = [float(read_input.read_lammps_bin(rfile, 
                                    rpath, 
                                    skiprows=3, 
                                    max_rows=1)) for rfile in r_fnames]
    
    recirc_rtime = [(n_recirc*i*dt)*(p/(n_recirc*dt*recirc_hz)) 
                     for i, p in enumerate(P)]
    recirc_uids = [[int(uid) for uid in group] 
                   for group in [read_input.read_lammps_bin(rfile, 
                                                            rpath, 
                                                            skiprows=9) 
                   for rfile in r_fnames]]

    cpath = path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(path.join(cpath, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])
    n_files = len(c_fnames)
    #n_files = 30

    adj_scale = P[0]/(n_recirc*dt*recirc_hz)
    adj_simtime = n_dump*dt
    adj_dt = (adj_simtime*adj_scale)//1
    for i, starting_uids in enumerate(recirc_uids[0:-1]):
        uids = starting_uids
        recirc = {}
        raw_points = {u:[] for u in uids}
        raw_rts = []
        range_start = int(i + (0.5*n_recirc)//n_dump)

        init_simt = range_start*n_dump*dt
        init_i_r = int((range_start*n_dump)//n_recirc)
        init_tscale = P[init_i_r]/(n_recirc*dt*recirc_hz)
        init_rt = init_simt*init_tscale
        for tstep in range(n_files)[range_start:]:
            if len(uids) == 0:
                break
            i_r = int((tstep*(n_dump))//n_recirc)
            t_scale = P[i_r]/(n_recirc*dt*recirc_hz)
            
            if tstep == range_start:
                reactor_time = init_rt
            else:
                reactor_time += (n_dump*dt*t_scale)
            if reactor_time < recirc_rtime[i+1]:
                pass
            elif reactor_time >= recirc_rtime[i+1]:
                raw_rts.append(reactor_time)
                data = read_input.read_lammps_bin(c_fnames[tstep], cpath)
                for d in data:
                    uid = int(d[0])
                    recirc_n = int(d[5])
                    if uid in uids:
                        if len(raw_points[uid]) == 0:
                            recirc[uid] = recirc_n
                            x = d[-3]
                            y = d[-2]
                            r = 100*(x**2+y**2)**0.5
                            z = 100*d[-1]
                            raw_points[uid] = [np.array([r, z])]
                        else:
                            if recirc_n != recirc[uid]:
                                uids.remove(uid)
                            elif recirc_n == recirc[uid]:
                                x = d[-3]
                                y = d[-2]
                                r = 100*(x**2+y**2)**0.5
                                z = 100*d[-1]
                                raw_points[uid].append(np.array([r, z]))
                                next_rt = raw_rts[-1] + (n_dump*dt*t_scale)
        points = {}
        for uid, raw_pts in raw_points.items():
            points[uid] = []
            for j, pt in enumerate(raw_pts):
                if j == 0:
                    points[uid] = [[float(pt[0]), float(pt[1])]]
                else:
                    adj_t = raw_rts[0] + j*adj_dt
                    for t_index, rt in enumerate(raw_rts):
                        if rt > adj_t:
                            rt0 = raw_rts[t_index-1]
                            rt1 = rt
                            break

                    adj_r = (points[uid][-1][0]*((rt1-adj_t)/(rt1-rt0)) + 
                              pt[0]*((adj_t-rt0)/(rt1-rt0)))
                    adj_z = (points[uid][-1][1]*((rt1-adj_t)/(rt1-rt0)) + 
                              pt[1]*((adj_t-rt0)/(rt1-rt0)))
                    points[uid].append([float(adj_r), float(adj_z)])
        pad = 7 - len(str(i))
        fname = pad*'0'+str(i)+'-transit.csv'
        fields = [uidkey for uidkey in points.keys()]
        with open(fname, "w", newline='') as f:
            writer = csv.writer(f)
            writer.writerow(points.keys())
            writer.writerows(zip(*points.values()))
                

    print(adj_dt)
    return adj_dt

def plot_velocity(transitpath, peb_d, height_lim, adj_dt, 
                  sort_int = -4 ):
    '''
    '''
    tpath = path.expanduser(transitpath)
    unsorted_t_fnames = glob.glob(path.join(tpath, "*transit.csv"))
    transitfiles = sorted(unsorted_t_fnames, key=lambda x:x[:sort_int])
    binsize = 4*peb_d
    r_bins = [(0, 4*peb_d), 
              (4*peb_d, 8*peb_d), 
              (8*peb_d, 12*peb_d),
              (12*peb_d, 15*peb_d),
              (15*peb_d, 18*peb_d),
              (18*peb_d, 20*peb_d)] # [lower, upper) r_bin bounds
    r_labels = [bound[0]+0.5*(bound[1]-bound[0]) for bound in r_bins]
    ax_layers = {}
    for transitfile in transitfiles:
        csvcols = pd.read_csv(transitfile, nrows=0).columns.tolist()
        converter = {col:ast.literal_eval for col in csvcols}
        df = pd.read_csv(transitfile, converters=converter)
        strdict = df.to_dict()
        transit = {}
        for uidstr in strdict.keys():
            transit[int(uidstr)] = list(strdict[uidstr].values())
        for uid, points in transit.items():
            for point in points:
                if point[1] >= height_lim:
                    continue
                ax_bin = point[1]//binsize
                r_criter = [(point[0]>=r_bin[0] and point[0]<r_bin[1])
                                  for r_bin in r_bins]
                r_binkey = r_labels[np.where(r_criter)[0][0]]
                if ax_bin in ax_layers:
                    if r_binkey in ax_layers[ax_bin]:
                        if uid in ax_layers[ax_bin][r_binkey]:
                            ax_layers[ax_bin][r_binkey][uid].append(point)
                        else:
                            ax_layers[ax_bin][r_binkey][uid] = [point]
                    else:
                        ax_layers[ax_bin][r_binkey] = {uid:[point]}
                else:
                    ax_layers[ax_bin] = {r_binkey:{uid:[point]}}
    vel_z = {}
    for ax_bin, bin_data in ax_layers.items():
        vel_z[ax_bin] = {}
        for r_binkey, pebbles in bin_data.items():
            v_z = []
            for uid, points in pebbles.items():
                if len(points) == 1:
                    continue
                for i, point in list(enumerate(points))[1:]:
                    z_disp = -abs(point[1] - points[i-1][1])
                    v_z.append(z_disp/adj_dt)
            if len(v_z) == 0:
                continue
            v_z_avg = sum(v_z)/len(v_z)
            vel_z[ax_bin][r_binkey] = v_z_avg


    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for ax_bin, bin_data in vel_z.items():
        if ax_bin < 0:
            continue
        r_s = np.array(list(bin_data.keys()))
        vz_s = np.array(list(bin_data.values()))
        ax.plot(r_s, vz_s,  
                zs=(ax_bin*binsize+(binsize/2)), zdir='y', c=(0,0,0))
    ax.set_zlabel('Z-Velocity [cm/s]', labelpad=15)
    ax.set_xlabel('Radius [cm]', labelpad=15)
    ax.set_ylabel('Z [cm]', labelpad=15)
    fig.savefig('vel_profile.png')
    plt.close()

def plot_velocityv2(transitpath, peb_d, height_lim, adj_dt, 
                  sort_int = -4 ):
    '''
    '''
    #this is a temp function to try another method for making this
    #plot so I can compare - the change is too in-depth to just add a flag to
    #the original version.  In the end, there should only be one velocity
    #function that handles 
    tpath = path.expanduser(transitpath)
    unsorted_t_fnames = glob.glob(path.join(tpath, "*transit.csv"))
    transitfiles = sorted(unsorted_t_fnames, key=lambda x:x[:sort_int])
    binsize = 5*peb_d
    r_bins = [(0, 4*peb_d), 
              (4*peb_d, 8*peb_d), 
              (8*peb_d, 12*peb_d),
              (12*peb_d, 15*peb_d),
              (15*peb_d, 18*peb_d),
              (18*peb_d, 20*peb_d)] # [lower, upper) r_bin bounds
    r_keys = [bound[0]+0.5*(bound[1]-bound[0]) for bound in r_bins]

    #note to self - because you track by uid, when the pebbles recirculate,
    #i think you are lumping those recircs together even if they aren't in the
    #same r  r binning would help
    vz_unavg = {}
    for transitfile in transitfiles:
        csvcols = pd.read_csv(transitfile, nrows=0).columns.tolist()
        converter = {col:ast.literal_eval for col in csvcols}
        df = pd.read_csv(transitfile, converters=converter)
        strdict = df.to_dict()
        transit = {}
        for uidstr in strdict.keys():
            transit[int(uidstr)] = list(strdict[uidstr].values())
        for uid, points in transit.items():
            bin_points = {}
            for point in points:
                if point[1] >= height_lim:
                    continue
                ax_bin = point[1]//binsize
                r_criter = [(point[0]>=bound[0] and point[0]<bound[1])
                                  for bound in r_bins]
                r_binkey = r_keys[np.where(r_criter)[0][0]]
                if ax_bin in bin_points:
                    if r_binkey in bin_points[ax_bin]:
                        bin_points[ax_bin][r_binkey].append(point)
                    else:
                        bin_points[ax_bin][r_binkey] = [point]
                else:
                    bin_points[ax_bin] = {r_binkey:[point]}

            for ax_bin, r_points in bin_points.items():
                if ax_bin not in vz_unavg:
                    vz_unavg[ax_bin] = {}
                for r_bin, pts in r_points.items():
                    zs = [pt[1] for pt in pts]
                    z_disp = min(zs)-max(zs)
                    peb_vz = z_disp/adj_dt#this is wrong, need to consider
                    #that it's not just one dt
                    if r_bin in vz_unavg:
                        vz_unavg[ax_bin][r_bin].append(peb_vz)
                    else:
                        vz_unavg[ax_bin][r_bin] = [peb_vz]
    print(vz_unavg)
    vel_z = {}
    for ax_bin, bin_data in vz_unavg.items():
        vel_z[ax_bin] = {'r':[], 'vz':[]}
        for r_bin, vzs in bin_data.items():
            if len(vzs) == 0:
                continue
            vel_z[ax_bin]['r'].append(r_bin)
            vz = sum(vzs)/len(vzs)
            vel_z[ax_bin]['vz'].append(vz)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for ax_bin, bin_data in vel_z.items():
        ax.plot(bin_data['r'], bin_data['vz'],  
                zs=(ax_bin*binsize+(binsize/2)), zdir='x', c=(0,0,0))
    ax.set_zlabel('velocity')
    ax.set_xlabel('ax_bin')
    ax.set_ylabel('radius')
    fig.savefig('vel_profilesv2_test.png')
    plt.close()


                

def plot_streamlines(transitpath, peb_d, height_lim, adj_dt, sort_int = -4):
    '''
    '''

    tpath = path.expanduser(transitpath)
    unsorted_t_fnames = glob.glob(path.join(tpath, "*transit.csv"))
    transitfiles = sorted(unsorted_t_fnames, key=lambda x:x[:sort_int])
    r_bins = [(0, 4*peb_d), 
              (4*peb_d, 8*peb_d), 
              (8*peb_d, 12*peb_d),
              (12*peb_d, 15*peb_d),
              (15*peb_d, 18*peb_d),
              (18*peb_d, 20*peb_d)] # [lower, upper) r_bin bounds
    r_keys = [bound[0]+0.5*(bound[1]-bound[0]) for bound in r_bins]

    r_streamlines = {r_key:[] for r_key in r_keys}
    for transitfile in transitfiles:
        csvcols = pd.read_csv(transitfile, nrows=0).columns.tolist()
        converter = {col:ast.literal_eval for col in csvcols}
        df = pd.read_csv(transitfile, converters=converter)
        strdict = df.to_dict()
        streamlines = {}
        for uidstr in strdict.keys():
            streamlines[int(uidstr)] = list(strdict[uidstr].values())
        for uid, points in streamlines.items():
            for point in points:
                if point[1] < height_lim:
                    r_criter = [(point[0]>=r_bin[0] and point[0]<r_bin[1])
                                  for r_bin in r_bins]
                    r_binkey = r_keys[np.where(r_criter)[0][0]]
                    r_streamlines[r_binkey].append(points)
                    break
    
    for r_bin, all_lines in r_streamlines.items():
        for streamline in all_lines:
            r = list(zip(*streamline))[0]
            z = list(zip(*streamline))[1]
            plt.plot(r, z, c=(0,0,0))
        fname = 'streamlines_rbin_'+str(r_bin)+'.png'
        plt.savefig(fname)
        plt.close()


def para_bed_plot(csvpath, sort_int = -19):
    '''
    '''

    cpath = path.expanduser(csvpath)

    unsorted_csvnames = glob.glob(path.join(cpath, '*.csv'))
    csv_fnames = sorted(unsorted_csvnames, key=lambda x:x[sort_int:])

    for file in csv_fnames:
        csvreader = CSVReader(FileName=file)
        csvreader.HaveHeaders=0
        points = TableToPoints(Input=csvreader)
        points.XColumn = 'Field 1' 
        points.YColumn = 'Field 2' 
        points.ZColumn = 'Field 3'
        points.UpdatePipeline()
        pebbles = Sphere()
        pebbles.Radius=3.0
        pebbles.PhiResolution = 10
        pebbles.ThetaResolution = 10

    










            











    


