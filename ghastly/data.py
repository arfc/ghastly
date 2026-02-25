import numpy as np
from ghastly import read_input
from ghastly import pebble
from os import path
import sys
import subprocess
import itertools
import glob
import string
import h5py
import XDMFWrite_h5py as xh
import pathlib

rng = np.random.default_rng()

def _recirc_data(recircpath, sort_int):
    '''
    pull lammps recirc data given recirc directory, return P = N pebs recirc'd
    each recirc
    '''
    rpath = path.expanduser(recircpath)
    unsorted_r_fnames = glob.glob(path.join(rpath, "*.bin"))
    r_fnames = sorted(unsorted_r_fnames, key=lambda x:x[sort_int:])
    P = [float(read_input.read_lammps_bin(rfile, 
                                    rpath, 
                                    skiprows=3, 
                                    max_rows=1)) for rfile in r_fnames]
    P_uids = [[int(uid) for uid in group] 
              for group in [read_input.read_lammps_bin(rfile, 
                                                       rpath, 
                                                       skiprows=9) 
                            for rfile in r_fnames]]
    return P, P_uids


def create_data_hdf5(inputfile, coordpath, recircpath, 
                     data_name="pebble_data", sort_int=-19, n_skip=1, 
                     delimiter=' ', skiprows=9, n_recirc=2500000, 
                     n_dump=250000, dt=2.4790e-07, recirc_hz=0.014, 
                     pattern=[0, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8], units = 'm'):
    '''
    create peb coord and reactor-time converted velocity hdf5 w/ extra data 
    needed for paraview plotting. values in cm for later openmc use. 
    also time step value csv
    '''

    inp_path = path.expanduser(inputfile)
    inp_block = read_input.InputBlock(inp_path)
    sim_block = inp_block.create_obj()
    P, _ = _recirc_data(recircpath, sort_int)
    

    cpath = path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(path.join(cpath, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])
    n_files = len(c_fnames)
    match units:
        case 'm':
            c = 100
        case 'cm':
            c = 1
        case 'in':
            c = 2.54
        case 'ft':
            c = 30.48
    time_key = []
    data_file = data_name + ".h5"
    xdmf_file = data_name + "xdmf"
    for tstep in range(n_files):
        i_r = int((tstep*(n_dump))//n_recirc)
        t_scale = P[i_r]/(n_recirc*dt*recirc_hz)
        sim_time = tstep*n_dump*dt
        reactor_time = sim_time*t_scale
        time_key.append([tstep, reactor_time])
        data = read_input.read_lammps_bin(c_fnames[tstep], cpath)
        pebbles = []
        for d in data:
            uid = int(d[pattern[0]])
            x = c*d[pattern[1]]
            y = c*d[pattern[2]]
            z = c*d[pattern[3]]
            zone = int(d[pattern[4]])
            layer = int(d[pattern[5]])
            pass_n = int(d[pattern[6]])
            recirc = int(d[pattern[7]])
            vx = c*d[pattern[8]]
            vx_adj = vx/t_scale
            vy = c*d[pattern[9]]
            vy_adj = vy/t_scale
            vz = c*d[pattern[10]]
            vz_adj = vz/t_scale
            pebbles.append([uid, x, y, z, vx, vy, vz, vx_adj, vy_adj, vz_adj,
                            c*sim_block.r_pebble, zone, layer, pass_n, recirc])


        pad = 20 - len(str(tstep))
        set_name = "/"+str(tstep)
        with h5py.File(data_file, mode='a') as datah5:
            try:
                dset = datah5.create_dataset(set_name, data=pebbles)
                dset.attrs['sim_step'] = tstep
                dset.attrs['reactor_time'] = reactor_time
                dset.attrs['time_scale'] = t_scale
                dset.attrs['dt'] = dt
                dset.attrs['recirc_hz'] = recirc_hz
            except:
                pass

            #with xh.Grid(pathlib.Path(xdmf_file)) as xhf:
            #    xhf += xh.Unstructured(datah5[set_name], element_type="Quadrilateral")
            #    xhf += xh.Attribute(datah5['sim_step'], "Cell")
            #    xhf += xh.Attribute(datah5['reactor_time'], "Cell")
            #    xhf += xh.Attribute(datah5['time_scale'], "Cell")
            #    xhf += xh.Attribute(datah5['dt'], "Cell")
            #    xhf += xh.Attribute(datah5['recirc_hz'], "Cell")


def interp_coordsremakelater(directory, coordpath, recircpath,
                     sort_int=-19, delimiter=' ', skiprows=9,
                     n_recirc=2500000, n_dump=372366, 
                     dt=2.6855e-07, recirc_hz=0.014):
    '''
    '''
    P, P_uids = _recirc_data(recircpath, sort_int)
    
    recirc_rtime = [(n_recirc*i*dt)*(p/(n_recirc*dt*recirc_hz)) 
                     for i, p in enumerate(P)]

    cpath = path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(path.join(cpath, "*.bin"))
    c_fnames = sorted(unsorted_c_fnames, key=lambda x:x[sort_int:])
    n_files = len(c_fnames)

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

