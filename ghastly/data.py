import numpy as np
from ghastly import read_input
from ghastly import pebble
from ghastly.main import _templater
from jinja2 import Environment, PackageLoader
import os
import glob
import string
import h5py

env = Environment(loader=PackageLoader('ghastly', 'templates'))

rng = np.random.default_rng()

def _recirc_data(recircpath, sort_int):
    '''
    pull lammps recirc data given recirc directory, return P = N pebs recirc'd
    each recirc
    '''
    rpath = os.path.expanduser(recircpath)
    unsorted_r_fnames = glob.glob(os.path.join(rpath, "*.bin"))
    r_fnames = sorted(unsorted_r_fnames, key=lambda x:x[sort_int:])
    P = [float(read_input.read_lammps_bin(rfile, rpath, skiprows=3, 
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

    inp_path = os.path.expanduser(inputfile)
    inp_block = read_input.InputBlock(inp_path)
    sim_block = inp_block.create_obj()
    P, _ = _recirc_data(recircpath, sort_int)
    

    cpath = os.path.expanduser(coordpath)
    unsorted_c_fnames = glob.glob(os.path.join(cpath, "*.bin"))
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

    data_file = data_name + ".h5"
    xmf_file = data_name + ".xmf"
    try:
        os.remove(data_file)
    except OSError:
        pass
    try:
        os.remove(xmf_file)
    except OSError:
        pass

    
    for tstep in range(n_files):
        i_r = int((tstep*(n_dump))//n_recirc)
        t_scale = P[i_r]/(n_recirc*dt*recirc_hz)
        sim_time = tstep*n_dump*dt
        reactor_time = sim_time*t_scale
        data = read_input.read_lammps_bin(c_fnames[tstep], cpath)
        uid = []
        coord = []
        vx, vy, vz, vmag = [], [], [], []
        zone = []
        layer = []
        pass_n = []
        recirc_n = []

        for d in data:
            uid.append(int(d[pattern[0]]))
            coord.append([c*d[pattern[1]], c*d[pattern[2]], c*d[pattern[3]]])
            zone.append(int(d[pattern[4]]))
            layer.append(int(d[pattern[5]]))
            pass_n.append(int(d[pattern[6]]))
            recirc_n.append(int(d[pattern[7]]))
            vxi, vyi, vzi = c*d[pattern[8]], c*d[pattern[9]], c*d[pattern[10]]
            vx.append(vxi)
            vy.append(vyi)
            vz.append(vzi)
            vmag.append(sum([vxi**2 + vyi**2 + vzi**2])**0.5)
            


        pad = 10 - len(str(tstep))
        group_name = pad*'0'+str(tstep)
        with h5py.File(data_file, mode='a') as h5f:
            if 'peb_r' not in h5f.attrs:
                h5f.attrs['peb_r'] = c*sim_block.r_pebble
            if 'dt' not in h5f.attrs:
                h5f.attrs['dt'] = dt
            if 'recirc_hz' not in h5f.attrs:
                h5f.attrs['recirc_hz'] = recirc_hz
            group = h5f.require_group(group_name)

            group.create_dataset('xyz', data=coord)

            group.create_dataset('vx', data=vx)
            group.create_dataset('vx_adj', 
                                 data=np.array(vx)/t_scale)

            group.create_dataset('vy', data=vy)
            group.create_dataset('vy_adj', 
                                 data=np.array(vy)/t_scale)

            group.create_dataset('vz', data=vz)
            group.create_dataset('vz_adj', 
                                 data=np.array(vz)/t_scale)

            group.create_dataset('vmag', data=vmag)
            group.create_dataset('vmag_adj', 
                                 data=np.array(vmag)/t_scale)

            group.create_dataset('uid', data=uid, dtype=np.int_)
            group.create_dataset('zone',
                                 data=zone, dtype=np.int_)
            group.create_dataset('layer', 
                                 data=layer, dtype=np.int_)
            group.create_dataset('pass_n', 
                                 data=pass_n, dtype=np.int_)
            group.create_dataset('recirc_n', 
                                 data=recirc_n, dtype=np.int_)
            group.attrs['sim_step'] = tstep
            group.attrs['reactor_time'] = reactor_time
            group.attrs['time_scale'] = t_scale

    _write_data_xmf(data_file, xmf_file)

def _write_data_xmf(data_file, xmf_file):
    """
    Write an XMF file to accompany the HDF5 data file.

    Parameters
    ----------
    data_file : str
        Input .h5 file absolute directory
    xmf_file : str
        Name of output xmf file.

    Returns
    ----------
    """

    with h5py.File(data_file, mode='r') as h5f:
        h5keys = list(h5f.keys())
        shape_n = len(h5f[h5keys[0]]['xyz'])
        n_dset = len(h5f)
    
    params = {
        "data_file":data_file,
        "shape_n":shape_n,
        "h5keys":h5keys
    }

    _templater(params, xmf_file, "data/xmf_template.txt")
