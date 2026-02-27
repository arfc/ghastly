import numpy as np
from ghastly import read_input
from ghastly import pebble
import os
import glob
import string
import h5py

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
        print(tstep)
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

    _write_data_xdmf(data_file, xmf_file)

def _write_data_xdmf(data_file, xmf_file):
    """
    Write an XDMF file to accompany the HDF5 data file.
    """

    with h5py.File(data_file, mode='r') as h5f:
        h5keys = list(h5f.keys())
        shape_n = len(h5f[h5keys[0]]['xyz'])
        n_dset = len(h5f)


    with open(xmf_file, mode='w') as xmf:
        xmf.write('<?xml version="1.0" ?>\n')
        xmf.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        xmf.write('<Xdmf Version="3.0">\n')
        xmf.write('  <Domain>\n')
        # CollectionType="Temporal" for paraview
        xmf.write('    <Grid Name="Pebbles" '
                  'GridType="Collection" CollectionType="Temporal">\n')
        for i, h5k in enumerate(h5keys):
            xmf.write('      <Grid Name="Pebbles" GridType="Uniform">\n')
            xmf.write(f'        <Time Value="{h5k}" />\n')

            xmf.write(f'        <Topology TopologyType="Polyvertex" '
                      'NumberOfElements="{shape_n}" />\n')
            xmf.write('        <Geometry GeometryType="XYZ">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n} 3">\n')
            xmf.write(f'            {data_file}:/{h5k}/xyz\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Geometry>\n')

            xmf.write('        <Attribute Name="Vx" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vx\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')
            xmf.write('        <Attribute Name="Vy" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vy\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')
            xmf.write('        <Attribute Name="Vz" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vz\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Vx Adjusted" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vx_adj\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')
            xmf.write('        <Attribute Name="Vy Adjusted" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vy_adj\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')
            xmf.write('        <Attribute Name="Vz Adjusted" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vz_adj\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="V Magnitude" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vmag\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')
            xmf.write('        <Attribute Name="V Magnitude Adjusted" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/vmag_adj\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="UID" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/uid\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Zone" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/zone\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Layer" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/layer\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Pass Number" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/pass_n\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Number of Recirculations" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="{shape_n}">\n')
            xmf.write(f'            {data_file}:/{h5k}/recirc_n\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Sim Step" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="1">\n')
            xmf.write(f'            {data_file}:/{h5k}/sim_step\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Reactor Time" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="1">\n')
            xmf.write(f'            {data_file}:/{h5k}/reactor_time\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('        <Attribute Name="Time Scaling Factor" AttributeType="Scalar" Center="Node">\n')
            xmf.write(f'          <DataItem Format="HDF" Dimensions="1">\n')
            xmf.write(f'            {data_file}:/{h5k}/time_scale\n')
            xmf.write('          </DataItem>\n')
            xmf.write('        </Attribute>\n')

            xmf.write('      </Grid>\n')
        
        xmf.write('        <Attribute Name="Pebble Radius" AttributeType="Scalar" Center="Node">\n')
        xmf.write(f'          <DataItem Format="HDF" Dimensions="1">\n')
        xmf.write(f'            {data_file}:/peb_r\n')
        xmf.write('          </DataItem>\n')
        xmf.write('        </Attribute>\n')

        xmf.write('        <Attribute Name="Sim dt" AttributeType="Scalar" Center="Node">\n')
        xmf.write(f'          <DataItem Format="HDF" Dimensions="1">\n')
        xmf.write(f'            {data_file}:/dt\n')
        xmf.write('          </DataItem>\n')
        xmf.write('        </Attribute>\n')

        xmf.write('        <Attribute Name="Recirculation Frequency" AttributeType="Scalar" Center="Node">\n')
        xmf.write(f'          <DataItem Format="HDF" Dimensions="1">\n')
        xmf.write(f'            {data_file}:/recirc_hz\n')
        xmf.write('          </DataItem>\n')
        xmf.write('        </Attribute>\n')
        
        xmf.write('    </Grid>\n')
        xmf.write('  </Domain>\n')
        xmf.write('</Xdmf>\n')



