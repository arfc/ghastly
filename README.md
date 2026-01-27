# Ghastly

Ghastly is a python package that externally couples Discrete Element Method
(DEM) physics simulation tools to neutron transport and depletion tools
to model flowing pebbles in pebble-bed reactors (PBRs).

Ghastly is compatible with [LAMMPS](https://www.lammps.org/) and 
[OpenMC](https://openmc.org/). Ghastly automatically generates input files and
can read output files from compatible codes.  The `examples/` directory contains
Ghastly input files for a generic HTGR design.

Each template has some hard-coded settings, some of which are required to
create output files needed by Ghastly.  Users unfamiliar with compatible tools
should read the documentation carefully before making changes to templates.

Currently, Ghastly supports coupling LAMMPS to OpenMC only, but contributions via a Pull Request are welcome to expand capabilities.

Note that the paraview-plot.py file is meant to be loaded into your Paraview
GUI python script editor.

## Ghastly can:
- Fill both upwards and downwards flowing PBRs with an initial bed 
configuration
- Recirculate downwards flowing (HTGR) pebble beds.

## Ghastly cannot:
- Model PBRs with a central reflector
- Model PBRs with multiple inlet and outlet chutes

# Installation:

Ghastly is not packaged with any source code for compatible tools.  It is 
recommended, but optional, that users install Ghastly and the DEM and Monte Carlo
tools of choice inside a conda environment.  The installation tutorial below
assumes that you are using a conda environment.  If you are not, then you should
follow whatever the source documentation recommends for installing their individual
pieces of software.

Activate the conda environment before installing any dependencies.  The 
setup-example script is only meant to serve as an example, and was not written 
with the intention of running on an arbitrary system.

## Dependencies and Compatible Tool installation:

For using OpenMC and LAMMPS, install the following with conda before moving on using conda.
- python
- numpy
- scipy
- matplotlib
- cython
- pytest
- mpich, libmpich-dev, mpich-mpicc, mpich-mpicxx, mpich-mpifort
- zlib
- libpng (optional, see LAMMPS below)
- libjpeg-turbo (optional, see LAMMPS below)
- ffmpeg (optional, see LAMMPS below)
- vtk
- uncertainties
- lxml

### OpenMC:
The quick install guide for OpenMC is [here](https://docs.openmc.org/en/stable/quickinstall.html).
You will generally need to follow the instructions to build from source.  First,
you will need to install OpenMC's dependencies.
- Install hdf5 from source.  Download the hdf5 archive from the developer; extract it.
```
cd YOUR_HDF5_DIRECTORY/
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DZLIB_LIBRARY:FILEPATH=${CONDA_PREFIX}/lib/libz.so -DZLIB_INCLUDE_DIR:PATH=${CONDA_PREFIX}/include -DZLIB_USE_EXTERNAL=OFF -DHDF5_ENABLE_SZIP_SUPPORT=OFF ..
cd ..
cmake --build ./build --config Release -j 4
cmake --install ./build --config Release
conda deactivate
conda activate YOUR-ENV
h5dump --version

```

- Install mpi4py using `python -m pip install mpi4py`
- Install h5py using `cc=${CONDA_PREFIX} HDF5_MPI=ON HDF5_DIR=${CONDA_PREFIX} python -m pip install --no-binary=h5py h5py`.  This will install with support for MPI, and the `${CONDA_PREFIX}` options will
direct your system to the environment you've installed into so far.  Double 
check the installation worked with `python -c "import h5py; print(h5py.version.hdf5_version)"`

Now you are ready to install OpenMC. Follow the steps in the tutorial for 
installing from source, but, when first calling cmake in `build/`, set
`HDF5_ROOT=${CONDA_PREFIX}`, `-DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX}`, 
`-DOPENMC_USE_OPENMP=on`, and `-DOPENMC_USE_MPI=on`.  After `make install`,
return to the top level of the OpenMC source directory, and use `python -m pip install .`
to finish installing OpenMC.  Remember to use the `-e` flag if you want to install
in dev mode.

### LAMMPS:
The documentation page for building LAMMPS begins [here](https://docs.lammps.org/Build.html).
LAMMPS has many optional packages and build options which will not be discussed
in depth here beyond what should be installed if working with Ghastly.  Users 
should set `-D BUILD_SHARED_LIBS=yes`, `-D BUILD_MPI=yes`, `-D BUILD_OMP=yes`.  

The setup script example also turns the JPEG, PNG, and FFMPEG options on.  
These are not needed for Ghastly to work, though users may prefer to have the option.
If using JPEG, PNG, and FFMPEG, those should be installed before installing
LAMMPS.  Even if you are installing without FFMPEG support, you may still find
it a useful tool to have later when processing images.

As always, set`-D CMAKE_INSTALL_PREFIX` and `-D CMAKE_PREFIX_PATH` to
`${CONDA_PREFIX}`.  To use LAMMPS to model pebble flow, users must install the 
GRANULAR package with `PKG_GRANULAR=yes`.  Users should also include the 
OPENMP package with `PKG_OPENMP=yes`.

## Ghastly:

Ghastly can be installed from PyPi for general use.  If you plan to contribute
to Ghastly, you should instead create a fork of the arfc/ghastly repository, 
then clone and install Ghastly in editable mode:

```
git clone git@github.com:YourUsername/ghastly.git
cd ghastly/
python -m pip install --editable .
```

## Notes:

### Installing inside conda environments:
Installing inside a conda environment makes uninstalling Ghastly and its
compatible tools easier, which is useful if anything goes wrong (especially
if you would like to contribute to Ghastly).  Installing packages available
with conda-forge is simple.  For installing software from source into a conda
environment, remember:

- When a conda environment is active, the environment variable ${CONDA_PREFIX}
points to the directory in which the environment is installed.  This can be
used to direct your system to the install directories of dependencies and point
cmake to install into your conda environment.
    - `-D CMAKE_INSTALL_PREFIX=${CONDA_PREFIX}` will install to your conda
    environment
    - `-D CMAKE_PREFIX_PATH=${CONDA_PREFIX}` will help cmake find most packages
    installed with conda

