# Ghastly

Ghastly is a python package that externally couples Discrete Element Method
(DEM) physics simulation tools to neutron transport and depletion tools in
order to model flowing, pebble-bed High-Temperature Gas-Cooled Reactors (HTGRs).

Ghastly is compatible with [LAMMPS](https://www.lammps.org/) and 
[OpenMC](https://openmc.org/). Ghastly automatically generates input files and
can read output files from compatible codes.  The examples/ directory contains
Ghastly input files for a generic HTGR design.

Each template has some hard-coded settings, some of which are required to
create output files needed by Ghastly.  Users unfamiliar with compatible tools
should read the documentation carefully before making changes to templates.

Currently, Ghastly supports coupling LAMMPS to OpenMC only, but developing 
and submitting an update to add further code support via PR is welcome.

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
tools of choice inside a conda environment.

Activate the conda environment if you choose to use one before installing 
any dependencies. (See the Notes section below for tips about using conda
environments while installing from source with cmake)

## Dependencies and Compatible Tool installation:

### OpenMC:
The quick install guide for OpenMC is [here](https://docs.openmc.org/en/stable/quickinstall.html).
Follow the instructions to build from source.  The first step is to install HDF5.
The setup script included in this level of the repository is an example of how
to install mpi4py, parallel h5py, and parallel HDF5 inside a conda environment.


If you want to use OpenMC with OpenMP or MPI, set 
`-DOpenMC_USE_OPENMP=on` and `-DOPENMC_USE_MPI=on`.  If you are installing 
inside a conda environment, also set `-DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX}` 
and call `HDF5_ROOT=${CONDA_PREFIX}` before calling cmake.

### LAMMPS:
The documentation page for building LAMMPS begins [here](https://docs.lammps.org/Build.html) here.
LAMMPS has many optional packages and build options which will not be discussed
in depth here beyond what should be installed if working with Ghastly.  Users 
should set `-D BUILD_SHARED_LIBS=yes`, `-D BUILD_MPI=yes`, `-D BUILD_OMP=yes`.  
The setup script example also turns the JPEG, PNG, and FFMPEG options on.  
These are not needed for Ghastly to work, though users may prefer to have the option.
As always, setting `-D CMAKE_INSTALL_PREFIX` and `-D CMAKE_PREFIX_PATH` to
`${CONDA_PREFIX}` is required if using conda environments.  In order to use
LAMMPS to model pebble flow, users must install the GRANULAR package with 
`PKG_GRANULAR=yes`.  Users should also include the OPENMP package with 
`PKG_OPENMP=yes`.

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

- when a conda environment is active, the environment variable ${CONDA_PREFIX}
points to the directory in which the environment is installed.  This can be
used to direct your system to the install directories of dependencies and point
cmake to install into your conda environment.
    - `-D CMAKE_INSTALL_PREFIX=${CONDA_PREFIX}` will install to your conda
    environment
    - `-D CMAKE_PREFIX_PATH=${CONDA_PREFIX}` will help cmake find most packages
    installed with conda

