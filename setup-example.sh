#!/bin/bash

#the following is an example setup script for installing ghastly, openmc, and
#LAMMPS into a conda environment.  It will need to be edited for your system.

#before running this script, you should have all the zip files/tarballs/etc
#this references downloaded and extracted, and you will need build essentials, conda, git,
# cmake, and a conda env called ghastly-dev created

#you need to download:
#hdf5 zip

source activate base
conda activate ghastly-dev
conda config --add channels conda-forge
conda install mpich
conda install python
conda install numpy
conda install scipy
conda install matplotlib
conda install cython
conda install libmpich-dev
conda install mpich-mpicc mpich-mpicxx mpich-mpifort
conda install zlib
conda install libpng
conda install vtk
conda install pytest
conda install ffmpeg
conda install libjpeg-turbo
conda install uncertainties
conda install lxml
cd hdf5_version/ # update this line with the correct hdf5 directory name
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DZLIB_LIBRARY:FILEPATH=${CONDA_PREFIX}/lib/libz.so -DZLIB_INCLUDE_DIR:PATH=${CONDA_PREFIX}/include -DZLIB_USE_EXTERNAL=OFF -DHDF5_ENABLE_SZIP_SUPPORT=OFF ..
cd ..
cmake --build ./build --config Release -j 4
cmake --install ./build --config Release
conda deactivate
conda activate ghastly-dev
h5dump --version
python -m pip --version
python -m pip install mpi4py
cc=${CONDA_PREFIX} HDF5_MPI=ON HDF5_DIR=${CONDA_PREFIX} python -m pip install --no-binary=h5py h5py
python -c "import h5py; print(h5py.version.hdf5_version)"
git clone --recurse-submodules https://github.com/openmc-dev/openmc.git
cd ~/openmc/
mkdir build && cd build
HDF5_ROOT=${CONDA_PREFIX} cmake -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX} -DOPENMC_USE_OPENMP=on -DOPENMC_USE_MPI=on ..
make
make install
cd ..
python -m pip install -e .
cd ~
git clone -b release https://github.com/lammps/lammps.git lammps_install
cd lammps_install/
mkdir build && cd build
cmake -C ../cmake/presets/most.cmake -D CMAKE_INSTALL_PREFIX=${CONDA_PREFIX} -D CMAKE_PREFIX_PATH=${CONDA_PREFIX} -D BUILD_SHARED_LIBS=yes -D BUILD_MPI=yes -D BUILD_OMP=yes -D WITH_JPEG=yes -D WITH_PNG=yes -D WITH_FFMPEG=yes -D PKG_GRANULAR=yes -D PKG_PYTHON=yes -D PKG_OMP=yes -D PKG_INTEL=yes -D PKG_OPT=yes ../cmake
make
make install
cd ~
git clone git@github.com:arfc/ghastly.git
cd ghastly/
python -m pip install --editable .
echo "Done!"
