# MSPIN Package
This package contains a `fix rigid/nvt/mspin` command that updates nanoparticle
dynamics subjected to external magnetic field and mangetic dipolar interactions.

It also contains commands to compute the externel field interaction energy,
dipolar interaction energy, and interparticle distance during simulation.

See the doc page for the `fix rigid/nvt/mspin` or the `compute mspin/energy`
or `compute mspin/distance` commands for detailed usage instructions.

Use of this package requires LAMMPS to be built with the RIGID package.

There are example scripts for using commands in this package in the
examples/mspin directory.

The authors of the package is Akhlak U. Mahmood (amahmoo3 at ncsu dot edu)
and Yaroslava G. Yingling (yara_yingling at ncsu dot edu) at North Carolina
State University, USA. Contact the authors directly if you have questions.

# Installation
CMake based installation preset.

```cmake
# Enable required packages
set(ALL_PACKAGES KSPACE MOLECULE RIGID MSPIN)
foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

# Update as necessary
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/mspin" CACHE PATH "Default install path" FORCE)
set(LAMMPS_MACHINE serial CACHE STRING "" FORCE)

# Turn on MPI support
# Make sure you installed openmpi or mpich
# apt-get install libopenmpi-dev
# set(MPI_CXX "icpx" CACHE STRING "" FORCE)
# set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
# set(BUILD_MPI ON CACHE BOOL "" FORCE)
# set(LAMMPS_MACHINE mpi CACHE STRING "" FORCE)
```
