# MSPIN Package

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