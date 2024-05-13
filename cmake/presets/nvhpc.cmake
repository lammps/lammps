# preset that will enable Nvidia HPC SDK compilers with support for MPI and OpenMP (on Linux boxes)

set(CMAKE_CXX_COMPILER "nvc++" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "nvc" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "nvfortran" CACHE STRING "" FORCE)
set(MPI_CXX "nvc++" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)

