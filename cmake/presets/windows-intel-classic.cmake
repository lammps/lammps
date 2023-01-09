# preset that will enable Intel compilers with support for MPI and OpenMP (on Linux boxes)

set(CMAKE_CXX_COMPILER "icl" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "icl" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "ifort" CACHE STRING "" FORCE)

unset(HAVE_OMP_H_INCLUDE CACHE)

