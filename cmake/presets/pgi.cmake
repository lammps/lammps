# preset that will enable PGI (Nvidia) compilers with support for MPI and OpenMP (on Linux boxes)

set(CMAKE_CXX_COMPILER "pgc++" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "pgcc" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "pgfortran" CACHE STRING "" FORCE)
set(MPI_CXX "pgc++" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)

