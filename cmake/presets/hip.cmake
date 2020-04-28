# preset that will enable hipcc plus gcc with support for MPI and OpenMP (on Linux boxes)

set(CMAKE_CXX_COMPILER "hipcc" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "gcc" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -O2 -DNDEBG" CACHE STRING "" FORCE)
set(MPI_CXX "hipcc" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)

set(OpenMP_CXX "hipcc" CACHE STRING "" FORCE)
set(OpenMP_CXX_FLAGS "-fopenmp" CACHE STRING "" FORCE)
set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
set(OpenMP_omp_LIBRARY "libomp.so" CACHE PATH "" FORCE)

