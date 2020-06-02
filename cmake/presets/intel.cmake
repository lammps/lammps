# preset that will enable clang/clang++ with support for MPI and OpenMP (on Linux boxes)

set(CMAKE_CXX_COMPILER "icpc" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "icc" CACHE STRING "" FORCE)
set(MPI_CXX "icpc" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)

set(OpenMP_C "icc" CACHE STRING "" FORCE)
set(OpenMP_C_FLAGS "-qopenmp" CACHE STRING "" FORCE)
set(OpenMP_C_LIB_NAMES "omp" CACHE STRING "" FORCE)
set(OpenMP_CXX "icpc" CACHE STRING "" FORCE)
set(OpenMP_CXX_FLAGS "-qopenmp" CACHE STRING "" FORCE)
set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "" FORCE)
set(OpenMP_omp_LIBRARY "libiomp5.so" CACHE PATH "" FORCE)

