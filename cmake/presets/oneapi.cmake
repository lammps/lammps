# preset that will enable the LLVM based Intel compilers with support for MPI and OpenMP and Fortran (on Linux boxes)

set(CMAKE_CXX_COMPILER "icpx" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "icx" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "ifx" CACHE STRING "" FORCE)
set(MPI_CXX "icpx" CACHE STRING "" FORCE)
set(MPI_CXX_COMPILER "mpicxx" CACHE STRING "" FORCE)

# force using internal BLAS/LAPCK since external ones may not be ABI compatible
set(USE_INTERNAL_LINALG ON CACHE BOOL "" FORCE)

