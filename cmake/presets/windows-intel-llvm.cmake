# preset that will enable Intel compilers with support for MPI and OpenMP (on Linux boxes)

set(CMAKE_CXX_COMPILER "icx" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "icx" CACHE STRING "" FORCE)
set(CMAKE_Fortran_COMPILER "ifx" CACHE STRING "" FORCE)
set(INTEL_LRT_MODE "C++11" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)
set(CMAKE_TUNE_FLAGS -Wno-unused-command-line-argument)
