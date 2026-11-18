# preset that will enable Intel compilers with support for MPI and OpenMP (on Linux and Windows boxes)

set(CMAKE_CXX_COMPILER "icx" CACHE STRING "" FORCE)
set(CMAKE_C_COMPILER "icx" CACHE STRING "" FORCE)
# prefer ifx over flang and gfortran, if available
find_program(INTEL_FORTRAN NAMES ifx flang gfortran)
if(INTEL_FORTRAN)
  set(CMAKE_Fortran_COMPILER ${INTEL_FORTRAN} CACHE STRING "" FORCE)
endif()
set(INTEL_LRT_MODE "c++17" CACHE STRING "" FORCE)
unset(HAVE_OMP_H_INCLUDE CACHE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-command-line-argument" CACHE STRING "" FORCE)
