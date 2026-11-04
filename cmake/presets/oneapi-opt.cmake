# preset that will enable Intel compilers with support for MPI and OpenMP (on Linux and Windows boxes)
# like oneapi.cmake, but with more aggressive optimization: -O3 and code for the CPU of the
# machine where LAMMPS is compiled (-xHost), which also lets the INTEL package styles use
# AVX-512 instructions where available.  LAMMPS compiled this way is typically about 10%
# faster, but may not run on a different kind of CPU, and its results differ from those with
# the default flags in the last digits.  Please see the "Compiler optimization" section in
# the manual.

include(${CMAKE_CURRENT_LIST_DIR}/oneapi.cmake)

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3 -xHost -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O3 -xHost -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -DNDEBUG" CACHE STRING "" FORCE)
