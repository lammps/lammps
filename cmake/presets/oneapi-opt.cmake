# preset that will enable Intel compilers with support for MPI and OpenMP (on Linux and Windows boxes)
# like oneapi.cmake, but with more aggressive optimization: -O3 and code for the x86-64-v3
# instruction set (AVX2 with fused multiply-add, available since Intel Haswell and AMD Zen).
# LAMMPS compiled this way is typically about 10% faster, but does not run on older CPUs,
# and its results differ from those with the default flags in the last digits.  Since a
# target CPU is selected, the INTEL package does not add -xHost and its styles will not use
# AVX-512 instructions.  Please see the "Compiler optimization" section in the manual.

include(${CMAKE_CURRENT_LIST_DIR}/oneapi.cmake)

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3 -march=x86-64-v3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=x86-64-v3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -O3 -march=x86-64-v3 -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=x86-64-v3 -DNDEBUG" CACHE STRING "" FORCE)
