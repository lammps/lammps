# preset that will enable clang/clang++ with support for MPI and OpenMP (on Linux boxes)
# like clang.cmake, but with more aggressive optimization: -O3 and on x86_64 CPUs code for
# the x86-64-v3 instruction set (AVX2 with fused multiply-add, available since Intel
# Haswell and AMD Zen).  LAMMPS compiled this way is typically 5-10% faster, but does not
# run on older CPUs, and its results differ from those with the default flags in the last
# digits.  Please see the "Compiler optimization" section in the manual.

include(${CMAKE_CURRENT_LIST_DIR}/clang.cmake)

cmake_host_system_information(RESULT _host_platform QUERY OS_PLATFORM)
if(_host_platform MATCHES "^(x86_64|AMD64|amd64)$")
  set(_opt_flags "-O3 -march=x86-64-v3")
else()
  set(_opt_flags "-O3")
endif()

set(_warn_flags "-Wall -Wextra -Wno-bitwise-instead-of-logical")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${_warn_flags} -g ${_opt_flags} -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "${_opt_flags} -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${_warn_flags} -g ${_opt_flags} -DNDEBUG" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "${_opt_flags} -DNDEBUG" CACHE STRING "" FORCE)
unset(_warn_flags)
unset(_opt_flags)
unset(_host_platform)
