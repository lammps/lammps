# Kokkos minimally requires 3.16 right now,
# but your project can set it higher
cmake_minimum_required(VERSION 3.16)

# Projects can safely mix languages - must have C++ support
# Kokkos flags will only apply to C++ files
project(Example CXX Fortran)

# You need this for using Kokkos_ROOT variable
message(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
cmake_policy(SET CMP0074 NEW)

# Look for an installed Kokkos but force using the compiler launcher
# to ensure that targets depending on Kokkos use the same compiler
# as when kokkos was installed, e.g. if kokkos was built with
# g++ and the CMAKE_CXX_COMPILER=clang++ then example_with_kokkos
# will be compiled and linked with g++ whereas example_no_kokkos
# will be compiled and linked with clang++
find_package(Kokkos REQUIRED COMPONENTS launch_compiler)

add_executable(example_no_kokkos bar.cpp)
add_executable(example_with_kokkos foo.cpp)

# This is the only thing required to set up compiler/linker flags
target_link_libraries(example_with_kokkos Kokkos::kokkos)

enable_testing()
add_test(NAME KokkosLauncher_NoKokkos_Verify COMMAND example_no_kokkos 10)
add_test(NAME KokkosLauncher_WithKokkos_Verify COMMAND example_with_kokkos 10)
