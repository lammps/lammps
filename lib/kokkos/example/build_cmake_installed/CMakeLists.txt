# Kokkos minimally requires 3.16 right now,
# but your project can set it higher
cmake_minimum_required(VERSION 3.16)

# Projects can safely mix languages - must have C++ support
# Kokkos flags will only apply to C++ files
project(Example CXX Fortran)

# Look for an installed Kokkos
find_package(Kokkos REQUIRED)

add_executable(example cmake_example.cpp foo.f)
if(CMAKE_Fortran_COMPILER_ID STREQUAL LLVMFlang)
  set_target_properties(example PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_options(example PRIVATE -fno-fortran-main)
endif()

# This is the only thing required to set up compiler/linker flags
target_link_libraries(example Kokkos::kokkos)

enable_testing()
add_test(NAME KokkosInTree_Verify COMMAND example 10)
