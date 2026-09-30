cmake_minimum_required(VERSION 3.28.2)

project(Example CXX)

# Look for an installed Kokkos
find_package(Kokkos REQUIRED)

add_executable(example cmake_example.cpp)

# This is the only thing required to set up compiler/linker flags
# FIXME_MODULES We should be able to link against Kokkos::kokkos
target_link_libraries(example Kokkos::kokkoscore)

enable_testing()
add_test(NAME KokkosInTree_Verify COMMAND example 10)
