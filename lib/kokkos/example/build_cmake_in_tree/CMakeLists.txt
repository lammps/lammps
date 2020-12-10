# Kokkos minimally requires 3.10 right now,
# but your project can set it higher
cmake_minimum_required(VERSION 3.10)

# Project can mix languages - must have C++ support
# Kokkos flags are only applied to C++ files
project(Example CXX)

# We build kokkos as a subdirectory of our project
add_subdirectory(${Example_SOURCE_DIR}/../.. ${Example_BINARY_DIR}/kokkos)

add_executable(example cmake_example.cpp)

# This is the only line required to set up all compiler/linker flags
target_link_libraries(example Kokkos::kokkos)

# Adds a test for the executable
enable_testing()
add_test(NAME KokkosInTree_Verify COMMAND example 10)
