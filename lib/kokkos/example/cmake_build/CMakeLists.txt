# Kokkos requires CMake version 3.1 or higher and that you have the following
# line with a version of 3.1 or higher as the first line of your project:
#   cmake_minimum_required(VERSION 3.1)
#
# The other CMake commands required to build Kokkos as part of your application
# are:
#   add_subdirectory(path/to/kokkos)
#   target_link_libraries(executable or library)
#
# If Kokkos is not a subdirectory of your project, you will also need to pass a
# binary directory to add_subdirectory().  We had to pass the binary directory
# for this example for that reason.  Note that target_link_libraries() can be
# called on a target added by add_executable(), add_library(), or another
# similar command.
#
# All the flags, etc. required to build using the Kokkos library are
# transitively added to targets which depend on the library.
#
# The CMake variables CMAKE_CXX_STANDARD and CMAKE_CXX_EXTENSIONS are
# respected.  We recommend that you set CMAKE_CXX_EXTENSIONS to OFF.
# Otherwise, CMake defaults to using extensions for the C++ standard, and the
# GNU extensions (-std=gnu++11) will be used for compilers that support it
# instead of standard C++11 (-std=c++11).
#
# A bunch of build options are added as variables (all starting with KOKKOS_)
# to the build.  Check them out using ccmake or the CMake GUI.
#
# Building this example:
#   1. Create a build directory.
#   2. cd /path/to/build/directory
#   3. cmake /path/to/example
#   4. make

cmake_minimum_required(VERSION 3.1)
project(Example CXX C Fortran)

list(APPEND CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -O3)

add_subdirectory(${Example_SOURCE_DIR}/../.. ${Example_BINARY_DIR}/kokkos)

include_directories(${Kokkos_INCLUDE_DIRS_RET})

add_executable(example cmake_example.cpp foo.f)
target_link_libraries(example kokkos)
