cmake_minimum_required(VERSION 3.10)
project(simple CXX)

set(LAMMPS_SRC_DIRECTORY "" CACHE PATH "Path for lammps source")
if(NOT LAMMPS_SRC_DIRECTORY STREQUAL "" AND EXISTS ${LAMMPS_SRC_DIRECTORY}/cmake/CMakeLists.txt)
  option(BUILD_LIB "Build LAMMPS library" ON)
  option(LAMMPS_CREATE_IN_SOURCES_INCLUDE_DIR "include in-source lammps/ include dir" ON)
  add_subdirectory(${LAMMPS_SRC_DIRECTORY}/cmake lammps)
else()
  find_package(LAMMPS REQUIRED)
endif()

add_executable(simple simple.cpp)
target_link_libraries(simple LAMMPS::lammps)
