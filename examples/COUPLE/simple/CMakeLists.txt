cmake_minimum_required(VERSION 3.10)
project(simple CXX)

set(LAMMPS_SRC_DIRECTORY "" CACHE PATH "Path for lammps source")
if(NOT LAMMPS_SRC_DIRECTORY STREQUAL "" AND EXISTS ${LAMMPS_SRC_DIRECTORY}/cmake/CMakeLists.txt)
  option(BUILD_LIB "Build LAMMPS library" ON)
  add_subdirectory(${LAMMPS_SRC_DIRECTORY}/cmake lammps)
else()
  find_package(LAMMPS REQUIRED)
endif()

add_executable(simpleCC simple.cpp)
target_link_libraries(simpleCC LAMMPS::lammps)

enable_language(C)
add_executable(simpleC simple.c)
target_link_libraries(simpleC LAMMPS::lammps)
