cmake_minimum_required(VERSION 3.10)
project(simple CXX)

set(LAMMPS_SRC_DIRECTORY "" CACHE PATH "Path for lammps source")
if(NOT LAMMPS_SRC_DIRECTORY STREQUAL "" AND EXISTS ${LAMMPS_SRC_DIRECTORY}/cmake/CMakeLists.txt)
  option(BUILD_LIB "Build LAMMPS library" ON)
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${LAMMPS_SRC_DIRECTORY}/src ${LAMMPS_SRC_DIRECTORY}/src/lammps)
  add_subdirectory(${LAMMPS_SRC_DIRECTORY}/cmake lammps)
else()
  find_package(LAMMPS REQUIRED)
endif()

add_executable(simple simple.cpp)
target_link_libraries(simple LAMMPS::lammps MPI::MPI_CXX)
