# -*- CMake -*- configuration for building the fix symmetry helper tool: find-symmetry

cmake_minimum_required(VERSION 3.27)

# set timestamp of downloaded files to that of archive
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

project(symmetry
  DESCRIPTION "Space group operators for fix symmetry"
  LANGUAGES CXX C)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

if(MSVC)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

# standalone build: find the LAMMPS CMake modules
if(NOT LAMMPS_DIR)
  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules)
  set(LAMMPS_THIRDPARTY_URL "https://download.lammps.org/thirdparty")
endif()
# the tool requires spglib, which phana can be told not to use
if(DEFINED USE_SPGLIB AND NOT USE_SPGLIB)
  message(STATUS "Skipping the find-symmetry tool since USE_SPGLIB is off")
  return()
endif()

include(GNUInstallDirs)
include(BuildSpglib)

add_executable(find-symmetry find-symmetry.cpp)
target_compile_features(find-symmetry PRIVATE cxx_std_17)
target_link_libraries(find-symmetry PRIVATE SPGLIB::SYMSPG)
install(TARGETS find-symmetry DESTINATION ${CMAKE_INSTALL_BINDIR})
