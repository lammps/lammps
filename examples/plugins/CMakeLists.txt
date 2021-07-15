##########################################
# CMake build system for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.10)

# enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. You must create and use a build directory. "
    "Please remove CMakeCache.txt and CMakeFiles first.")
endif()

project(plugins VERSION 1.0 LANGUAGES CXX)

# NOTE: the next line should be commented out when used outside of the LAMMPS package
get_filename_component(LAMMPS_SOURCE_DIR ${PROJECT_SOURCE_DIR}/../../src ABSOLUTE)
set(LAMMPS_HEADER_DIR ${LAMMPS_SOURCE_DIR} CACHE PATH "Location of LAMMPS headers")
if(NOT LAMMPS_HEADER_DIR)
  message(FATAL_ERROR "Must set LAMMPS_HEADER_DIR")
endif()

# by default, install into $HOME/.local (not /usr/local),
# so that no root access (and sudo) is needed
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
endif()

# C++11 is required
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Need -restrict with Intel compilers
if((CMAKE_CXX_COMPILER_ID STREQUAL "Intel") OR (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM"))
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict")
endif()

# bail out on windows
if(CMAKE_SYSTEM_NAME STREQUAL Windows)
  message(FATAL_ERROR "LAMMPS plugins are currently not supported on Windows")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
include(CheckIncludeFileCXX)
include(LAMMPSInterfaceCXX)

##########################
# building the plugins

add_library(morse2plugin MODULE morse2plugin.cpp pair_morse2.cpp pair_morse2_omp.cpp)
target_include_directories(morse2plugin PRIVATE "${LAMMPS_HEADER_DIR}/OPENMP")
target_link_libraries(morse2plugin PRIVATE lammps)

add_library(nve2plugin MODULE nve2plugin.cpp fix_nve2.cpp)
target_link_libraries(nve2plugin PRIVATE lammps)

add_library(helloplugin MODULE helloplugin.cpp)
target_link_libraries(helloplugin PRIVATE lammps)

add_library(zero2plugin MODULE zero2plugin.cpp pair_zero2.cpp bond_zero2.cpp
                               angle_zero2.cpp dihedral_zero2.cpp improper_zero2.cpp)
target_link_libraries(zero2plugin PRIVATE lammps)

set_target_properties(morse2plugin nve2plugin helloplugin zero2plugin PROPERTIES
                      PREFIX ""
                      LINK_FLAGS "-rdynamic")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(morse2plugin nve2plugin helloplugin zero2plugin PROPERTIES
                        LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
endif()
