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

# when this file is included as subdirectory in the LAMMPS build, many settings are directly imported
if(LAMMPS_DIR)
  set(LAMMPS_HEADER_DIR ${LAMMPS_SOURCE_DIR})
else()
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
  # ugly hacks for MSVC which by default always reports an old C++ standard in the __cplusplus macro
  # and prints lots of pointless warnings about "unsafe" functions
  if(MSVC)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
      add_compile_options(/Zc:__cplusplus)
      add_compile_options(/wd4244)
      add_compile_options(/wd4267)
      if(LAMMPS_EXCEPTIONS)
        add_compile_options(/EHsc)
      endif()
    endif()
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
  endif()
endif()

# C++11 is required
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Need -restrict with Intel compilers
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(CheckIncludeFileCXX)
if(NOT LAMMPS_DIR)
  include(LAMMPSInterfaceCXX)
endif()

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

set_target_properties(morse2plugin nve2plugin helloplugin zero2plugin PROPERTIES PREFIX "" SUFFIX ".so")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(morse2plugin nve2plugin helloplugin zero2plugin
    PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
# tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
  set_target_properties(morse2plugin nve2plugin helloplugin zero2plugin
    PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(morse2plugin nve2plugin helloplugin zero2plugin
      PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
  endif()
else()
  set_target_properties(morse2plugin nve2plugin helloplugin zero2plugin PROPERTIES
    LINK_FLAGS "-rdynamic")
endif()

add_custom_target(plugins ALL ${CMAKE_COMMAND} -E echo "Building Plugins"
  DEPENDS morse2plugin nve2plugin helloplugin zero2plugin morse2plugin)
