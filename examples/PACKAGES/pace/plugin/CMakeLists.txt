##########################################
# CMake build system for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.10)

project(paceplugin VERSION 1.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(CheckIncludeFileCXX)
include(LAMMPSInterfacePlugin)
include(ML-PACE)

##########################
# building the plugins

add_library(paceplugin MODULE paceplugin.cpp ${LAMMPS_SOURCE_DIR}/ML-PACE/pair_pace.cpp)
target_link_libraries(paceplugin PRIVATE pace)
target_link_libraries(paceplugin PRIVATE lammps)
target_include_directories(paceplugin PRIVATE ${LAMMPS_SOURCE_DIR}/ML-PACE)
set_target_properties(paceplugin PROPERTIES PREFIX "" SUFFIX ".so")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(paceplugin PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
# tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
  set_target_properties(paceplugin PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(paceplugin  PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
  endif()

  get_lammps_version(${LAMMPS_SOURCE_DIR}/version.h LAMMPS_VERSION)
  find_program(MAKENSIS_PATH makensis)
  if(MAKENSIS_PATH)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/lammps.ico
      ${CMAKE_SOURCE_DIR}/lammps-text-logo-wide.bmp ${CMAKE_SOURCE_DIR}/paceplugin.nsis ${CMAKE_BINARY_DIR})
    if(BUILD_MPI)
      add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION}-MPI paceplugin.nsis
        DEPENDS paceplugin
        BYPRODUCTS LAMMPS-ML-PACE-plugin-${LAMMPS_VERSION}-MPI.exe)
    else()
      add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION} paceplugin.nsis
        COMMAND ${CMAKE_COMMAND} -E echo ${PWD}
        DEPENDS paceplugin lammps.ico lammps-text-logo-wide.bmp paceplugin.nsis
        BYPRODUCTS LAMMPS-ML-PACE-plugin-${LAMMPS_VERSION}.exe)
    endif()
  endif()
else()
  set_target_properties(paceplugin PROPERTIES LINK_FLAGS "-rdynamic")
endif()
