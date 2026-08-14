# -*- CMake -*- build system for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.20)

project(mbxplugin VERSION 1.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(CheckIncludeFileCXX)
include(LAMMPSInterfacePlugin)
list(APPEND CMAKE_MODULE_PATH ${LAMMPS_SOURCE_DIR}/../cmake/Modules)
include(MBX)

##########################
# building the plugins

add_library(mbxplugin MODULE mbxplugin.cpp ${LAMMPS_SOURCE_DIR}/MBX/fix_mbx.cpp
  ${LAMMPS_SOURCE_DIR}/MBX/pair_mbx.cpp)
target_link_libraries(mbxplugin PRIVATE LAMMPS::MBX)
target_link_libraries(mbxplugin PRIVATE lammps)
target_include_directories(mbxplugin PRIVATE ${LAMMPS_SOURCE_DIR}/MBX)
set_target_properties(mbxplugin PROPERTIES PREFIX "" SUFFIX ".so")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(mbxplugin PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
# tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
  set_target_properties(mbxplugin PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(mbxplugin  PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
  endif()

  get_lammps_version(${LAMMPS_SOURCE_DIR}/version.h LAMMPS_VERSION)
  find_program(MAKENSIS_PATH makensis)
  if(MAKENSIS_PATH)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/lammps.ico
      ${CMAKE_SOURCE_DIR}/lammps-text-logo-wide.bmp ${CMAKE_SOURCE_DIR}/mbxplugin.nsis
      ${CMAKE_BINARY_DIR})
    if(BUILD_MPI)
      if(USE_MSMPI)
        add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION}-MSMPI mbxplugin.nsis
          DEPENDS mbxplugin lammps.ico lammps-text-logo-wide.bmp mbxplugin.nsis
          BYPRODUCTS LAMMPS-MBX-plugin-${LAMMPS_VERSION}-MSMPI.exe)
      else()
        add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION}-MPI mbxplugin.nsis
          DEPENDS mbxplugin lammps.ico lammps-text-logo-wide.bmp mbxplugin.nsis
          BYPRODUCTS LAMMPS-MBX-plugin-${LAMMPS_VERSION}-MPI.exe)
      endif()
    else()
      add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION} mbxplugin.nsis
        COMMAND ${CMAKE_COMMAND} -E echo ${PWD}
        DEPENDS mbxplugin lammps.ico lammps-text-logo-wide.bmp mbxplugin.nsis
        BYPRODUCTS LAMMPS-MBX-plugin-${LAMMPS_VERSION}.exe)
    endif()
  endif()
else()
  set_target_properties(mbxplugin PROPERTIES LINK_FLAGS "-rdynamic")
endif()
