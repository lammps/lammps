# -*- CMake -*- build system for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.16)

project(plumedplugin VERSION 1.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(CheckIncludeFileCXX)
include(LAMMPSInterfacePlugin)
include(PLUMED)

##########################
# building the plugins

add_library(plumedplugin MODULE plumedplugin.cpp ${LAMMPS_SOURCE_DIR}/PLUMED/fix_plumed.cpp)
target_link_libraries(plumedplugin PRIVATE LAMMPS::PLUMED)
target_link_libraries(plumedplugin PRIVATE lammps)
target_include_directories(plumedplugin PRIVATE ${LAMMPS_SOURCE_DIR}/PLUMED)
set_target_properties(plumedplugin PROPERTIES PREFIX "" SUFFIX ".so")

# MacOS seems to need this
if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set_target_properties(plumedplugin PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
# tell CMake to export all symbols to a .dll on Windows with special case for MinGW cross-compilers
  set_target_properties(plumedplugin PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
  if(CMAKE_CROSSCOMPILING)
    set_target_properties(plumedplugin  PROPERTIES LINK_FLAGS "-Wl,--export-all-symbols")
  endif()

  get_lammps_version(${LAMMPS_SOURCE_DIR}/version.h LAMMPS_VERSION)
  find_program(MAKENSIS_PATH makensis)
  if(MAKENSIS_PATH)
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_SOURCE_DIR}/lammps.ico
      ${CMAKE_SOURCE_DIR}/lammps-text-logo-wide.bmp ${CMAKE_SOURCE_DIR}/plumedplugin.nsis
      ${CMAKE_BINARY_DIR})
    if(BUILD_MPI)
      if(USE_MSMPI)
        add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION}-MSMPI plumedplugin.nsis
          DEPENDS plumedplugin plumed_copy lammps.ico lammps-text-logo-wide.bmp plumedplugin.nsis
          BYPRODUCTS LAMMPS-PLUMED-plugin-${LAMMPS_VERSION}-MSMPI.exe)
      else()
        add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION}-MPI plumedplugin.nsis
          DEPENDS plumedplugin plumed_copy lammps.ico lammps-text-logo-wide.bmp plumedplugin.nsis
          BYPRODUCTS LAMMPS-PLUMED-plugin-${LAMMPS_VERSION}-MPI.exe)
      endif()
    else()
      add_custom_target(package ${MAKENSIS_PATH} -V1 -DVERSION=${LAMMPS_VERSION} plumedplugin.nsis
        COMMAND ${CMAKE_COMMAND} -E echo ${PWD}
        DEPENDS plumedplugin plumed_copy lammps.ico lammps-text-logo-wide.bmp plumedplugin.nsis
        BYPRODUCTS LAMMPS-PLUMED-plugin-${LAMMPS_VERSION}.exe)
    endif()
  endif()
else()
  set_target_properties(plumedplugin PROPERTIES LINK_FLAGS "-rdynamic")
endif()
