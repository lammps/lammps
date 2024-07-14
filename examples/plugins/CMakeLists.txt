# -*- CMake -*- file for plugin examples.
# The is meant to be used as a template for plugins that are
# distributed independent from the LAMMPS package.
##########################################

cmake_minimum_required(VERSION 3.16)

# enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. You must create and use a build directory. "
    "Please remove CMakeCache.txt and CMakeFiles first.")
endif()

project(plugins VERSION 1.0 LANGUAGES CXX)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
if(NOT LAMMPS_DIR)
  include(LAMMPSBuildMPI)
endif()

find_package(LAMMPS REQUIRED)

##########################
# building the plugins

lammps_add_plugin(morse2plugin)
target_sources(morse2plugin PRIVATE morse2plugin.cpp pair_morse2.cpp pair_morse2_omp.cpp)

lammps_add_plugin(nve2plugin)
target_sources(nve2plugin PRIVATE nve2plugin.cpp fix_nve2.cpp)

lammps_add_plugin(helloplugin)
target_sources(helloplugin PRIVATE helloplugin.cpp)


lammps_add_plugin(zero2plugin)
target_sources(zero2plugin PRIVATE zero2plugin.cpp pair_zero2.cpp bond_zero2.cpp
                                   angle_zero2.cpp dihedral_zero2.cpp improper_zero2.cpp)

add_custom_target(plugins ALL ${CMAKE_COMMAND} -E echo "Building Plugins"
  DEPENDS morse2plugin nve2plugin helloplugin zero2plugin morse2plugin)

include(GNUInstallDirs)
install(TARGETS morse2plugin nve2plugin helloplugin zero2plugin morse2plugin
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/${LAMMPS_PLUGIN_DIR})
