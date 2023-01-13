cmake_minimum_required(VERSION 3.10)

# enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. You must create and use a build directory. "
    "Please remove CMakeCache.txt and CMakeFiles first.")
endif()

project(couple-simple VERSION 1.0 LANGUAGES C CXX)

# by default, install into $HOME/.local (not /usr/local),
# so that no root access (and sudo) is needed
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
endif()

# ugly hacks for MSVC which by default always reports an old C++ standard in the __cplusplus macro
# and prints lots of pointless warnings about "unsafe" functions
if(MSVC)
  if(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    if(LAMMPS_EXCEPTIONS)
      add_compile_options(/EHsc)
    endif()
  endif()
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

find_package(MPI QUIET)
# do not include the (obsolete) MPI C++ bindings which makes
# for leaner object files and avoids namespace conflicts
set(MPI_CXX_SKIP_MPICXX TRUE)

##########################

# build within LAMMPS build system
if(NOT LAMMPS_SOURCE_DIR)
  find_package(LAMMPS REQUIRED)
endif()

add_executable(simpleCC simple.cpp)
target_link_libraries(simpleCC PRIVATE LAMMPS::lammps MPI::MPI_CXX)

add_executable(simpleC simple.c)
target_link_libraries(simpleC PRIVATE LAMMPS::lammps MPI::MPI_C)

if(LAMMPS_SOURCE_DIR)
  include(CheckGeneratorSupport)
  if(NOT CMAKE_GENERATOR_SUPPORT_FORTRAN)
    message(STATUS "Skipping Test for the LAMMPS Fortran Module Coupling: no Fortran support in build tool")
    return()
  endif()

  include(CheckLanguage)
  check_language(Fortran)

  if(CMAKE_Fortran_COMPILER AND MPI_Fortran_FOUND AND MPI_Fortran_HAVE_F90_MODULE)
    enable_language(Fortran)

    if(NOT CMAKE_Fortran_COMPILER_ID)
      message(STATUS "Skipping Tests for the LAMMPS Fortran Module: cannot identify Fortran compiler")
      return()
    endif()

    # GNU Fortran 4.8.x on RHEL/CentOS 7.x is not sufficient to compile the Fortran module.
    # Work around flang being detected as GNU
    get_filename_component(_tmp_fc ${CMAKE_Fortran_COMPILER} NAME)
    if((CMAKE_Fortran_COMPILER_ID STREQUAL "GNU") AND (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 6.0) AND NOT (_tmp_fc STREQUAL "flang"))
      message(FATAL_ERROR "Need GNU Fortran compiler version 6.x or later for unit testing")
    endif()

    get_filename_component(LAMMPS_FORTRAN_MODULE ${LAMMPS_SOURCE_DIR}/../fortran/lammps.f90 ABSOLUTE)
    add_executable(simpleF90 ${LAMMPS_FORTRAN_MODULE} simple.f90)
    target_link_libraries(simpleF90 PRIVATE MPI::MPI_Fortran LAMMPS::lammps)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      target_compile_options(simpleF90 PRIVATE -fallow-argument-mismatch)
    endif()

    get_filename_component(_tmp_fc ${CMAKE_Fortran_COMPILER} NAME)
    if (_tmp_fc STREQUAL "flang")
      target_link_libraries(simpleF90 PRIVATE gfortran)
    endif()
  endif()
endif()
