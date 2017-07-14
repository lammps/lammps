cmake_minimum_required(VERSION 3.1)

project(lammps)
set(SOVERSION 0)
set(LAMMPS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../src)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_C_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_C_FLAGS)

enable_language(CXX)

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCCompilerFlag)

########################################################################
# User input options                                                   #
########################################################################
option(BUILD_SHARED_LIBS "Build shared libs" OFF)
include(GNUInstallDirs)

option(ENABLE_MPI "Build MPI version" OFF)
if(ENABLE_MPI)
  find_package(MPI)
  include_directories(${MPI_C_INCLUDE_PATH})
  set(MPI_SOURCES)
else()
  file(GLOB MPI_SOURCES ${LAMMPS_SOURCE_DIR}/STUBS/mpi.c)
  include_directories(${LAMMPS_SOURCE_DIR}/STUBS)
  set(MPI_CXX_LIBRARIES)
endif()

find_package(UnixCommands)

option(CMAKE_VERBOSE_MAKEFILE "Verbose makefile" OFF)

set(PACKAGES ASPHERE REAX)
foreach(PKG ${PACKAGES})
  option(ENABLE_${PKG} "Build ${PKG} Package" OFF)
endforeach()

find_package(JPEG)
if(JPEG_FOUND)
  add_definitions(-DLAMMPS_JPEG)
  include_directories(${JPEG_INCLUDE_DIR})
else()
  set(JPEG_LIBRARIES)
endif()

find_package(PNG)
if(PNG_FOUND)
  include_directories(${PNG_INCLUDE_DIR})
  add_definitions(-DLAMMPS_PNG)
else(PNG_FOUND)
  set(PNG_LIBRARIES)
endif(PNG_FOUND)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFile)
foreach(HEADER math.h)
  check_include_file(${HEADER} FOUND_${HEADER})
  if(NOT FOUND_${HEADER})
    message(FATAL_ERROR "Could not find needed header - ${HEADER}")
  endif(NOT FOUND_${HEADER})
endforeach(HEADER)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
include(CheckLibraryExists)
foreach(FUNC sin cos)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)

######################################
# Include the following subdirectory # 
######################################

#Do NOT go into src to not conflict with old Makefile build system
#add_subdirectory(src)

file(GLOB LIB_SOURCES ${LAMMPS_SOURCE_DIR}/*.cpp)
file(GLOB LMP_SOURCES ${LAMMPS_SOURCE_DIR}/main.cpp)
list(REMOVE_ITEM LIB_SOURCES ${LMP_SOURCES})

foreach(PKG ${PACKAGES})
  if(ENABLE_${PKG})
    file(GLOB ${PKG}_SOURCES ${LAMMPS_SOURCE_DIR}/${PKG}/*.cpp)
    list(APPEND LIB_SOURCES ${${PKG}_SOURCES})
    include_directories(${LAMMPS_SOURCE_DIR}/${PKG})
  endif()
endforeach()

if(ENABLE_REAX)
  enable_language(Fortran)
  file(GLOB REAX_SOURCES ${CMAKE_SOURCE_DIR}/../lib/reax/*.F)
  list(APPEND LIB_SOURCES ${REAX_SOURCES})
  include_directories(${CMAKE_SOURCE_DIR}/../lib/reax)
endif()
include_directories(${LAMMPS_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Headers)
configure_file(${CMAKE_SOURCE_DIR}/Headers/package.h.cmakein ${CMAKE_BINARY_DIR}/cmake/package.h)
include_directories(${CMAKE_BINARY_DIR}/cmake)

add_library(lammps ${LIB_SOURCES} ${MPI_SOURCES})
target_link_libraries(lammps ${MPI_CXX_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${MATH_LIBRARIES})
set_target_properties(lammps PROPERTIES SOVERSION ${SOVERSION})
install(TARGETS lammps LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(lmp ${LMP_SOURCES})
target_link_libraries(lmp lammps)
install(TARGETS lammps DESTINATION ${CMAKE_INSTALL_BINDIR})
