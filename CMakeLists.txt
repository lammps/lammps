cmake_minimum_required(VERSION 3.1)

project(lammps)
set(SOVERSION 0)

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
  file(GLOB MPI_SOURCES src/STUBS/mpi.c)
  include_directories(src/STUBS)
  set(MPI_CXX_LIBRARIES)
endif()

find_package(UnixCommands)

option(CMAKE_VERBOSE_MAKEFILE "Verbose makefile" OFF)

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

file(GLOB LIB_SOURCES src/*.cpp)
file(GLOB LMP_SOURCES src/main.cpp)
list(REMOVE_ITEM LIB_SOURCES ${LMP_SOURCES})

add_custom_target(style COMMAND ${BASH} Make.sh style WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
add_library(lammps ${LIB_SOURCES} ${MPI_SOURCES})
add_dependencies(lammps style)
# better but slower
#add_custom_command(TARGET lammps PRE_BUILD COMMAND ${BASH} Make.sh style WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(lammps ${MPI_CXX_LIBRARIES} ${MATH_LIBRARIES})
set_target_properties(lammps PROPERTIES SOVERSION ${SOVERSION})
install(TARGETS lammps LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(lmp ${LMP_SOURCES})
target_link_libraries(lmp lammps)
install(TARGETS lammps DESTINATION ${CMAKE_INSTALL_BINDIR})
