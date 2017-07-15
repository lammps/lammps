cmake_minimum_required(VERSION 3.1)

project(lammps)
set(SOVERSION 0)
set(LAMMPS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../src)
set(LAMMPS_LIB_SOURCE_DIR ${CMAKE_SOURCE_DIR}/../lib)
set(LAMMPS_LIB_BINARY_DIR ${CMAKE_BINARY_DIR}/lib)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/Modules)

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
option(INSTALL_LIB "Install lammps library and header" ON)
include(GNUInstallDirs)

set(LAMMPS_LINK_LIBS)
option(ENABLE_MPI "Build MPI version" OFF)
if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_C_INCLUDE_PATH})
  list(APPEND LAMMPS_LINK_LIBS ${MPI_CXX_LIBRARIES})
endif()

find_package(UnixCommands)

option(CMAKE_VERBOSE_MAKEFILE "Verbose makefile" OFF)

option(ENABLE_ALL "Build all packages" OFF)
set(PACKAGES ASPHERE BODY CLASS2 COLLOID COMPRESS CORESHELL DIPOLE GRANULAR
	KSPACE MANYBODY MC MEAM MISC MOLECULE RIGID REAX)
foreach(PKG ${PACKAGES})
  option(ENABLE_${PKG} "Build ${PKG} Package" ${ENABLE_ALL})
endforeach()

set(ACCEL_PACKAGES USER_OMP KOKKOS)
foreach(PKG ${ACCEL_PACKAGES})
  option(ENABLE_${PKG} "Build ${PKG} Package" OFF)
endforeach()

if(ENABLE_KSPACE)
  find_package(FFTW3)
  if(FFTW3_FOUND)
    add_definitions(-DFFT_FFTW3)
    include_directories(${FFTW3_INCLUDE_DIRS})
    list(APPEND LAMMPS_LINK_LIBS ${FFTW3_LIBRARIES})
  endif()
endif()

if(ENABLE_KOKKOS)
  # starting with CMake 3.1 this is all you have to do to enforce C++11
  set (CMAKE_CXX_STANDARD 11)
  set(LAMMPS_LIB_KOKKOS_SRC_DIR ${LAMMPS_LIB_SOURCE_DIR}/kokkos)
  set(LAMMPS_LIB_KOKKOS_BIN_DIR ${LAMMPS_LIB_BINARY_DIR}/kokkos)
  add_definitions(-DLMP_KOKKOS)
  add_subdirectory(${LAMMPS_LIB_KOKKOS_SRC_DIR} ${LAMMPS_LIB_KOKKOS_BIN_DIR})
  message("KOKKOS_DIRS: ${Kokkos_INCLUDE_DIRS}")

  # TODO there probably is a better way
  set(Kokkos_INCLUDE_DIRS ${LAMMPS_LIB_KOKKOS_SRC_DIR}/core/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/containers/src
                          ${LAMMPS_LIB_KOKKOS_SRC_DIR}/algorithms/src
                          ${LAMMPS_LIB_KOKKOS_BIN_DIR})
  include_directories(${Kokkos_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${Kokkos_LIBRARIES})
endif()

find_package(JPEG)
if(JPEG_FOUND)
  add_definitions(-DLAMMPS_JPEG)
  include_directories(${JPEG_INCLUDE_DIR})
  list(APPEND LAMMPS_LINK_LIBS ${JPEG_LIBRARIES})
endif()

find_package(PNG)
find_package(ZLIB)
if(PNG_FOUND AND ZLIB_FOUND)
  include_directories(${PNG_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${PNG_LIBRARIES} ${ZLIB_LIBRARIES})
  add_definitions(-DLAMMPS_PNG)
endif()

find_program(GZIP gzip)
if(GZIP)
  add_definitions(-DLAMMPS_GZIP)
endif()

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
list(APPEND LAMMPS_LINK_LIBS ${MATH_LIBRARIES})

######################################
# Include the following subdirectory # 
######################################

#Do NOT go into src to not conflict with old Makefile build system
#add_subdirectory(src)

file(GLOB LIB_SOURCES ${LAMMPS_SOURCE_DIR}/*.cpp)
file(GLOB LMP_SOURCES ${LAMMPS_SOURCE_DIR}/main.cpp)
list(REMOVE_ITEM LIB_SOURCES ${LMP_SOURCES})

if(NOT ENABLE_MPI)
  file(GLOB MPI_SOURCES ${LAMMPS_SOURCE_DIR}/STUBS/mpi.c)
  list(APPEND LIB_SOURCES ${MPI_SOURCES})
  include_directories(${LAMMPS_SOURCE_DIR}/STUBS)
endif()

include(StyleHeaderUtils)
RegisterStyles(${LAMMPS_SOURCE_DIR})

# packages which include entire content when enabled

foreach(PKG ${PACKAGES})
  if(ENABLE_${PKG})
    set(${PKG}_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/${PKG})

    # detects styles in package and adds them to global list
    RegisterStyles(${${PKG}_SOURCES_DIR})

    file(GLOB ${PKG}_SOURCES ${${PKG}_SOURCES_DIR}/*.cpp)
    list(APPEND LIB_SOURCES ${${PKG}_SOURCES})
    include_directories(${${PKG}_SOURCES_DIR})
  endif()
endforeach()

# packages which selectively include variants based on enabled styles
# e.g. accelerator packages

if(ENABLE_USER_OMP)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        message(FATAL_ERROR "USER-OMP requires a compiler with OpenMP support")
    endif()

    set(USER_OMP_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-OMP)
    set(USER_OMP_SOURCES ${USER_OMP_SOURCES_DIR}/thr_data.cpp
                         ${USER_OMP_SOURCES_DIR}/thr_omp.cpp
                         ${USER_OMP_SOURCES_DIR}/fix_nh_omp.cpp
                         ${USER_OMP_SOURCES_DIR}/fix_nh_sphere_omp.cpp)
    set_property(GLOBAL PROPERTY "OMP_SOURCES" "${USER_OMP_SOURCES}")

    # detects styles which have USER-OMP version
    RegisterStylesExt(${USER_OMP_SOURCES_DIR} omp OMP_SOURCES)

    get_property(USER_OMP_SOURCES GLOBAL PROPERTY OMP_SOURCES)

    list(APPEND LIB_SOURCES ${USER_OMP_SOURCES})
    include_directories(${USER_OMP_SOURCES_DIR})
endif()

if(ENABLE_KOKKOS)
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else()
        message(FATAL_ERROR "USER-OMP requires a compiler with OpenMP support")
    endif()

    set(KOKKOS_PKG_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/KOKKOS)
    set(KOKKOS_PKG_SOURCES ${KOKKOS_PKG_SOURCES_DIR}/kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/atom_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/atom_vec_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/comm_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/comm_tiled_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/neighbor_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/neigh_list_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/neigh_bond_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/fix_nh_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/domain_kokkos.cpp
                           ${KOKKOS_PKG_SOURCES_DIR}/modify_kokkos.cpp)
    set_property(GLOBAL PROPERTY "KOKKOS_PKG_SOURCES" "${KOKKOS_PKG_SOURCES}")

    # detects styles which have KOKKOS version
    RegisterStylesExt(${KOKKOS_PKG_SOURCES_DIR} kokkos KOKKOS_PKG_SOURCES)

    get_property(KOKKOS_PKG_SOURCES GLOBAL PROPERTY KOKKOS_PKG_SOURCES)

    list(APPEND LIB_SOURCES ${KOKKOS_PKG_SOURCES})
    include_directories(${KOKKOS_PKG_SOURCES_DIR})
endif()


if(ENABLE_REAX OR ENABLE_MEAM)
  enable_language(Fortran)
endif()

if(ENABLE_REAX)
  file(GLOB REAX_SOURCES ${LAMMPS_LIB_SOURCE_DIR}/reax/*.F)
  list(APPEND LIB_SOURCES ${REAX_SOURCES})
  include_directories(${LAMMPS_LIB_SOURCE_DIR}/reax)
endif()

if(ENABLE_MEAM)
  file(GLOB MEAM_SOURCES ${LAMMPS_LIB_SOURCE_DIR}/meam/*.F ${LAMMPS_LIB_SOURCE_DIR}/meam/*.c)
  list(APPEND LIB_SOURCES ${MEAM_SOURCES})
  include_directories(${LAMMPS_LIB_SOURCE_DIR}/meam)
endif()


######################################################
# Generate style headers based on global list of
# styles registered during package selection           
######################################################
set(LAMMPS_STYLE_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR}/styles)

GenerateStyleHeaders(${LAMMPS_STYLE_HEADERS_DIR})

include_directories(${LAMMPS_SOURCE_DIR})
include_directories(${LAMMPS_STYLE_HEADERS_DIR})


add_library(lammps ${LIB_SOURCES})
target_link_libraries(lammps ${LAMMPS_LINK_LIBS})
set_target_properties(lammps PROPERTIES SOVERSION ${SOVERSION})

if(ENABLE_KOKKOS)
    target_link_libraries(lammps kokkos)
endif()

if(INSTALL_LIB)
  install(TARGETS lammps LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(FILES ${LAMMPS_SOURCE_DIR}/lammps.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
elseif(NOT BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Shared library has to install, use -DBUILD_SHARED_LIBS=OFF to install lammps with a a library")
endif()

add_executable(lmp ${LMP_SOURCES})
target_link_libraries(lmp lammps)

install(TARGETS lmp DESTINATION ${CMAKE_INSTALL_BINDIR})

foreach(PKG ${PACKAGES} ${ACCEL_PACKAGES})
  if(ENABLE_${PKG})
    message(STATUS "Building package: ${PKG}")
  endif()
endforeach()
