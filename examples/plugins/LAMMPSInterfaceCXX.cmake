# Cmake script code to define the LAMMPS C++ interface
# settings required for building LAMMPS plugins

################################################################################
# helper function
function(validate_option name values)
  string(TOLOWER ${${name}} needle_lower)
  string(TOUPPER ${${name}} needle_upper)
  list(FIND ${values} ${needle_lower} IDX_LOWER)
  list(FIND ${values} ${needle_upper} IDX_UPPER)
  if(${IDX_LOWER} LESS 0 AND ${IDX_UPPER} LESS 0)
    list_to_bulletpoints(POSSIBLE_VALUE_LIST ${${values}})
    message(FATAL_ERROR "\n########################################################################\n"
      "Invalid value '${${name}}' for option ${name}\n"
      "\n"
      "Possible values are:\n"
      "${POSSIBLE_VALUE_LIST}"
      "########################################################################")
  endif()
endfunction(validate_option)

#################################################################################
# LAMMPS C++ interface. We only need the header related parts.
add_library(lammps INTERFACE)
target_include_directories(lammps INTERFACE ${LAMMPS_HEADER_DIR})

################################################################################
# MPI configuration
if(NOT CMAKE_CROSSCOMPILING)
  set(MPI_CXX_SKIP_MPICXX TRUE)
  find_package(MPI QUIET)
  option(BUILD_MPI "Build MPI version" ${MPI_FOUND})
else()
  option(BUILD_MPI "Build MPI version" OFF)
endif()

if(BUILD_MPI)
  find_package(MPI REQUIRED)
  option(LAMMPS_LONGLONG_TO_LONG "Workaround if your system or MPI version does not recognize 'long long' data types" OFF)
  if(LAMMPS_LONGLONG_TO_LONG)
    target_compile_definitions(lammps INTERFACE -DLAMMPS_LONGLONG_TO_LONG)
  endif()
  target_link_libraries(lammps INTERFACE MPI::MPI_CXX)
else()
  target_include_directories(lammps INTERFACE "${LAMMPS_SOURCE_DIR}/STUBS")
endif()

set(LAMMPS_SIZES "smallbig" CACHE STRING "LAMMPS integer sizes (smallsmall: all 32-bit, smallbig: 64-bit #atoms #timesteps, bigbig: also 64-bit imageint, 64-bit atom ids)")
set(LAMMPS_SIZES_VALUES smallbig bigbig smallsmall)
set_property(CACHE LAMMPS_SIZES PROPERTY STRINGS ${LAMMPS_SIZES_VALUES})
validate_option(LAMMPS_SIZES LAMMPS_SIZES_VALUES)
string(TOUPPER ${LAMMPS_SIZES} LAMMPS_SIZES)
target_compile_definitions(lammps INTERFACE -DLAMMPS_${LAMMPS_SIZES})

################################################################################
# detect if we may enable OpenMP support by default
set(BUILD_OMP_DEFAULT OFF)
find_package(OpenMP QUIET)
if(OpenMP_FOUND)
  check_include_file_cxx(omp.h HAVE_OMP_H_INCLUDE)
  if(HAVE_OMP_H_INCLUDE)
    set(BUILD_OMP_DEFAULT ON)
  endif()
endif()

option(BUILD_OMP "Build with OpenMP support" ${BUILD_OMP_DEFAULT})

if(BUILD_OMP)
  find_package(OpenMP REQUIRED)
  check_include_file_cxx(omp.h HAVE_OMP_H_INCLUDE)
  if(NOT HAVE_OMP_H_INCLUDE)
    message(FATAL_ERROR "Cannot find the 'omp.h' header file required for full OpenMP support")
  endif()

  if (((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.0)) OR
      (CMAKE_CXX_COMPILER_ID STREQUAL "PGI") OR
      ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)) OR
      ((CMAKE_CXX_COMPILER_ID STREQUAL "Intel") AND (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.0)))
    # GCC 9.x and later plus Clang 10.x and later implement strict OpenMP 4.0 semantics for consts.
    # Intel 18.0 was tested to support both, so we switch to OpenMP 4+ from 19.x onward to be safe.
    target_compile_definitions(lammps INTERFACE -DLAMMPS_OMP_COMPAT=4)
  else()
    target_compile_definitions(lammps INTERFACE -DLAMMPS_OMP_COMPAT=3)
  endif()
  target_link_libraries(lammps INTERFACE OpenMP::OpenMP_CXX)
endif()
