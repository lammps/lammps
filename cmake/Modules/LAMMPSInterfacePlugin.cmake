# CMake script code to define LAMMPS settings required for building LAMMPS plugins

# enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. You must create and use a build directory. "
    "Please remove CMakeCache.txt and CMakeFiles first.")
endif()

set(LAMMPS_THIRDPARTY_URL "https://download.lammps.org/thirdparty"
  CACHE STRING "URL for thirdparty package downloads")

# global LAMMPS/plugin build settings
set(LAMMPS_SOURCE_DIR ""  CACHE PATH "Location of LAMMPS sources folder")
if(NOT LAMMPS_SOURCE_DIR)
  message(FATAL_ERROR "Must set LAMMPS_SOURCE_DIR")
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

# C++11 is required
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Need -restrict with Intel compilers
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -restrict")
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)

#######
# helper functions from LAMMPSUtils.cmake
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

# helper function for getting the most recently modified file or folder from a glob pattern
function(get_newest_file path variable)
  file(GLOB _dirs ${CONFIGURE_DEPENDS} ${path})
  set(_besttime 2000-01-01T00:00:00)
  set(_bestfile "<unknown>")
  foreach(_dir ${_dirs})
    file(TIMESTAMP ${_dir} _newtime)
    if(_newtime IS_NEWER_THAN _besttime)
      set(_bestfile ${_dir})
      set(_besttime ${_newtime})
    endif()
  endforeach()
  if(_bestfile STREQUAL "<unknown>")
    message(FATAL_ERROR "Could not find valid path at: ${path}")
  endif()
  set(${variable} ${_bestfile} PARENT_SCOPE)
endfunction()

# get LAMMPS version date
function(get_lammps_version version_header variable)
    file(STRINGS ${version_header} line REGEX LAMMPS_VERSION)
    string(REGEX REPLACE "#define LAMMPS_VERSION \"([0-9]+) ([A-Za-z]+) ([0-9]+)\"" "\\1\\2\\3" date "${line}")
    set(${variable} "${date}" PARENT_SCOPE)
endfunction()

#################################################################################
# LAMMPS C++ interface. We only need the header related parts except on windows.
add_library(lammps INTERFACE)
target_include_directories(lammps INTERFACE ${LAMMPS_SOURCE_DIR})
if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND CMAKE_CROSSCOMPILING)
  target_link_libraries(lammps INTERFACE ${CMAKE_BINARY_DIR}/../liblammps.dll.a)
endif()

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
  # do not include the (obsolete) MPI C++ bindings which makes
  # for leaner object files and avoids namespace conflicts
  set(MPI_CXX_SKIP_MPICXX TRUE)
  # We use a non-standard procedure to cross-compile with MPI on Windows
  if((CMAKE_SYSTEM_NAME STREQUAL "Windows") AND CMAKE_CROSSCOMPILING)
    # Download and configure MinGW compatible MPICH development files for Windows
    option(USE_MSMPI "Use Microsoft's MS-MPI SDK instead of MPICH2-1.4.1" OFF)
    if(USE_MSMPI)
      message(STATUS "Downloading and configuring MS-MPI 10.1 for Windows cross-compilation")
      set(MPICH2_WIN64_DEVEL_URL "${LAMMPS_THIRDPARTY_URL}/msmpi-win64-devel.tar.gz" CACHE STRING "URL for MS-MPI (win64) tarball")
      set(MPICH2_WIN64_DEVEL_MD5 "86314daf1bffb809f1fcbefb8a547490" CACHE STRING "MD5 checksum of MS-MPI (win64) tarball")
      mark_as_advanced(MPICH2_WIN64_DEVEL_URL)
      mark_as_advanced(MPICH2_WIN64_DEVEL_MD5)

      include(ExternalProject)
      if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        ExternalProject_Add(mpi4win_build
          URL     ${MPICH2_WIN64_DEVEL_URL}
          URL_MD5 ${MPICH2_WIN64_DEVEL_MD5}
          CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
          BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libmsmpi.a)
      else()
        message(FATAL_ERROR "Only x86 64-bit builds are supported with MS-MPI")
      endif()

      ExternalProject_get_property(mpi4win_build SOURCE_DIR)
      file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
      add_library(MPI::MPI_CXX UNKNOWN IMPORTED)
      set_target_properties(MPI::MPI_CXX PROPERTIES
        IMPORTED_LOCATION "${SOURCE_DIR}/lib/libmsmpi.a"
        INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include"
        INTERFACE_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX")
      add_dependencies(MPI::MPI_CXX mpi4win_build)

      # set variables for status reporting at the end of CMake run
      set(MPI_CXX_INCLUDE_PATH "${SOURCE_DIR}/include")
      set(MPI_CXX_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX")
      set(MPI_CXX_LIBRARIES "${SOURCE_DIR}/lib/libmsmpi.a")
    else()
      # Download and configure custom MPICH files for Windows
      message(STATUS "Downloading and configuring MPICH-1.4.1 for Windows")
      set(MPICH2_WIN64_DEVEL_URL "${LAMMPS_THIRDPARTY_URL}/mpich2-win64-devel.tar.gz" CACHE STRING "URL for MPICH2 (win64) tarball")
      set(MPICH2_WIN64_DEVEL_MD5 "4939fdb59d13182fd5dd65211e469f14" CACHE STRING "MD5 checksum of MPICH2 (win64) tarball")
      mark_as_advanced(MPICH2_WIN64_DEVEL_URL)
      mark_as_advanced(MPICH2_WIN64_DEVEL_MD5)

      include(ExternalProject)
      if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
        ExternalProject_Add(mpi4win_build
          URL     ${MPICH2_WIN64_DEVEL_URL}
          URL_MD5 ${MPICH2_WIN64_DEVEL_MD5}
          CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
          BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libmpi.a)
      else()
        ExternalProject_Add(mpi4win_build
          URL     ${MPICH2_WIN32_DEVEL_URL}
          URL_MD5 ${MPICH2_WIN32_DEVEL_MD5}
          CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND ""
          BUILD_BYPRODUCTS <SOURCE_DIR>/lib/libmpi.a)
      endif()

      ExternalProject_get_property(mpi4win_build SOURCE_DIR)
      file(MAKE_DIRECTORY "${SOURCE_DIR}/include")
      add_library(MPI::MPI_CXX UNKNOWN IMPORTED)
      set_target_properties(MPI::MPI_CXX PROPERTIES
        IMPORTED_LOCATION "${SOURCE_DIR}/lib/libmpi.a"
        INTERFACE_INCLUDE_DIRECTORIES "${SOURCE_DIR}/include"
        INTERFACE_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX")
      add_dependencies(MPI::MPI_CXX mpi4win_build)

      # set variables for status reporting at the end of CMake run
      set(MPI_CXX_INCLUDE_PATH "${SOURCE_DIR}/include")
      set(MPI_CXX_COMPILE_DEFINITIONS "MPICH_SKIP_MPICXX")
      set(MPI_CXX_LIBRARIES "${SOURCE_DIR}/lib/libmpi.a")
    endif()
  else()
    find_package(MPI REQUIRED)
    option(LAMMPS_LONGLONG_TO_LONG "Workaround if your system or MPI version does not recognize 'long long' data types" OFF)
    if(LAMMPS_LONGLONG_TO_LONG)
      target_compile_definitions(lammps INTERFACE -DLAMMPS_LONGLONG_TO_LONG)
    endif()
  endif()
  target_link_libraries(lammps INTERFACE MPI::MPI_CXX)
else()
  add_library(mpi_stubs INTERFACE)
  target_include_directories(mpi_stubs INTERFACE $<BUILD_INTERFACE:${LAMMPS_SOURCE_DIR}/STUBS>)
  target_link_libraries(lammps INTERFACE mpi_stubs)
endif()

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

################
# integer size selection
set(LAMMPS_SIZES "smallbig" CACHE STRING "LAMMPS integer sizes (smallsmall: all 32-bit, smallbig: 64-bit #atoms #timesteps, bigbig: also 64-bit imageint, 64-bit atom ids)")
set(LAMMPS_SIZES_VALUES smallbig bigbig smallsmall)
set_property(CACHE LAMMPS_SIZES PROPERTY STRINGS ${LAMMPS_SIZES_VALUES})
validate_option(LAMMPS_SIZES LAMMPS_SIZES_VALUES)
string(TOUPPER ${LAMMPS_SIZES} LAMMPS_SIZES)
target_compile_definitions(lammps INTERFACE -DLAMMPS_${LAMMPS_SIZES})
