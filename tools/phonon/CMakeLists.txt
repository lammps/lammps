
# Support Linux from Ubuntu 20.04LTS onward, CentOS 7.x (with EPEL),
# macOS, MSVC 2019 (=Version 16)
cmake_minimum_required(VERSION 3.10)

# set timestamp of downloaded files to that of archive
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# set up project
set(PHANA_MINOR_VERSION 48)
project(phonon VERSION ${PHANA_MINOR_VERSION}
  DESCRIPTION "Fix phonon post-processor"
  LANGUAGES CXX C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# hacks for MSVC to prevent lots of pointless warnings about "unsafe" functions,
# padding and Spectre mitigation
if(MSVC)
  add_compile_options(/wd4244)
  add_compile_options(/wd4267)
  add_compile_options(/wd4711)
  add_compile_options(/wd4820)
  add_compile_options(/wd5045)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

set(CMAKE_PROJECT_VERSION ${PHANA_MINOR_VERSION})
configure_file(version.h.in version.h @ONLY)
add_executable(phana
  main.cpp
  disp.cpp
  dynmat.cpp
  green.cpp
  input.cpp
  interpolate.cpp
  kpath.cpp
  memory.cpp
  phonon.cpp
  phonopy.cpp
  qnodes.cpp
  timer.cpp
)
target_include_directories(phana PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

if(NOT LAMMPS_DIR)
  set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules)
  set(LAMMPS_THIRDPARTY_URL "https://download.lammps.org/thirdparty")
endif()
find_package(FFTW3)
if(FFTW3_FOUND)
  target_compile_definitions(phana PRIVATE FFTW3)
  target_link_libraries(phana PRIVATE FFTW3::FFTW3)
endif()

# build bundeled libraries
add_subdirectory(tricubic)

# standalone build must build our own version of linalg
if(NOT LAMMPS_DIR)
  if(NOT USE_INTERNAL_LINALG)
    find_package(LAPACK)
    find_package(BLAS)
  endif()
  if(NOT LAPACK_FOUND OR NOT BLAS_FOUND OR USE_INTERNAL_LINALG)
    file(GLOB LINALG_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/linalg/[^.]*.cpp)
    add_library(linalg STATIC ${LINALG_SOURCES})
    set(BLAS_LIBRARIES "$<TARGET_FILE:linalg>")
    set(LAPACK_LIBRARIES "$<TARGET_FILE:linalg>")
  else()
    list(APPEND LAPACK_LIBRARIES ${BLAS_LIBRARIES})
  endif()
endif()

option(USE_SPGLIB "Download and use spglib for phonon DOS and other optional properties" ON)
if(USE_SPGLIB)
  set(SPGLIB_URL "https://github.com/spglib/spglib/archive/refs/tags/v1.11.2.1.tar.gz" CACHE STRING "URL for spglib v1.x tarball")
  set(SPGLIB_MD5 "3089782bc85b5034dd4765a18ee70bc7" CACHE STRING "MD5 checksum for spglib tarball")
  mark_as_advanced(SPGLIB_URL)
  mark_as_advanced(SPGLIB_MD5)
  include(LAMMPSUtils)
  GetFallbackURL(SPGLIB_URL SPGLIB_FALLBACK)

  include(ExternalProject)
  ExternalProject_Add(spglib_build
    URL ${SPGLIB_URL} ${SPGLIB_FALLBACK}
    URL_MD5 ${SPGLIB_MD5}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_MAKE_PROGRAM=${CMAKE_MAKE_PROGRAM}
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    UPDATE_COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.spglib ${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext/src/spglib_build/CMakeLists.txt
    INSTALL_COMMAND ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext/src/spglib_build-build --target install
    BUILD_BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext/lib/${CMAKE_STATIC_LIBRARY_PREFIX}symspg${CMAKE_STATIC_LIBRARY_SUFFIX}"
  )

  # workaround for older CMake versions
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext/lib)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext/include)

  add_library(SPGLIB::SYMSPG UNKNOWN IMPORTED)
  add_dependencies(SPGLIB::SYMSPG spglib_build)
  set_target_properties(SPGLIB::SYMSPG PROPERTIES
    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext/lib/${CMAKE_STATIC_LIBRARY_PREFIX}symspg${CMAKE_STATIC_LIBRARY_SUFFIX}"
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR}/spglib_build_ext/include
  )
    
  target_compile_definitions(phana PRIVATE UseSPG)
  target_link_libraries(phana PRIVATE SPGLIB::SYMSPG)
endif()

target_link_libraries(phana PRIVATE tricubic ${LAPACK_LIBRARIES})
install(TARGETS phana EXPORT LAMMPS_Targets DESTINATION ${CMAKE_INSTALL_BINDIR})
