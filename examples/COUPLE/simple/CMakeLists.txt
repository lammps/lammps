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
target_link_libraries(simpleCC LAMMPS::lammps)

add_executable(simpleC simple.c)
target_link_libraries(simpleC LAMMPS::lammps)
