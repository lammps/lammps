# -*- CMake -*- file for example programs that use the LAMMPS library
# where the library is loaded dynamically at runtime.
##########################################

cmake_minimum_required(VERSION 3.16)

# enforce out-of-source build
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. You must create and use a build directory. "
    "Please remove CMakeCache.txt and CMakeFiles first.")
endif()

project(liblammpsplugin VERSION 1.0 LANGUAGES C)

# by default, install into $HOME/.local (not /usr/local),
# so that no root access (and sudo) is needed
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local" CACHE PATH "Default install path" FORCE)
endif()

# ugly hacks for MSVC which by default always reports an old C++ standard in the __cplusplus macro
# and prints lots of pointless warnings about "unsafe" functions
if(MSVC)
  if(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(/EHsc)
  endif()
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
endif()

find_package(MPI REQUIRED)
# do not include the (obsolete) MPI C++ bindings which makes
# for leaner object files and avoids namespace conflicts
set(MPI_CXX_SKIP_MPICXX TRUE)

##########################

add_executable(simple-plugin simple.c liblammpsplugin.c)
target_link_libraries(simple-plugin PRIVATE MPI::MPI_C)
target_compile_definitions(simple-plugin PRIVATE LAMMPS_LIB_MPI)

# link with -ldl or equivalent for plugin loading; except on Windows
if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
   target_link_libraries(simple-plugin PRIVATE ${CMAKE_DL_LIBS})
endif()
