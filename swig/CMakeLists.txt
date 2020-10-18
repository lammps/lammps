# CMake configuration for generating script language interfaces with swig

# set minimum CMake version required and switch to new policies
cmake_minimum_required(VERSION 3.14)
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()

if(NOT BUILD_SHARED_LIBS)
  error(FATAL_ERROR "Option BUILD_SHARED_LIBS must be enabled to use SWIG wrappers")
endif()

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

option(BUILD_SWIG_JAVA "Build Java JNI wrapper with SWIG" OFF)
option(BUILD_SWIG_LUA "Build Lua wrapper with SWIG" OFF)
option(BUILD_SWIG_PERL5 "Build Perl5 wrapper with SWIG" OFF)
option(BUILD_SWIG_PYTHON "Build Python wrapper with SWIG" OFF)
option(BUILD_SWIG_TCL "Build Tcl wrapper with SWIG" OFF)

if(BUILD_SWIG_JAVA)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME javalammps)
  swig_add_library(javalammps TYPE MODULE LANGUAGE java SOURCES lammps.i)
  find_package(JNI REQUIRED)
  target_include_directories(javalammps PRIVATE ${JNI_INCLUDE_DIRS})
  swig_link_libraries(javalammps PRIVATE lammps ${JNI_LIBRARIES})
endif()

if(BUILD_SWIG_PYTHON)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME pylammps)
  swig_add_library(pylammps TYPE MODULE LANGUAGE python SOURCES lammps.i)
  find_package(Python REQUIRED COMPONENTS Development)
  swig_link_libraries(pylammps PRIVATE lammps Python::Python)
endif()

if(BUILD_SWIG_PERL5)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME pllammps)
  swig_add_library(pllammps TYPE MODULE LANGUAGE perl5 SOURCES lammps.i)
  find_package(PerlLibs REQUIRED 5.0)
  target_compile_definitions(pllammps PRIVATE _LARGEFILE64_SOURCE)
  target_include_directories(pllammps PRIVATE ${PERL_INCLUDE_PATH})
  swig_link_libraries(pllammps PRIVATE lammps ${PERL_LIBRARY})
endif()

if(BUILD_SWIG_LUA)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME lualammps)
  swig_add_library(lualammps TYPE MODULE LANGUAGE lua SOURCES lammps.i)
  find_package(Lua REQUIRED)
  target_include_directories(lualammps PRIVATE ${LUA_INCLUDE_PATH})
  swig_link_libraries(lualammps PRIVATE lammps ${LUA_LIBRARY})
endif()

if(BUILD_SWIG_TCL)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME tcllammps)
  swig_add_library(tcllammps TYPE MODULE LANGUAGE tcl SOURCES lammps.i)
  find_package(TCL REQUIRED)
  target_include_directories(tcllammps PRIVATE ${TCL_INCLUDE_PATH})
  swig_link_libraries(tcllammps PRIVATE lammps ${TCL_LIBRARY})
endif()

