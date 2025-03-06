# -*- CMake -*- configuration for generating LAMMPS script language wrappers with SWIG

# set minimum CMake version required and switch to new policies for SWIG
cmake_minimum_required(VERSION 3.16)
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 NEW)
endif()
enable_language(C)

# some of the find_package() scripts trigger developer warnings
# even though they are bundled with CMake; Suppress them
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE)

if(NOT BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Option BUILD_SHARED_LIBS must be enabled to use SWIG wrappers")
endif()

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

option(BUILD_SWIG_JAVA "Build Java JNI wrapper with SWIG" OFF)
option(BUILD_SWIG_LUA "Build Lua wrapper with SWIG" OFF)
option(BUILD_SWIG_PERL5 "Build Perl5 wrapper with SWIG" OFF)
option(BUILD_SWIG_PYTHON "Build Python wrapper with SWIG" OFF)
option(BUILD_SWIG_RUBY "Build Ruby wrapper with SWIG" OFF)
option(BUILD_SWIG_TCL "Build Tcl wrapper with SWIG" OFF)

set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
        pllammps.pm __pycache__ log.lammps
        example.class example.java example.lua example.pl example.py example.rb example.tcl example.tcllmp
        javalammps.class javalammpsJNI.class javalammps.java javalammpsJNI.java
        _LMP_DATATYPE_CONST.class _LMP_STYLE_CONST.class _LMP_TYPE_CONST.class
        _LMP_DATATYPE_CONST.java _LMP_STYLE_CONST.java _LMP_TYPE_CONST.java
        SWIGTYPE_p_double.class SWIGTYPE_p_int.class SWIGTYPE_p_p_char.class SWIGTYPE_p_p_double.class
        SWIGTYPE_p_p_int.class SWIGTYPE_p_p_void.class SWIGTYPE_p_void.class
        SWIGTYPE_p_double.java SWIGTYPE_p_int.java SWIGTYPE_p_p_char.java SWIGTYPE_p_p_double.java
        SWIGTYPE_p_p_int.java SWIGTYPE_p_p_void.java SWIGTYPE_p_void.java)

if(BUILD_SWIG_JAVA)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME javalammps)
  swig_add_library(javalammps TYPE MODULE LANGUAGE java SOURCES lammps.i)
  find_package(JNI REQUIRED)
  find_package(Java COMPONENTS Development REQUIRED)
  target_include_directories(javalammps PRIVATE ${JNI_INCLUDE_DIRS})
  swig_link_libraries(javalammps PRIVATE lammps ${JNI_LIBRARIES})
  configure_file(run_java_example.sh.in run_java_example.sh @ONLY)
endif()

if(BUILD_SWIG_LUA)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME lualammps)
  swig_add_library(lualammps TYPE MODULE LANGUAGE lua SOURCES lammps.i)
  find_package(Lua REQUIRED)
  target_include_directories(lualammps PRIVATE ${LUA_INCLUDE_PATH})
  swig_link_libraries(lualammps PRIVATE lammps ${LUA_LIBRARY})
  configure_file(run_lua_example.sh.in run_lua_example.sh @ONLY)
endif()

if(BUILD_SWIG_PERL5)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME pllammps)
  swig_add_library(pllammps TYPE MODULE LANGUAGE perl5 SOURCES lammps.i)
  find_package(Perl)
  find_package(PerlLibs REQUIRED 5.0)
  target_compile_definitions(pllammps PRIVATE _LARGEFILE64_SOURCE)
  target_include_directories(pllammps PRIVATE ${PERL_INCLUDE_PATH})
  swig_link_libraries(pllammps PRIVATE lammps ${PERL_LIBRARY})
  configure_file(run_perl_example.sh.in run_perl_example.sh @ONLY)
endif()

if(BUILD_SWIG_PYTHON)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME pylammps)
  swig_add_library(pylammps TYPE MODULE LANGUAGE python SOURCES lammps.i)
  find_package(Python3 REQUIRED COMPONENTS Development Interpreter)
  swig_link_libraries(pylammps PRIVATE lammps Python3::Python)
  configure_file(run_python_example.sh.in run_python_example.sh @ONLY)
endif()

if(BUILD_SWIG_RUBY)
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME rubylammps)
  swig_add_library(rubylammps TYPE MODULE LANGUAGE ruby SOURCES lammps.i)
  find_package(Ruby REQUIRED)
  target_include_directories(rubylammps PRIVATE ${RUBY_INCLUDE_DIRS})
  swig_link_libraries(rubylammps PRIVATE lammps ${RUBY_LIBRARY})
  configure_file(run_ruby_example.sh.in run_ruby_example.sh @ONLY)
endif()

if(BUILD_SWIG_TCL)
  # build loadable Tcl module
  set_property(SOURCE lammps.i PROPERTY SWIG_MODULE_NAME tcllammps)
  swig_add_library(tcllammps TYPE MODULE LANGUAGE tcl SOURCES lammps.i)
  find_package(TCL)
  if(NOT TCL_FOUND)
    message(FATAL_ERROR "Tcl development headers and libraries are required")
  endif()
  find_package(TclStub)
  if(TCL_STUB_LIBRARY)
    target_compile_definitions(tcllammps PRIVATE USE_TCL_STUBS)
    target_link_libraries(tcllammps PRIVATE ${TCL_STUB_LIBRARY})
  endif()
  target_include_directories(tcllammps PRIVATE ${TCL_INCLUDE_PATH})
  swig_link_libraries(tcllammps PRIVATE lammps ${TCL_LIBRARY})
  # build extended Tcl shell binary
  set_property(SOURCE lammps_shell.i PROPERTY SWIG_MODULE_NAME tcllmpsh)
  swig_add_library(libtcllmpsh TYPE STATIC LANGUAGE tcl SOURCES lammps_shell.i)
  add_executable(tcllmpsh tcldummy.c)
  target_include_directories(libtcllmpsh PRIVATE ${TCL_INCLUDE_PATH})
  target_link_libraries(tcllmpsh PRIVATE libtcllmpsh lammps ${TCL_LIBRARY})
  configure_file(run_tcl_example.sh.in run_tcl_example.sh @ONLY)
endif()

