# AutotoolsTouch.cmake -- restore the timestamp order that GNU autotools generated files require
#
# Usage: ${CMAKE_COMMAND} -D SOURCE_DIR=<path to extracted sources> -P AutotoolsTouch.cmake
#
# With policy CMP0135 set to NEW, files extracted from a downloaded archive get the time of
# their extraction as timestamp (in the order they are stored in the archive) instead of the
# timestamp recorded in the archive.  This is desirable, since it guarantees that all
# dependent objects are rebuilt when the archive is updated.  But for source trees prepared
# with GNU autotools (autoconf, autoheader, automake), the generated Makefiles contain rules
# to re-run those tools whenever their inputs (configure.ac, Makefile.am, m4/*.m4) are newer
# than their outputs (aclocal.m4, configure, config.h.in, Makefile.in).  Depending on the
# order of the files in the archive this can happen after extraction.  It then fails on
# systems without (the matching version of) the autotools installed and is not wanted anyway.
#
# This script is meant to be run as the PATCH_COMMAND of an ExternalProject_Add() call.
# It touches the generated files in dependency order so that afterward
#   configure.ac, Makefile.am, m4/*.m4  <  aclocal.m4  <  configure, config.h.in, Makefile.in
# holds.  All lookups are recursive so that nested configure scripts (e.g. in ScaFaCoS) are
# handled as well.  Source trees without automake (e.g. PLUMED) only get their configure
# script touched, which is harmless.

if(NOT SOURCE_DIR)
  message(FATAL_ERROR "AutotoolsTouch.cmake requires -D SOURCE_DIR=<path to extracted sources>")
endif()
if(NOT IS_DIRECTORY "${SOURCE_DIR}")
  message(FATAL_ERROR "AutotoolsTouch.cmake: SOURCE_DIR ${SOURCE_DIR} is not a directory")
endif()

message(STATUS "Restoring autotools timestamp order in ${SOURCE_DIR}")

# step 1: output of aclocal, which is an input to the other tools
file(GLOB_RECURSE _aclocal_outputs LIST_DIRECTORIES false "${SOURCE_DIR}/aclocal.m4")
if(_aclocal_outputs)
  file(TOUCH ${_aclocal_outputs})
endif()

# step 2: outputs of autoconf, autoheader, and automake
file(GLOB_RECURSE _autotools_outputs LIST_DIRECTORIES false
     "${SOURCE_DIR}/configure" "${SOURCE_DIR}/Makefile.in" "${SOURCE_DIR}/*.h.in")
if(_autotools_outputs)
  file(TOUCH ${_autotools_outputs})
endif()
