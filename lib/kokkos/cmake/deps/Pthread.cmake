#@HEADER
# ************************************************************************
#
#                        Kokkos v. 4.0
#       Copyright (2022) National Technology & Engineering
#               Solutions of Sandia, LLC (NTESS).
#
# Under the terms of Contract DE-NA0003525 with NTESS,
# the U.S. Government retains certain rights in this software.
#
# Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
#
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
#
# ************************************************************************
# @HEADER


SET(USE_THREADS FALSE)

IF(NOT TPL_Pthread_INCLUDE_DIRS AND NOT TPL_Pthread_LIBRARY_DIRS AND NOT TPL_Pthread_LIBRARIES)
  # Use CMake's Thread finder since it is a bit smarter in determining
  # whether pthreads is already built into the compiler and doesn't need
  # a library to link.
  FIND_PACKAGE(Threads)
  #If Threads found a copy of pthreads make sure it is one of the cases the tribits
  #tpl system cannot handle.
  IF(Threads_FOUND AND CMAKE_USE_PTHREADS_INIT)
    IF(CMAKE_THREAD_LIBS_INIT STREQUAL "" OR CMAKE_THREAD_LIBS_INIT STREQUAL "-pthread")
      SET(USE_THREADS TRUE)
    ENDIF()
  ENDIF()
ENDIF()

IF(USE_THREADS)
  SET(TPL_Pthread_INCLUDE_DIRS "")
  SET(TPL_Pthread_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
  SET(TPL_Pthread_LIBRARY_DIRS "")
  KOKKOS_CREATE_IMPORTED_TPL_LIBRARY(Pthread)
ELSE()
  KOKKOS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Pthread
    REQUIRED_HEADERS pthread.h
    REQUIRED_LIBS_NAMES pthread
      )
ENDIF()
