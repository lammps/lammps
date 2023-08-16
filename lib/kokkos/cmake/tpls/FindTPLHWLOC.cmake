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


#-----------------------------------------------------------------------------
#  Hardware locality detection and control library.
#
#  Acquisition information:
#    Date checked:  November 2011
#    Checked by:    H. Carter Edwards <hcedwar AT sandia.gov>
#    Source:        http://www.open-mpi.org/projects/hwloc/
#    Version:       1.3
#

KOKKOS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( HWLOC
  REQUIRED_HEADERS hwloc.h
  REQUIRED_LIBS_NAMES "hwloc"
  )

