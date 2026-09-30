# SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception

#-----------------------------------------------------------------------------
#  Hardware locality detection and control library.
#
#  Acquisition information:
#    Date checked:  November 2011
#    Checked by:    H. Carter Edwards <hcedwar AT sandia.gov>
#    Source:        http://www.open-mpi.org/projects/hwloc/
#    Version:       1.3
#

kokkos_tpl_find_include_dirs_and_libraries(HWLOC REQUIRED_HEADERS hwloc.h REQUIRED_LIBS_NAMES "hwloc")
