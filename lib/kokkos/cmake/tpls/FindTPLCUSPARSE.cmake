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
#@HEADER

# Check for CUDA support

IF (NOT TPL_ENABLE_CUDA)
  MESSAGE(FATAL_ERROR "\nCUSPARSE requires CUDA")
ELSE()
  GLOBAL_SET(TPL_CUSPARSE_LIBRARY_DIRS)
  GLOBAL_SET(TPL_CUSPARSE_INCLUDE_DIRS ${TPL_CUDA_INCLUDE_DIRS})
  GLOBAL_SET(TPL_CUSPARSE_LIBRARIES    ${CUDA_cusparse_LIBRARY})
ENDIF()

