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

#include(${TRIBITS_DEPS_DIR}/CUDA.cmake)

#IF (TPL_ENABLE_CUDA)
#  GLOBAL_SET(TPL_CUSPARSE_LIBRARY_DIRS)
#  GLOBAL_SET(TPL_CUSPARSE_INCLUDE_DIRS ${TPL_CUDA_INCLUDE_DIRS})
#  GLOBAL_SET(TPL_CUSPARSE_LIBRARIES    ${CUDA_cusparse_LIBRARY})
#  KOKKOS_CREATE_IMPORTED_TPL_LIBRARY(CUSPARSE)
#ENDIF()

