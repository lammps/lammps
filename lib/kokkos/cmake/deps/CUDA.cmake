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

# Check for CUDA support

set(_CUDA_FAILURE OFF)

# Have CMake find CUDA
if(NOT _CUDA_FAILURE)
  find_package(CUDA 3.2)
  if(NOT CUDA_FOUND)
    set(_CUDA_FAILURE ON)
  endif()
endif()

if(NOT _CUDA_FAILURE)
  # if we haven't met failure
  macro(PACKAGE_ADD_CUDA_LIBRARY cuda_target)
    tribits_add_library(${cuda_target} ${ARGN} CUDALIBRARY)
  endmacro()
  global_set(TPL_CUDA_LIBRARY_DIRS)
  global_set(TPL_CUDA_INCLUDE_DIRS ${CUDA_TOOLKIT_INCLUDE})
  global_set(TPL_CUDA_LIBRARIES ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cufft_LIBRARY})
else()
  set(TPL_ENABLE_CUDA OFF)
endif()
