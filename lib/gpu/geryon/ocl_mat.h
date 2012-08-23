/***************************************************************************
                                  ocl_mat.h
                             -------------------
                               W. Michael Brown

  OpenCL Specific Vector/Matrix Containers, Memory Management, and I/O

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Wed Jan 13 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

/*! \file */
   
#ifndef OCL_MAT_H
#define OCL_MAT_H

#include "ocl_memory.h"

/// Namespace for OpenCL routines
namespace ucl_opencl {

#define _UCL_MAT_ALLOW
#define _UCL_DEVICE_PTR_MAT
#define _OCL_MAT
#include "ucl_basemat.h"
#include "ucl_h_vec.h"
#include "ucl_h_mat.h"
#include "ucl_d_vec.h"
#include "ucl_d_mat.h"
#include "ucl_s_obj_help.h"
#include "ucl_vector.h"
#include "ucl_matrix.h"
#undef _UCL_DEVICE_PTR_MAT
#undef _OCL_MAT
#undef _UCL_MAT_ALLOW

#define UCL_COPY_ALLOW
#include "ucl_copy.h"
#undef UCL_COPY_ALLOW

#define UCL_PRINT_ALLOW
#include "ucl_print.h"
#undef UCL_PRINT_ALLOW

} // namespace ucl_cudart 

#endif
