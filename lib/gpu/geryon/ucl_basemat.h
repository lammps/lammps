/***************************************************************************
                                 ucl_basemat.h
                             -------------------
                               W. Michael Brown

  Vector/Matrix Base Container

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Thu Jun 25 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

// Only allow this file to be included by CUDA and OpenCL specific headers
#ifdef _UCL_MAT_ALLOW

#include "ucl_types.h"

#define UCL_H_VecT UCL_H_Vec<numtyp>
#define UCL_H_VecD UCL_H_Vec<double>
#define UCL_H_VecS UCL_H_Vec<float>
#define UCL_H_VecI UCL_H_Vec<int>

#define UCL_D_VecT UCL_D_Vec<numtyp>
#define UCL_D_VecD UCL_D_Vec<double>
#define UCL_D_VecS UCL_D_Vec<float>
#define UCL_D_VecI UCL_D_Vec<int>
#define UCL_D_VecI2 UCL_D_Vec<int2>
#define UCL_D_VecU2 UCL_D_Vec<uint2>

#define UCL_D_MatT UCL_D_Mat<numtyp>
#define UCL_D_MatD UCL_D_Mat<double>
#define UCL_D_MatS UCL_D_Mat<float>
#define UCL_D_MatI UCL_D_Mat<int>

#define UCL_ConstMatT UCL_ConstMat<numtyp>
#define UCL_ConstMatD UCL_ConstMat<double>
#define UCL_ConstMatS UCL_ConstMat<float>
#define UCL_ConstMatI UCL_ConstMat<int>
#define UCL_ConstMatD2 UCL_ConstMat<double2>

/// Base class for vector/matrix containers
/** All containers are associated with a default command queue.
  * For CUDA, this is the default stream.
  * 
  * The default queue is used for asynchonrous operations on the container 
  * that do not specify a queue. For OpenCL, this queue is also used in
  * calls for reserving and copying memory **/ 
class UCL_BaseMat {
 public:
  UCL_BaseMat() : _cq(0) { }
  virtual ~UCL_BaseMat() { }
  /// Return the default command queue/stream associated with this data
  inline command_queue & cq() { return _cq; }
  /// Block until command_queue associated with matrix is complete
  inline void sync() { ucl_sync(_cq); }
  
  #ifdef UCL_DEBUG
  // Returns the type of host allocation
  virtual inline enum UCL_MEMOPT kind() const { return UCL_NOT_PINNED; }
  #endif 
 protected:
  command_queue _cq;
};

#endif

