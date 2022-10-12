/***************************************************************************
                                ocl_texture.h
                             -------------------
                               W. Michael Brown

  Utilities for dealing with OpenCL textures

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Fri Jul 2 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef OCL_TEXTURE
#define OCL_TEXTURE

#include "ocl_kernel.h"
#include "ocl_mat.h"

namespace ucl_opencl {

/// Class storing a texture reference
class UCL_Texture {
 public:
  UCL_Texture() {}
  ~UCL_Texture() {}
  /// Construct with a specified texture reference
  inline UCL_Texture(UCL_Program &prog, const char *texture_name) { }
  /// Set the texture reference for this object
  inline void get_texture(UCL_Program &prog, const char *texture_name) { }

  /// Bind a float array where each fetch grabs a vector of length numel
  template<class mat_typ>
  inline void bind_float(mat_typ &vec, const unsigned numel) { }

  /// Unbind the texture reference from the memory allocation
  inline void unbind() { }

  /// Make a texture reference available to kernel
  inline void allow(UCL_Kernel &kernel) { }

 private:
  friend class UCL_Kernel;
};

/// Class storing a const global memory reference
class UCL_Const {
 public:
   UCL_Const() : _global_bytes(0), _active(false) {}
  ~UCL_Const() { clear(); }
  /// Construct with a specified global reference
  inline UCL_Const(UCL_Program &prog, const char *global_name)
    { get_global(prog,global_name); }
  /// Set the global reference for this object
  inline void get_global(UCL_Program &prog, const char *global_name) {
    if (_active) {
      CL_DESTRUCT_CALL(clReleaseContext(_context));
      CL_DESTRUCT_CALL(clReleaseCommandQueue(_cq));
    }
    _active = true;
    _context = prog._context;
    _cq = prog._cq;
    CL_SAFE_CALL(clRetainContext(_context));
    CL_SAFE_CALL(clRetainCommandQueue(_cq));
  }
  /// Copy from array on host to const memory
  template <class numtyp>
  inline void update_device(UCL_H_Vec<numtyp> &src, const int numel) {
    const int bytes=numel*sizeof(numtyp);
    if (_global_bytes < bytes) {
      if (_global_bytes) CL_SAFE_CALL(clReleaseMemObject(_global));
      cl_int e;
      _global = clCreateBuffer(_context, CL_MEM_READ_ONLY, bytes, NULL, &e);
      CL_SAFE_CALL(e);
    }
    CL_SAFE_CALL(clEnqueueWriteBuffer(_cq, _global, CL_FALSE, 0, bytes,
				      (void *)src.begin(), 0, NULL, NULL));
  }
  /// Get device ptr associated with object
  inline const cl_mem * begin() const { return &_global; }
  inline void clear() {
    if (_global_bytes) CL_SAFE_CALL(clReleaseMemObject(_global));
    if (_active) {
      CL_DESTRUCT_CALL(clReleaseContext(_context));
      CL_DESTRUCT_CALL(clReleaseCommandQueue(_cq));
    }
    _global_bytes=0;
    _active=false;
  }

 private:
  cl_mem _global;
  size_t _global_bytes;
  cl_context _context;
  cl_command_queue _cq;
  bool _active;
};

} // namespace

#endif

