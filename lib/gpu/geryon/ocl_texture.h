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

} // namespace

#endif

