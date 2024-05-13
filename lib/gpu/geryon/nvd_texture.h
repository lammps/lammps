/***************************************************************************
                                nvd_texture.h
                             -------------------
                               W. Michael Brown

  Utilities for dealing with CUDA Driver textures

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

#ifndef NVD_TEXTURE
#define NVD_TEXTURE

#include "nvd_kernel.h"
#include "nvd_mat.h"

namespace ucl_cudadr {

/// Class storing a texture reference
class UCL_Texture {
 public:
  UCL_Texture() {}
  ~UCL_Texture() {}
  /// Construct with a specified texture reference
  inline UCL_Texture(UCL_Program &prog, const char *texture_name)
    { get_texture(prog,texture_name); }
  /// Set the texture reference for this object
  inline void get_texture(UCL_Program &prog, const char *texture_name) {
    #if (CUDA_VERSION < 11000)
    CU_SAFE_CALL(cuModuleGetTexRef(&_tex, prog._module, texture_name));
    #endif
  }

  /// Bind a float array where each fetch grabs a vector of length numel
  template<class numtyp>
  inline void bind_float(UCL_D_Vec<numtyp> &vec, const unsigned numel)
    { _bind_float(vec,numel); }

  /// Bind a float array where each fetch grabs a vector of length numel
  template<class numtyp>
  inline void bind_float(UCL_D_Mat<numtyp> &vec, const unsigned numel)
    { _bind_float(vec,numel); }

  /// Bind a float array where each fetch grabs a vector of length numel
  template<class numtyp, class devtyp>
  inline void bind_float(UCL_Vector<numtyp, devtyp> &vec, const unsigned numel)
    { _bind_float(vec.device,numel); }

  /// Bind a float array where each fetch grabs a vector of length numel
  template<class numtyp, class devtyp>
  inline void bind_float(UCL_Matrix<numtyp, devtyp> &vec, const unsigned numel)
    { _bind_float(vec.device,numel); }

  /// Unbind the texture reference from the memory allocation
  inline void unbind() { }

  /// Make a texture reference available to kernel
  inline void allow(UCL_Kernel &kernel) {
  }

 private:
  #if (CUDA_VERSION < 11000)
  CUtexref _tex;
  #endif
  friend class UCL_Kernel;

  template<class mat_typ>
  inline void _bind_float(mat_typ &vec, const unsigned numel) {
    #if (CUDA_VERSION < 11000)
    #ifdef UCL_DEBUG
    assert(numel!=0 && numel<5);
    #endif
    CU_SAFE_CALL(cuTexRefSetAddress(nullptr, _tex, vec.cbegin(),
                 vec.numel()*vec.element_size()));
    if (vec.element_size()==sizeof(float))
      CU_SAFE_CALL(cuTexRefSetFormat(_tex, CU_AD_FORMAT_FLOAT, numel));
    else {
      if (numel>2)
        CU_SAFE_CALL(cuTexRefSetFormat(_tex, CU_AD_FORMAT_SIGNED_INT32, numel));
      else
        CU_SAFE_CALL(cuTexRefSetFormat(_tex,CU_AD_FORMAT_SIGNED_INT32,numel*2));
    }
    #endif
  }

};

/// Class storing a const global memory reference
class UCL_Const {
 public:
  UCL_Const() {}
  ~UCL_Const() {}
  /// Construct with a specified global reference
  inline UCL_Const(UCL_Program &prog, const char *global_name)
    { get_global(prog,global_name); }
  /// Set the global reference for this object
  inline void get_global(UCL_Program &prog, const char *global_name) {
    _cq=prog.cq();
    CU_SAFE_CALL(cuModuleGetGlobal(&_global, &_global_bytes, prog._module,
                                   global_name));
  }
  /// Copy from array on host to const memory
  template <class numtyp>
  inline void update_device(UCL_H_Vec<numtyp> &src, const int numel) {
    CU_SAFE_CALL(cuMemcpyHtoDAsync(_global, src.begin(), numel*sizeof(numtyp),
                                   _cq));
  }
  /// Get device ptr associated with object
  inline const CUdeviceptr * begin() const { return &_global; }
  inline void clear() {}

 private:
  CUstream _cq;
  CUdeviceptr _global;
  size_t _global_bytes;
  friend class UCL_Kernel;
};

} // namespace

#endif

