/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef HIP_TEXTURE
#define HIP_TEXTURE


#include <hip/hip_runtime.h>
#include "hip_kernel.h"
#include "hip_mat.h"

namespace ucl_hip {

// With NVIDIA GPUs the kernels are compiled with CUDA 11 or later and then do
// not use texture references, but access the data directly (same as with the
// CUDA back end), so there is nothing to look up or bind on the host side.

/// Class storing a texture reference
class UCL_Texture {
 public:
  UCL_Texture() {}
  ~UCL_Texture() {}
  /// Construct with a specified texture reference
  inline UCL_Texture(UCL_Program &prog, const char *texture_name)
    { get_texture(prog,texture_name); }
  /// Set the texture reference for this object
  inline void get_texture(UCL_Program &prog, const char *texture_name)
    {
  #if !defined(__HIP_PLATFORM_NVCC__) && !defined(__HIP_PLATFORM_NVIDIA__)
      size_t _global_var_size;
      CU_SAFE_CALL(hipModuleGetGlobal(&_device_ptr_to_global_var, &_global_var_size, prog._module, texture_name));
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
  inline void allow(UCL_Kernel &) {
  }

 private:
#if !defined(__HIP_PLATFORM_NVCC__) && !defined(__HIP_PLATFORM_NVIDIA__)
  void* _device_ptr_to_global_var;
#endif
  friend class UCL_Kernel;

  template<class mat_typ>
  inline void _bind_float(mat_typ &vec, const unsigned numel) {
    #ifdef UCL_DEBUG
    assert(numel!=0 && numel<5);
    #endif

#if !defined(__HIP_PLATFORM_NVCC__) && !defined(__HIP_PLATFORM_NVIDIA__)
    void* data_ptr = (void*)vec.cbegin();
    CU_SAFE_CALL(hipMemcpyHtoD(hipDeviceptr_t(_device_ptr_to_global_var), &data_ptr, sizeof(void*)));
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
    CU_SAFE_CALL(hipModuleGetGlobal(&_global, &_global_bytes, prog._module,
                                    global_name));
  }
  /// Copy from array on host to const memory
  template <class numtyp>
  inline void update_device(UCL_H_Vec<numtyp> &src, const int numel) {
    CU_SAFE_CALL(hipMemcpyHtoDAsync(_global, src.begin(), numel*sizeof(numtyp),
                                    _cq));
  }
  /// Get device ptr associated with object
  inline const hipDeviceptr_t * begin() const { return &_global; }
  inline void clear() {}

 private:
  hipStream_t _cq;
  hipDeviceptr_t _global;
  size_t _global_bytes;
  friend class UCL_Kernel;
};

} // namespace

#endif

