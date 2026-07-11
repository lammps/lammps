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

#ifdef __HIP_PLATFORM_NVCC__
inline hipError_t hipModuleGetTexRef(CUtexref* texRef, hipModule_t hmod, const char* name){
  return hipCUResultTohipError(cuModuleGetTexRef(texRef, hmod, name));
}
inline hipError_t hipTexRefSetFormat(CUtexref tex, hipArray_Format fmt, int NumPackedComponents) {
    return hipCUResultTohipError(cuTexRefSetFormat(tex, (CUarray_format)fmt, NumPackedComponents ));
}
inline hipError_t hipTexRefSetAddress(size_t* offset, CUtexref tex, hipDeviceptr_t devPtr, size_t size) {
    return hipCUResultTohipError(cuTexRefSetAddress(offset, tex, devPtr, size));
}
#endif

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
  #ifdef __HIP_PLATFORM_NVCC__
      CU_SAFE_CALL(hipModuleGetTexRef(&_tex, prog._module, texture_name));
  #else
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
#ifdef __HIP_PLATFORM_NVCC__
  CUtexref _tex;
#else
  void* _device_ptr_to_global_var;
#endif
  friend class UCL_Kernel;

  template<class mat_typ>
  inline void _bind_float(mat_typ &vec, const unsigned numel) {
    #ifdef UCL_DEBUG
    assert(numel!=0 && numel<5);
    #endif

#ifdef __HIP_PLATFORM_NVCC__
    if (vec.element_size()==sizeof(float))
      CU_SAFE_CALL(hipTexRefSetFormat(_tex, HIP_AD_FORMAT_FLOAT, numel));
    else {
      if (numel>2)
        CU_SAFE_CALL(hipTexRefSetFormat(_tex, HIP_AD_FORMAT_SIGNED_INT32, numel));
      else
        CU_SAFE_CALL(hipTexRefSetFormat(_tex,HIP_AD_FORMAT_SIGNED_INT32,numel*2));
    }
    CU_SAFE_CALL(hipTexRefSetAddress(nullptr, _tex, vec.cbegin(), vec.numel()*vec.element_size()));
#else
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

