/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef HIP_MEMORY_H
#define HIP_MEMORY_H


#include <hip/hip_runtime.h>
#include <iostream>
#include <cassert>
#include <cstring>
#include "hip_macros.h"
#include "hip_device.h"
#include "ucl_types.h"

namespace ucl_hip {

// --------------------------------------------------------------------------
// - API Specific Types
// --------------------------------------------------------------------------
//typedef dim3 ucl_kernel_dim;

#ifdef __HIP_PLATFORM_NVCC__
typedef enum hipArray_Format {
    HIP_AD_FORMAT_UNSIGNED_INT8 = 0x01,
    HIP_AD_FORMAT_UNSIGNED_INT16 = 0x02,
    HIP_AD_FORMAT_UNSIGNED_INT32 = 0x03,
    HIP_AD_FORMAT_SIGNED_INT8 = 0x08,
    HIP_AD_FORMAT_SIGNED_INT16 = 0x09,
    HIP_AD_FORMAT_SIGNED_INT32 = 0x0a,
    HIP_AD_FORMAT_HALF = 0x10,
    HIP_AD_FORMAT_FLOAT = 0x20
}hipArray_Format;
#endif

// --------------------------------------------------------------------------
// - API SPECIFIC DEVICE POINTERS
// --------------------------------------------------------------------------
typedef hipDeviceptr_t device_ptr;

// --------------------------------------------------------------------------
// - HOST MEMORY ALLOCATION ROUTINES
// --------------------------------------------------------------------------
template <class mat_type, class copy_type>
inline int _host_alloc(mat_type &mat, copy_type &cm, const size_t n,
                       const enum UCL_MEMOPT kind, const enum UCL_MEMOPT kind2){
  hipError_t err=hipSuccess;
  if (kind==UCL_NOT_PINNED)
    *(mat.host_ptr())=(typename mat_type::data_type*)malloc(n);
  else if (kind==UCL_WRITE_ONLY)
    err=hipHostMalloc((void **)mat.host_ptr(),n,hipHostMallocWriteCombined);
  else
    err=hipHostMalloc((void **)mat.host_ptr(),n,hipHostMallocDefault);
  if (err!=hipSuccess || *(mat.host_ptr())==NULL)
    return UCL_MEMORY_ERROR;
  mat.cq()=cm.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _host_alloc(mat_type &mat, UCL_Device &dev, const size_t n,
                       const enum UCL_MEMOPT kind, const enum UCL_MEMOPT kind2){
  hipError_t err=hipSuccess;
  if (kind==UCL_NOT_PINNED)
    *(mat.host_ptr())=(typename mat_type::data_type*)malloc(n);
  else if (kind==UCL_WRITE_ONLY)
    err=hipHostMalloc((void **)mat.host_ptr(),n,hipHostMallocWriteCombined);
  else
    err=hipHostMalloc((void **)mat.host_ptr(),n,hipHostMallocDefault);
  if (err!=hipSuccess || *(mat.host_ptr())==NULL)
    return UCL_MEMORY_ERROR;
  mat.cq()=dev.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline void _host_free(mat_type &mat) {
  if (mat.kind()==UCL_VIEW)
    return;
  else if (mat.kind()!=UCL_NOT_PINNED)
    CU_DESTRUCT_CALL(hipHostFree(mat.begin()));
  else
    free(mat.begin());
}

template <class mat_type>
inline int _host_resize(mat_type &mat, const size_t n) {
  _host_free(mat);
  hipError_t err=hipSuccess;
  if (mat.kind()==UCL_NOT_PINNED)
    *(mat.host_ptr())=(typename mat_type::data_type*)malloc(n);
  else if (mat.kind()==UCL_WRITE_ONLY)
    err=hipHostMalloc((void **)mat.host_ptr(),n,hipHostMallocWriteCombined);
  else
    err=hipHostMalloc((void **)mat.host_ptr(),n,hipHostMallocDefault);
  if (err!=hipSuccess || *(mat.host_ptr())==NULL)
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}

// --------------------------------------------------------------------------
// - DEVICE MEMORY ALLOCATION ROUTINES
// --------------------------------------------------------------------------
template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, copy_type &cm, const size_t n,
                         const enum UCL_MEMOPT kind) {
  hipError_t err=hipMalloc((void**)&mat.cbegin(),n);
  if (err!=hipSuccess)
    return UCL_MEMORY_ERROR;
  mat.cq()=cm.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _device_alloc(mat_type &mat, UCL_Device &dev, const size_t n,
                         const enum UCL_MEMOPT kind) {
  hipError_t err=hipMalloc((void**)&mat.cbegin(),n);
  if (err!=hipSuccess)
    return UCL_MEMORY_ERROR;
  mat.cq()=dev.cq();
  return UCL_SUCCESS;
}

template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, copy_type &cm, const size_t rows,
                         const size_t cols, size_t &pitch,
                         const enum UCL_MEMOPT kind) {
  hipError_t err;
  size_t upitch;
  err=hipMallocPitch((void**)&mat.cbegin(),&upitch,
                      cols*sizeof(typename mat_type::data_type),rows);
  pitch=static_cast<size_t>(upitch);
  if (err!=hipSuccess)
    return UCL_MEMORY_ERROR;
  mat.cq()=cm.cq();
  return UCL_SUCCESS;
}

template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, UCL_Device &d, const size_t rows,
                         const size_t cols, size_t &pitch,
                         const enum UCL_MEMOPT kind) {
  hipError_t err;
  size_t upitch;
  err=hipMallocPitch((void**)&mat.cbegin(),&upitch,
                      cols*sizeof(typename mat_type::data_type),rows);
  pitch=static_cast<size_t>(upitch);
  if (err!=hipSuccess)
    return UCL_MEMORY_ERROR;
  mat.cq()=d.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline void _device_free(mat_type &mat) {
  if (mat.kind()!=UCL_VIEW){
    CU_DESTRUCT_CALL(hipFree((void*)mat.cbegin()));
  }
}

template <class mat_type>
inline int _device_resize(mat_type &mat, const size_t n) {
  _device_free(mat);
  hipError_t err=hipMalloc((void**)&mat.cbegin(),n);
  if (err!=hipSuccess)
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _device_resize(mat_type &mat, const size_t rows,
                          const size_t cols, size_t &pitch) {
  _device_free(mat);
  hipError_t err;
  size_t upitch;
  err=hipMallocPitch((void**)&mat.cbegin(),&upitch,
                      cols*sizeof(typename mat_type::data_type),rows);
  pitch=static_cast<size_t>(upitch);
  if (err!=hipSuccess)
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}

inline void _device_view(hipDeviceptr_t *ptr, hipDeviceptr_t &in) {
  *ptr=in;
}

template <class numtyp>
inline void _device_view(hipDeviceptr_t *ptr, numtyp *in) {
  *ptr=0;
}

inline void _device_view(hipDeviceptr_t *ptr, hipDeviceptr_t &in,
                         const size_t offset, const size_t numsize) {
  *ptr=(hipDeviceptr_t)(((char*)in)+offset*numsize);
}

template <class numtyp>
inline void _device_view(hipDeviceptr_t *ptr, numtyp *in,
                         const size_t offset, const size_t numsize) {
  *ptr=0;
}

// --------------------------------------------------------------------------
// - DEVICE IMAGE ALLOCATION ROUTINES
// --------------------------------------------------------------------------
template <class mat_type, class copy_type>
inline void _device_image_alloc(mat_type &mat, copy_type &cm, const size_t rows,
                                const size_t cols) {
  assert(0==1);
}

template <class mat_type, class copy_type>
inline void _device_image_alloc(mat_type &mat, UCL_Device &d, const size_t rows,
                                const size_t cols) {
  assert(0==1);
}

template <class mat_type>
inline void _device_image_free(mat_type &mat) {
  assert(0==1);
}

// --------------------------------------------------------------------------
// - ZERO ROUTINES
// --------------------------------------------------------------------------
inline void _host_zero(void *ptr, const size_t n) {
  memset(ptr,0,n);
}

template <class mat_type>
inline void _device_zero(mat_type &mat, const size_t n, command_queue &cq) {
    CU_SAFE_CALL(hipMemsetAsync((void*)mat.cbegin(),0,n,cq));
}


// --------------------------------------------------------------------------
// - MEMCPY ROUTINES
// --------------------------------------------------------------------------


template<class mat1, class mat2>
hipMemcpyKind _memcpy_kind(mat1 &dst, const mat2 &src){
  assert(mat1::MEM_TYPE < 2 && mat2::MEM_TYPE < 2);
  return (hipMemcpyKind)((1 - mat2::MEM_TYPE)*2 + (1 - mat1::MEM_TYPE));
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const mat2 &src, const size_t n) {
  CU_SAFE_CALL(hipMemcpy((void*)dst.begin(), (void*)src.begin(), n, _memcpy_kind(dst, src)));
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const mat2 &src, const size_t n, hipStream_t &cq) {
  CU_SAFE_CALL(hipMemcpyAsync((void*)dst.begin(), (void*)src.begin(), n, _memcpy_kind(dst, src), cq));
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const size_t dpitch, const mat2 &src,
                       const size_t spitch, const size_t cols,
                       const size_t rows) {
  CU_SAFE_CALL(hipMemcpy2D((void*)dst.begin(), dpitch, (void*)src.begin(), spitch, cols, rows, _memcpy_kind(dst, src)));
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const size_t dpitch, const mat2 &src,
                       const size_t spitch, const size_t cols,
                       const size_t rows,hipStream_t &cq) {
  CU_SAFE_CALL(hipMemcpy2DAsync((void*)dst.begin(), dpitch, (void*)src.begin(), spitch, cols, rows, _memcpy_kind(dst, src), cq));
}

} // namespace ucl_cudart

#endif

