/***************************************************************************
                                nvd_memory.h
                             -------------------
                               W. Michael Brown

  CUDA Driver Specific Memory Management and Vector/Matrix Containers

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Thu Jan 21 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

#ifndef NVD_MEMORY_H
#define NVD_MEMORY_H

#include <iostream>
#include <cassert>
#include <cstring>
#include "nvd_macros.h"
#include "ucl_types.h"

namespace ucl_cudadr {

// --------------------------------------------------------------------------
// - API Specific Types
// --------------------------------------------------------------------------
//typedef dim3 ucl_kernel_dim;

// --------------------------------------------------------------------------
// - API SPECIFIC DEVICE POINTERS
// --------------------------------------------------------------------------
typedef CUdeviceptr device_ptr;

// --------------------------------------------------------------------------
// - HOST MEMORY ALLOCATION ROUTINES
// --------------------------------------------------------------------------
template <class mat_type, class copy_type>
inline int _host_alloc(mat_type &mat, copy_type &cm, const size_t n,
                       const enum UCL_MEMOPT kind, const enum UCL_MEMOPT kind2){
  CUresult err=CUDA_SUCCESS;
  if (kind==UCL_NOT_PINNED)
    *(mat.host_ptr())=(typename mat_type::data_type*)malloc(n);
  else if (kind==UCL_WRITE_ONLY)
    err=cuMemHostAlloc((void **)mat.host_ptr(),n,CU_MEMHOSTALLOC_WRITECOMBINED);
  else
    err=cuMemAllocHost((void **)mat.host_ptr(),n);
  if (err!=CUDA_SUCCESS || *(mat.host_ptr())==nullptr)
    return UCL_MEMORY_ERROR;
  mat.cq()=cm.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _host_alloc(mat_type &mat, UCL_Device &dev, const size_t n,
                       const enum UCL_MEMOPT kind, const enum UCL_MEMOPT kind2){
  CUresult err=CUDA_SUCCESS;
  if (kind==UCL_NOT_PINNED)
    *(mat.host_ptr())=(typename mat_type::data_type*)malloc(n);
  else if (kind==UCL_WRITE_ONLY)
    err=cuMemHostAlloc((void **)mat.host_ptr(),n,CU_MEMHOSTALLOC_WRITECOMBINED);
  else
    err=cuMemAllocHost((void **)mat.host_ptr(),n);
  if (err!=CUDA_SUCCESS || *(mat.host_ptr())==nullptr)
    return UCL_MEMORY_ERROR;
  mat.cq()=dev.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline void _host_free(mat_type &mat) {
  if (mat.kind()==UCL_VIEW)
    return;
  else if (mat.kind()!=UCL_NOT_PINNED)
    CU_DESTRUCT_CALL(cuMemFreeHost(mat.begin()));
  else
    free(mat.begin());
}

template <class mat_type>
inline int _host_resize(mat_type &mat, const size_t n) {
  _host_free(mat);
  CUresult err=CUDA_SUCCESS;
  if (mat.kind()==UCL_NOT_PINNED)
    *(mat.host_ptr())=(typename mat_type::data_type*)malloc(n);
  else if (mat.kind()==UCL_WRITE_ONLY)
    err=cuMemHostAlloc((void **)mat.host_ptr(),n,CU_MEMHOSTALLOC_WRITECOMBINED);
  else
    err=cuMemAllocHost((void **)mat.host_ptr(),n);
  if (err!=CUDA_SUCCESS || *(mat.host_ptr())==nullptr)
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}

// --------------------------------------------------------------------------
// - DEVICE MEMORY ALLOCATION ROUTINES
// --------------------------------------------------------------------------
template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, copy_type &cm, const size_t n,
                         const enum UCL_MEMOPT kind) {
  CUresult err=cuMemAlloc(&mat.cbegin(),n);
  if (err!=CUDA_SUCCESS)
    return UCL_MEMORY_ERROR;
  mat.cq()=cm.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _device_alloc(mat_type &mat, UCL_Device &dev, const size_t n,
                         const enum UCL_MEMOPT kind) {
  CUresult err=cuMemAlloc(&mat.cbegin(),n);
  if (err!=CUDA_SUCCESS)
    return UCL_MEMORY_ERROR;
  mat.cq()=dev.cq();
  return UCL_SUCCESS;
}

template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, copy_type &cm, const size_t rows,
                         const size_t cols, size_t &pitch,
                         const enum UCL_MEMOPT kind) {
  CUresult err;
  CUDA_INT_TYPE upitch;
  err=cuMemAllocPitch(&mat.cbegin(),&upitch,
                      cols*sizeof(typename mat_type::data_type),rows,16);
  pitch=static_cast<size_t>(upitch);
  if (err!=CUDA_SUCCESS)
    return UCL_MEMORY_ERROR;
  mat.cq()=cm.cq();
  return UCL_SUCCESS;
}

template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, UCL_Device &d, const size_t rows,
                         const size_t cols, size_t &pitch,
                         const enum UCL_MEMOPT kind) {
  CUresult err;
  unsigned upitch;
  err=cuMemAllocPitch(&mat.cbegin(),&upitch,
                      cols*sizeof(typename mat_type::data_type),rows,16);
  pitch=static_cast<size_t>(upitch);
  if (err!=CUDA_SUCCESS)
    return UCL_MEMORY_ERROR;
  mat.cq()=d.cq();
  return UCL_SUCCESS;
}

template <class mat_type>
inline void _device_free(mat_type &mat) {
  if (mat.kind()!=UCL_VIEW)
    CU_DESTRUCT_CALL(cuMemFree(mat.cbegin()));
}

template <class mat_type>
inline int _device_resize(mat_type &mat, const size_t n) {
  _device_free(mat);
  CUresult err=cuMemAlloc(&mat.cbegin(),n);
  if (err!=CUDA_SUCCESS)
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _device_resize(mat_type &mat, const size_t rows,
                          const size_t cols, size_t &pitch) {
  _device_free(mat);
  CUresult err;
  CUDA_INT_TYPE upitch;
  err=cuMemAllocPitch(&mat.cbegin(),&upitch,
                      cols*sizeof(typename mat_type::data_type),rows,16);
  pitch=static_cast<size_t>(upitch);
  if (err!=CUDA_SUCCESS)
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}

inline void _device_view(CUdeviceptr *ptr, CUdeviceptr &in) {
  *ptr=in;
}

template <class numtyp>
inline void _device_view(CUdeviceptr *ptr, numtyp *in) {
  *ptr=0;
}

inline void _device_view(CUdeviceptr *ptr, CUdeviceptr &in,
                         const size_t offset, const size_t numsize) {
  *ptr=in+offset*numsize;
}

template <class numtyp>
inline void _device_view(CUdeviceptr *ptr, numtyp *in,
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
  if (n%32==0)
    CU_SAFE_CALL(cuMemsetD32Async(mat.cbegin(),0,n/4,cq));
  else if (n%16==0)
    CU_SAFE_CALL(cuMemsetD16Async(mat.cbegin(),0,n/2,cq));
  else
    CU_SAFE_CALL(cuMemsetD8Async(mat.cbegin(),0,n,cq));
}

// --------------------------------------------------------------------------
// - HELPER FUNCTIONS FOR MEMCPY ROUTINES
// --------------------------------------------------------------------------

inline void _nvd_set_2D_loc(CUDA_MEMCPY2D &ins, const size_t dpitch,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
  ins.srcXInBytes=0;
  ins.srcY=0;
  ins.srcPitch=spitch;
  ins.dstXInBytes=0;
  ins.dstY=0;
  ins.dstPitch=dpitch;
  ins.WidthInBytes=cols;
  ins.Height=rows;
}

template <int mem> struct _nvd_set_2D_mem;
template <> struct _nvd_set_2D_mem<1>
  { static CUmemorytype a() { return CU_MEMORYTYPE_HOST; } };
template <> struct _nvd_set_2D_mem<2>
  { static CUmemorytype a() { return CU_MEMORYTYPE_ARRAY; } };
template <int mem> struct _nvd_set_2D_mem
  { static CUmemorytype a() { return CU_MEMORYTYPE_DEVICE; } };


// --------------------------------------------------------------------------
// - MEMCPY ROUTINES
// --------------------------------------------------------------------------

template<int mem1, int mem2> struct _ucl_memcpy;

// Both are images
template<> struct _ucl_memcpy<2,2> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    assert(0==1);
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstArray=dst.cbegin();
    ins.srcArray=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstArray=dst.cbegin();
    ins.srcArray=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Destination is texture, source on device
template<> struct _ucl_memcpy<2,0> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    assert(0==1);
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstArray=dst.cbegin();
    ins.srcDevice=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstArray=dst.cbegin();
    ins.srcDevice=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Destination is texture, source on host
template<> struct _ucl_memcpy<2,1> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    assert(0==1);
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstArray=dst.cbegin();
    ins.srcHost=src.begin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstArray=dst.cbegin();
    ins.srcHost=src.begin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Source is texture, dest on device
template<> struct _ucl_memcpy<0,2> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    assert(0==1);
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstDevice=dst.cbegin();
    ins.srcArray=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstDevice=dst.cbegin();
    ins.srcArray=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Source is texture, dest on host
template<> struct _ucl_memcpy<1,2> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    assert(0==1);
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstHost=dst.begin();
    ins.srcArray=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstHost=dst.begin();
    ins.srcArray=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Neither are textures, destination on host
template <> struct _ucl_memcpy<1,0> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    CU_SAFE_CALL(cuMemcpyDtoH(dst.begin(),src.cbegin(),n));
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    CU_SAFE_CALL(cuMemcpyDtoHAsync(dst.begin(),src.cbegin(),n,cq));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstHost=dst.begin();
    ins.srcDevice=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstHost=dst.begin();
    ins.srcDevice=src.cbegin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Neither are textures, source on host
template <> struct _ucl_memcpy<0,1> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    CU_SAFE_CALL(cuMemcpyHtoD(dst.cbegin(),src.begin(),n));
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    CU_SAFE_CALL(cuMemcpyHtoDAsync(dst.cbegin(),src.begin(),n,cq));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstDevice=dst.cbegin();
    ins.srcHost=src.begin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstDevice=dst.cbegin();
    ins.srcHost=src.begin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Neither are textures, both on host
template <> struct _ucl_memcpy<1,1> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n)
    { memcpy(dst.begin(),src.begin(),n); }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq)
    { memcpy(dst.begin(),src.begin(),n); }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstHost=dst.begin();
    ins.srcHost=src.begin();
    CU_SAFE_CALL(cuMemcpy2D(&ins));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    CUDA_MEMCPY2D ins;
    _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
    ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
    ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
    ins.dstHost=dst.begin();
    ins.srcHost=src.begin();
    CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
  }
};

// Neither are textures, both on device
template <int mem1, int mem2> struct _ucl_memcpy {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n) {
    CU_SAFE_CALL(cuMemcpyDtoD(dst.cbegin(),src.cbegin(),n));
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        CUstream &cq) {
    CU_SAFE_CALL(cuMemcpyDtoDAsync(dst.cbegin(),src.cbegin(),n,cq));
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows) {
    if (p1::PADDED==0 || p2::PADDED==0) {
      size_t src_offset=0, dst_offset=0;
      for (size_t i=0; i<rows; i++) {
        CU_SAFE_CALL(cuMemcpyDtoD(dst.cbegin()+dst_offset,
                                  src.cbegin()+src_offset,cols));
        src_offset+=spitch;
        dst_offset+=dpitch;
      }
    } else {
      CUDA_MEMCPY2D ins;
      _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
      ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
      ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
      ins.dstDevice=dst.cbegin();
      ins.srcDevice=src.cbegin();
      CU_SAFE_CALL(cuMemcpy2D(&ins));
    }
  }
  template <class p1, class p2>
      static inline void mc(p1 &dst, const size_t dpitch, const p2 &src,
                            const size_t spitch, const size_t cols,
                            const size_t rows, CUstream &cq) {
    if (p1::PADDED==0 || p2::PADDED==0) {
      size_t src_offset=0, dst_offset=0;
      for (size_t i=0; i<rows; i++) {
        CU_SAFE_CALL(cuMemcpyDtoDAsync(dst.cbegin()+dst_offset,
                                       src.cbegin()+src_offset,cols,cq));
        src_offset+=spitch;
        dst_offset+=dpitch;
      }
    } else {
      CUDA_MEMCPY2D ins;
      _nvd_set_2D_loc(ins,dpitch,spitch,cols,rows);
      ins.dstMemoryType=_nvd_set_2D_mem<p1::MEM_TYPE>::a();
      ins.srcMemoryType=_nvd_set_2D_mem<p2::MEM_TYPE>::a();
      ins.dstDevice=dst.cbegin();
      ins.srcDevice=src.cbegin();
      CU_SAFE_CALL(cuMemcpy2DAsync(&ins,cq));
    }
  }
};

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const mat2 &src, const size_t n) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,src,n);
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const mat2 &src, const size_t n,
                       CUstream &cq) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,src,n,cq);
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const size_t dpitch, const mat2 &src,
                       const size_t spitch, const size_t cols,
                       const size_t rows) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,dpitch,src,spitch,cols,
                                                 rows);
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const size_t dpitch, const mat2 &src,
                       const size_t spitch, const size_t cols,
                       const size_t rows,CUstream &cq) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,dpitch,src,spitch,cols,
                                                 rows,cq);
}

} // namespace ucl_cudart

#endif

