/***************************************************************************
                                ocl_memory.h
                             -------------------
                               W. Michael Brown

  OpenCL Specific Memory Management and Vector/Matrix Containers

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

#ifndef OCL_MEMORY_H
#define OCL_MEMORY_H

#include <iostream>
#include <cassert>
#include <cstring>
#include "ucl_types.h"

namespace ucl_opencl {

// --------------------------------------------------------------------------
// - API Specific Types
// --------------------------------------------------------------------------
struct ocl_kernel_dim {
  size_t x,y,z;
  ocl_kernel_dim(size_t _x = 1, size_t _y = 1, size_t _z = 1) : 
    x(_x), y(_y), z(_z) {}
  operator size_t * () { return (size_t *)this; }
  operator const size_t * () const { return (const size_t *)this; } 
};
typedef ocl_kernel_dim ucl_kernel_dim;

// --------------------------------------------------------------------------
// - API SPECIFIC DEVICE POINTERS
// --------------------------------------------------------------------------
typedef cl_mem device_ptr;

// --------------------------------------------------------------------------
// - HOST MEMORY ALLOCATION ROUTINES
// --------------------------------------------------------------------------

template <class mat_type, class copy_type>
inline int _host_alloc(mat_type &mat, copy_type &cm, const size_t n,  
                       const enum UCL_MEMOPT kind, const enum UCL_MEMOPT kind2){
  cl_int error_flag;
  cl_context context;
  CL_SAFE_CALL(clGetMemObjectInfo(cm.cbegin(),CL_MEM_CONTEXT,sizeof(context),
                                  &context,NULL));
  
  cl_mem_flags buffer_perm;
  cl_map_flags map_perm;
  if (kind2==UCL_NOT_SPECIFIED) {
    if (kind==UCL_READ_ONLY) {
      #ifdef CL_VERSION_1_2
      buffer_perm=CL_MEM_HOST_READ_ONLY|CL_MEM_WRITE_ONLY|CL_MEM_ALLOC_HOST_PTR;
      #else
      buffer_perm=CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR;
      #endif
      map_perm=CL_MAP_READ;
    } else if (kind==UCL_WRITE_ONLY) {
      #ifdef CL_VERSION_1_2
      buffer_perm=CL_MEM_HOST_WRITE_ONLY|CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR;
      #else
      buffer_perm=CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR;
      #endif
      map_perm=CL_MAP_WRITE;
    } else {
      buffer_perm=CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR;
      map_perm=CL_MAP_READ | CL_MAP_WRITE;
    }
  } else {
    if (kind2==UCL_READ_ONLY)
      buffer_perm=CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR;
    else if (kind2==UCL_WRITE_ONLY)
      buffer_perm=CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR;
    else
      buffer_perm=CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR;
    
    if (kind==UCL_READ_ONLY) {
      #ifdef CL_VERSION_1_2
      buffer_perm=buffer_perm | CL_MEM_HOST_READ_ONLY;
      #endif
      map_perm=CL_MAP_READ;
    } else if (kind==UCL_WRITE_ONLY) {
      #ifdef CL_VERSION_1_2
      buffer_perm=buffer_perm | CL_MEM_HOST_WRITE_ONLY;
      #endif
      map_perm=CL_MAP_WRITE;
    } else
      map_perm=CL_MAP_READ | CL_MAP_WRITE;
  }
    
  mat.cbegin()=clCreateBuffer(context,buffer_perm,n,NULL,&error_flag);
  if (error_flag != CL_SUCCESS) 
    return UCL_MEMORY_ERROR;
    *mat.host_ptr() = (typename mat_type::data_type*)
                      clEnqueueMapBuffer(cm.cq(),mat.cbegin(),CL_TRUE,
                                         map_perm,0,n,0,NULL,NULL,NULL);

  mat.cq()=cm.cq();
  CL_SAFE_CALL(clRetainCommandQueue(mat.cq()));
  return UCL_SUCCESS;
}

template <class mat_type, class copy_type>
inline int _host_view(mat_type &mat, copy_type &cm, const size_t n) {
  cl_int error_flag;
  cl_context context;
  CL_SAFE_CALL(clGetMemObjectInfo(cm.cbegin(),CL_MEM_CONTEXT,sizeof(context),
                                  &context,NULL));
  cl_mem_flags orig_flags;
  CL_SAFE_CALL(clGetMemObjectInfo(cm.cbegin(),CL_MEM_FLAGS,sizeof(orig_flags),
                                  &orig_flags,NULL));
  orig_flags=orig_flags & ~CL_MEM_ALLOC_HOST_PTR;
  
  mat.cbegin()=clCreateBuffer(context, CL_MEM_USE_HOST_PTR | orig_flags, n,
                              *mat.host_ptr(), &error_flag);

  CL_CHECK_ERR(error_flag);
  CL_SAFE_CALL(clRetainCommandQueue(mat.cq()));
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _host_alloc(mat_type &mat, UCL_Device &dev, const size_t n,  
                       const enum UCL_MEMOPT kind, const enum UCL_MEMOPT kind2){
  cl_mem_flags buffer_perm;
  cl_map_flags map_perm;
  if (kind==UCL_READ_ONLY) {
    #ifdef CL_VERSION_1_2
    buffer_perm=CL_MEM_HOST_READ_ONLY|CL_MEM_WRITE_ONLY|CL_MEM_ALLOC_HOST_PTR;
    #else
    buffer_perm=CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR;
    #endif
    map_perm=CL_MAP_READ;
  } else if (kind==UCL_WRITE_ONLY) {
    #ifdef CL_VERSION_1_2
    buffer_perm=CL_MEM_HOST_WRITE_ONLY|CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR;
    #else
    buffer_perm=CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR;
    #endif
    map_perm=CL_MAP_WRITE;
  } else {
    buffer_perm=CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR;
    map_perm=CL_MAP_READ | CL_MAP_WRITE;
  }

  cl_int error_flag;
  mat.cbegin()=clCreateBuffer(dev.context(),buffer_perm,n,NULL,&error_flag);
  if (error_flag != CL_SUCCESS) 
    return UCL_MEMORY_ERROR;

  *mat.host_ptr() = (typename mat_type::data_type*)
                    clEnqueueMapBuffer(dev.cq(),mat.cbegin(),CL_TRUE,
                                       map_perm,0,n,0,NULL,NULL,NULL);
  mat.cq()=dev.cq();
  CL_SAFE_CALL(clRetainCommandQueue(mat.cq()));
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _host_view(mat_type &mat, UCL_Device &dev, const size_t n) {
  cl_int error_flag;
  mat.cbegin()=clCreateBuffer(dev.context(), CL_MEM_USE_HOST_PTR,
                              n,*mat.host_ptr(),&error_flag);
  CL_CHECK_ERR(error_flag);
  CL_SAFE_CALL(clRetainCommandQueue(mat.cq()));
  return UCL_SUCCESS;
}

template <class mat_type>
inline void _host_free(mat_type &mat) {
  if (mat.cols()>0) {
    CL_DESTRUCT_CALL(clReleaseMemObject(mat.cbegin()));
    CL_DESTRUCT_CALL(clReleaseCommandQueue(mat.cq()));
  }
}

template <class mat_type>
inline int _host_resize(mat_type &mat, const size_t n) {
  cl_int error_flag;
  cl_context context;
  CL_SAFE_CALL(clGetMemObjectInfo(mat.cbegin(),CL_MEM_CONTEXT,sizeof(context),
                                  &context,NULL));
  cl_mem_flags buffer_perm;
  CL_SAFE_CALL(clGetMemObjectInfo(mat.cbegin(),CL_MEM_FLAGS,sizeof(buffer_perm),
                                  &buffer_perm,NULL));

  CL_DESTRUCT_CALL(clReleaseMemObject(mat.cbegin()));

  cl_map_flags map_perm;
  if (mat.kind()==UCL_READ_ONLY)
    map_perm=CL_MAP_READ;
  else if (mat.kind()==UCL_WRITE_ONLY)
    map_perm=CL_MAP_WRITE;
  else
    map_perm=CL_MAP_READ | CL_MAP_WRITE;

  mat.cbegin()=clCreateBuffer(context,buffer_perm,n,NULL,&error_flag);
  if (error_flag != CL_SUCCESS) 
    return UCL_MEMORY_ERROR;
  *mat.host_ptr() = (typename mat_type::data_type*)
                    clEnqueueMapBuffer(mat.cq(),mat.cbegin(),CL_TRUE,
                                       map_perm,0,n,0,NULL,NULL,NULL);
  return UCL_SUCCESS;
}

// --------------------------------------------------------------------------
// - DEVICE MEMORY ALLOCATION ROUTINES
// --------------------------------------------------------------------------

template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, copy_type &cm, const size_t n,
                         const enum UCL_MEMOPT kind) {
  cl_int error_flag;

  cl_context context;
  CL_SAFE_CALL(clGetMemObjectInfo(cm.cbegin(),CL_MEM_CONTEXT,sizeof(context),
               &context,NULL));
  cl_mem_flags flag;
  if (kind==UCL_READ_WRITE)
    flag=CL_MEM_READ_WRITE;
  else if (kind==UCL_READ_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY;
    #else
    flag=CL_MEM_READ_ONLY;
    #endif
  else if (kind==UCL_WRITE_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY;
    #else
    flag=CL_MEM_WRITE_ONLY;
    #endif
  else
    assert(0==1);
  mat.cbegin()=clCreateBuffer(context,flag,n,NULL,&error_flag);
  if (error_flag != CL_SUCCESS) 
    return UCL_MEMORY_ERROR;
  mat.cq()=cm.cq();
  CL_SAFE_CALL(clRetainCommandQueue(mat.cq()));
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _device_alloc(mat_type &mat, UCL_Device &dev, const size_t n,
                         const enum UCL_MEMOPT kind) {
  cl_int error_flag;
  cl_mem_flags flag;
  if (kind==UCL_READ_WRITE)
    flag=CL_MEM_READ_WRITE;
  else if (kind==UCL_READ_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY;
    #else
    flag=CL_MEM_READ_ONLY;
    #endif
  else if (kind==UCL_WRITE_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY;
    #else
    flag=CL_MEM_WRITE_ONLY;
    #endif
  else
    assert(0==1);
  mat.cbegin()=clCreateBuffer(dev.context(),flag,n,NULL,
                              &error_flag);
  if (error_flag != CL_SUCCESS) 
    return UCL_MEMORY_ERROR;
  mat.cq()=dev.cq();
  CL_SAFE_CALL(clRetainCommandQueue(mat.cq()));
  return UCL_SUCCESS;
}

template <class mat_type, class copy_type>
inline int _device_alloc(mat_type &mat, copy_type &cm, const size_t rows,
                         const size_t cols, size_t &pitch,
                         const enum UCL_MEMOPT kind) {
  size_t padded_cols=cols;
  if (cols%256!=0)
    padded_cols+=256-cols%256;
  pitch=padded_cols*sizeof(typename mat_type::data_type);
  return _device_alloc(mat,cm,pitch*rows,kind);
}

template <class mat_type>
inline int _device_alloc(mat_type &mat, UCL_Device &dev, const size_t rows,
                         const size_t cols, size_t &pitch,
                         const enum UCL_MEMOPT kind) {
  size_t padded_cols=cols;
  if (dev.device_type()!=UCL_CPU && cols%256!=0)
    padded_cols+=256-cols%256;
  pitch=padded_cols*sizeof(typename mat_type::data_type);
  return _device_alloc(mat,dev,pitch*rows,kind);  
}

template <class mat_type>
inline void _device_free(mat_type &mat) {
  if (mat.cols()>0) {
    CL_DESTRUCT_CALL(clReleaseMemObject(mat.cbegin()));
    CL_DESTRUCT_CALL(clReleaseCommandQueue(mat.cq()));
  }
}

template <class mat_type>
inline int _device_resize(mat_type &mat, const size_t n) {
  cl_int error_flag;

  cl_context context;
  CL_SAFE_CALL(clGetMemObjectInfo(mat.cbegin(),CL_MEM_CONTEXT,sizeof(context),
               &context,NULL));
  CL_DESTRUCT_CALL(clReleaseMemObject(mat.cbegin()));

  cl_mem_flags flag;
  if (mat.kind()==UCL_READ_WRITE)
    flag=CL_MEM_READ_WRITE;
  else if (mat.kind()==UCL_READ_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY;
    #else
    flag=CL_MEM_READ_ONLY;
    #endif
  else if (mat.kind()==UCL_WRITE_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY;
    #else
    flag=CL_MEM_WRITE_ONLY;
    #endif
  else
    assert(0==1);
  mat.cbegin()=clCreateBuffer(context,flag,n,NULL,&error_flag);
  if (error_flag != CL_SUCCESS) 
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}

template <class mat_type>
inline int _device_resize(mat_type &mat, const size_t rows,
                         const size_t cols, size_t &pitch) {
  size_t padded_cols=cols;
  if (cols%256!=0)
    padded_cols+=256-cols%256;
  pitch=padded_cols*sizeof(typename mat_type::data_type);

  cl_int error_flag;

  cl_context context;
  CL_SAFE_CALL(clGetMemObjectInfo(mat.cbegin(),CL_MEM_CONTEXT,sizeof(context),
               &context,NULL));
  CL_DESTRUCT_CALL(clReleaseMemObject(mat.cbegin()));

  cl_mem_flags flag;
  if (mat.kind()==UCL_READ_WRITE)
    flag=CL_MEM_READ_WRITE;
  else if (mat.kind()==UCL_READ_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_READ_ONLY | CL_MEM_HOST_WRITE_ONLY;
    #else
    flag=CL_MEM_READ_ONLY;
    #endif
  else if (mat.kind()==UCL_WRITE_ONLY)
    #ifdef CL_VERSION_1_2
    flag=CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY;
    #else
    flag=CL_MEM_WRITE_ONLY;
    #endif
  else
    assert(0==1);
  mat.cbegin()=clCreateBuffer(context,flag,pitch*rows,NULL,&error_flag);
  if (error_flag != CL_SUCCESS) 
    return UCL_MEMORY_ERROR;
  return UCL_SUCCESS;
}


// --------------------------------------------------------------------------
// - ZERO ROUTINES
// --------------------------------------------------------------------------
inline void _host_zero(void *ptr, const size_t n) {
  memset(ptr,0,n);
}

inline void _ocl_build(cl_program &program, cl_device_id &device,
                       const char* options = "") {
  clBuildProgram(program,1,&device,options,NULL,NULL);
    
  cl_build_status build_status;
  CL_SAFE_CALL(clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 
                                     sizeof(cl_build_status),&build_status,
                                     NULL));
  if (build_status == CL_SUCCESS)
    return;
    
  size_t ms;
  CL_SAFE_CALL(clGetProgramBuildInfo(program, device,CL_PROGRAM_BUILD_LOG, 0, 
                                     NULL, &ms));
  char build_log[ms];                                     
  CL_SAFE_CALL(clGetProgramBuildInfo(program,device,CL_PROGRAM_BUILD_LOG,ms,
                                     build_log, NULL));
    
  std::cerr << std::endl
            << "----------------------------------------------------------\n"
            << " Error compiling OpenCL Program...\n"
            << "----------------------------------------------------------\n";
  std::cerr << build_log << std::endl;
}

inline void _ocl_kernel_from_source(cl_context &context, cl_device_id &device,
                                    const char **source, const size_t lines,
                                    cl_kernel &kernel, const char *function,
                                    const char *options="") {
  cl_int error_flag;
  
  cl_program program=clCreateProgramWithSource(context,lines,source,
                                               NULL,&error_flag);
  CL_CHECK_ERR(error_flag);                                               
  _ocl_build(program,device,options);
  kernel=clCreateKernel(program,function,&error_flag);
  CL_CHECK_ERR(error_flag);                                               
}

template <class mat_type>
inline void _device_zero(mat_type &mat, const size_t n, command_queue &cq) {
  #ifdef CL_VERSION_1_2
  #ifndef __APPLE__
  #define UCL_CL_ZERO
  #endif
  #endif

  #ifdef UCL_CL_ZERO
  cl_int zeroint=0;
  CL_SAFE_CALL(clEnqueueFillBuffer(cq,mat.begin(),&zeroint,sizeof(cl_int),
                                   mat.byteoff(),n,0,NULL,NULL));

  #else
  cl_context context;
  CL_SAFE_CALL(clGetMemObjectInfo(mat.cbegin(),CL_MEM_CONTEXT,sizeof(context),
                                  &context,NULL));
  cl_device_id device;
  CL_SAFE_CALL(clGetContextInfo(context,CL_CONTEXT_DEVICES,
               sizeof(cl_device_id),&device,NULL));
  
  const char * szero[3]={
    "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n",
    "__kernel void _device_zero(__global NUMTYP *a, const int offset)",
    "  { int gid=get_global_id(0)+offset; a[gid]=(NUMTYP)0; }"
  };
  
  cl_kernel kzero;
  _ocl_kernel_from_source(context,device,szero,3,kzero,"_device_zero",
                   _UCL_DATA_ID<typename mat_type::data_type>::numtyp_flag());
  
  cl_int offset=mat.offset();
  CL_SAFE_CALL(clSetKernelArg(kzero,0,sizeof(cl_mem),(void *)&mat.begin()));
  CL_SAFE_CALL(clSetKernelArg(kzero,1,sizeof(cl_int),(void *)&offset));
  size_t kn=n/sizeof(typename mat_type::data_type);
  CL_SAFE_CALL(clEnqueueNDRangeKernel(cq,kzero,1,0,&kn,0,0,0,0));
  #endif
}

// --------------------------------------------------------------------------
// - MEMCPY ROUTINES
// --------------------------------------------------------------------------

template<int mem1, int mem2> struct _ucl_memcpy;

// Both are images
template<> struct _ucl_memcpy<2,2> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq,
                        const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
};

// Destination is texture, source on device
template<> struct _ucl_memcpy<2,0> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq,
                        const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
};

// Destination is texture, source on host
template<> struct _ucl_memcpy<2,1> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq,
                        const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
};

// Source is texture, dest on device
template<> struct _ucl_memcpy<0,2> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq,
                        const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
};

// Source is texture, dest on host
template<> struct _ucl_memcpy<1,2> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq,
                        const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    assert(0==1);
  }
};

// Neither are textures, destination on host
template <> struct _ucl_memcpy<1,0> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    if (src.cbegin()==dst.cbegin()) {
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_COPY 1S\n";
      #endif
      if (block) ucl_sync(cq);
      return;
    }
    #ifdef UCL_DBG_MEM_TRACE
    std::cerr << "UCL_COPY 1NS\n";
    #endif
    CL_SAFE_CALL(clEnqueueReadBuffer(cq,src.cbegin(),block,src_offset,n,
                                     dst.begin(),0,NULL,NULL));
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq, 
                        const cl_bool block,
                        size_t dst_offset, size_t src_offset) {
    if (src.cbegin()==dst.cbegin()) {
      if (block) ucl_sync(cq);
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_COPY 2S\n";
      #endif
      return;
    }
    #ifdef UCL_DBG_MEM_TRACE
    std::cerr << "UCL_COPY 2NS\n";
    #endif
    if (spitch==dpitch && dst.cols()==src.cols() && 
        src.cols()==cols/src.element_size())
      CL_SAFE_CALL(clEnqueueReadBuffer(cq,src.cbegin(),block,src_offset,
                                       spitch*rows,
                                       (char *)dst.begin()+dst_offset,0,NULL,
                                       NULL));
    else
      for (size_t i=0; i<rows; i++) {                       
        CL_SAFE_CALL(clEnqueueReadBuffer(cq,src.cbegin(),block,src_offset,cols,
                                         (char *)dst.begin()+dst_offset,0,NULL,
                                         NULL));
        src_offset+=spitch;
        dst_offset+=dpitch;
      }                                       
  }
};

// Neither are textures, source on host
template <> struct _ucl_memcpy<0,1> {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    if (src.cbegin()==dst.cbegin()) {
      if (block) ucl_sync(cq);
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_COPY 3S\n";
      #endif
      return;                        
    }
    #ifdef UCL_DBG_MEM_TRACE
    std::cerr << "UCL_COPY 3NS\n";
    #endif
    CL_SAFE_CALL(clEnqueueWriteBuffer(cq,dst.cbegin(),block,dst_offset,n,
                                      src.begin(),0,NULL,NULL));
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq, 
                        const cl_bool block,
                        size_t dst_offset, size_t src_offset) {
    if (src.cbegin()==dst.cbegin()) {
      if (block) ucl_sync(cq);
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_COPY 4S\n";
      #endif
      return;                        
    }
    #ifdef UCL_DBG_MEM_TRACE
    std::cerr << "UCL_COPY 4NS\n";
    #endif
    if (spitch==dpitch && dst.cols()==src.cols() && 
        src.cols()==cols/src.element_size())
      CL_SAFE_CALL(clEnqueueWriteBuffer(cq,dst.cbegin(),block,dst_offset,
                                        spitch*rows,
                                        (char *)src.begin()+src_offset,0,NULL,
                                        NULL));
    else
      for (size_t i=0; i<rows; i++) {
        CL_SAFE_CALL(clEnqueueWriteBuffer(cq,dst.cbegin(),block,dst_offset,cols,
                                          (char *)src.begin()+src_offset,0,NULL,
                                          NULL));
        src_offset+=spitch;
        dst_offset+=dpitch;
      }                                       
  }
};

// Neither are textures, both on device
template <int mem1, int mem2> struct _ucl_memcpy {
  template <class p1, class p2>
  static inline void mc(p1 &dst, const p2 &src, const size_t n,
                        cl_command_queue &cq, const cl_bool block,
                        const size_t dst_offset, const size_t src_offset) {
    if (src.cbegin()!=dst.cbegin() || src_offset!=dst_offset) {
      CL_SAFE_CALL(clEnqueueCopyBuffer(cq,src.cbegin(),dst.cbegin(),src_offset,
                                       dst_offset,n,0,NULL,NULL));
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_COPY 6NS\n";
      #endif
    }
    #ifdef UCL_DBG_MEM_TRACE
    else std::cerr << "UCL_COPY 6S\n";
    #endif
    
    if (block==CL_TRUE) ucl_sync(cq);
  }
  template <class p1, class p2>
  static inline void mc(p1 &dst, const size_t dpitch, const p2 &src, 
                        const size_t spitch, const size_t cols,
                        const size_t rows, cl_command_queue &cq,
                        const cl_bool block,
                        size_t dst_offset, size_t src_offset) {
    if (src.cbegin()!=dst.cbegin() || src_offset!=dst_offset) {                        
      #ifdef UCL_DBG_MEM_TRACE
      std::cerr << "UCL_COPY 7NS\n";
      #endif
      if (spitch==dpitch && dst.cols()==src.cols() && 
          src.cols()==cols/src.element_size())
        CL_SAFE_CALL(clEnqueueCopyBuffer(cq,src.cbegin(),dst.cbegin(),src_offset,
                                         dst_offset,spitch*rows,0,NULL,NULL));
        
      else
        for (size_t i=0; i<rows; i++) {                       
          CL_SAFE_CALL(clEnqueueCopyBuffer(cq,src.cbegin(),dst.cbegin(),
                                           src_offset,dst_offset,cols,0,
                                           NULL,NULL));
          src_offset+=spitch;
          dst_offset+=dpitch;
        }                                       
    }                                 
    #ifdef UCL_DBG_MEM_TRACE
    else std::cerr << "UCL_COPY 7S\n";
    #endif

    if (block==CL_TRUE) ucl_sync(cq);
  }
};

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const mat2 &src, const size_t n) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,src,n,dst.cq(),CL_TRUE,
                                                 dst.byteoff(),src.byteoff());
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const mat2 &src, const size_t n,
                       cl_command_queue &cq) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,src,n,cq,CL_FALSE,
                                                 dst.byteoff(),src.byteoff());
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const size_t dpitch, const mat2 &src, 
                       const size_t spitch, const size_t cols, 
                       const size_t rows) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,dpitch,src,spitch,cols,
                                                 rows,dst.cq(),CL_TRUE,
                                                 dst.byteoff(),src.byteoff());
}

template<class mat1, class mat2>
inline void ucl_mv_cpy(mat1 &dst, const size_t dpitch, const mat2 &src, 
                           const size_t spitch, const size_t cols, 
                           const size_t rows,cl_command_queue &cq) {
  _ucl_memcpy<mat1::MEM_TYPE,mat2::MEM_TYPE>::mc(dst,dpitch,src,spitch,cols,
                                                 rows,cq,CL_FALSE,
                                                 dst.byteoff(),src.byteoff());
}

} // namespace ucl_cudart 

#endif

