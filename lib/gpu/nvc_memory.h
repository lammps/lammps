/***************************************************************************
                                nvc_memory.h
                             -------------------
                               W. Michael Brown

  Routines for memory management on CUDA devices

 __________________________________________________________________________
    This file is part of the NVC Library
 __________________________________________________________________________

    begin                : Thu Jun 25 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#ifndef NVC_MEMORY_H
#define NVC_MEMORY_H

#include <iostream>
#include "nvc_macros.h"

#define NVC_HostT NVC_Host<numtyp>
#define NVC_HostD NVC_Host<double>
#define NVC_HostS NVC_Host<float>
#define NVC_HostI NVC_Host<int>

#define NVC_VecT NVC_Vec<numtyp>
#define NVC_VecD NVC_Vec<double>
#define NVC_VecS NVC_Vec<float>
#define NVC_VecI NVC_Vec<int>
#define NVC_VecI2 NVC_Vec<int2>
#define NVC_VecU2 NVC_Vec<uint2>

#define NVC_MatT NVC_Mat<numtyp>
#define NVC_MatD NVC_Mat<double>
#define NVC_MatS NVC_Mat<float>
#define NVC_MatI NVC_Mat<int>

#define NVC_ConstMatT NVC_ConstMat<numtyp>
#define NVC_ConstMatD NVC_ConstMat<double>
#define NVC_ConstMatS NVC_ConstMat<float>
#define NVC_ConstMatI NVC_ConstMat<int>
#define NVC_ConstMatD2 NVC_ConstMat<double2>

namespace NVC {

// Get a channel for float array
template <class numtyp>
inline void cuda_gb_get_channel(cudaChannelFormatDesc &channel) {
  channel = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
}

// Get a channel for float2 array
template <>
inline void cuda_gb_get_channel<float2>(cudaChannelFormatDesc &channel) {
  channel = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
}

// Get a channel for double array
template <>
inline void cuda_gb_get_channel<double>(cudaChannelFormatDesc &channel) {
  channel = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindSigned);
}

// Get a channel for double array
template <>
inline void cuda_gb_get_channel<double2>(cudaChannelFormatDesc &channel) {
  channel = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindSigned);
}

// Get a channel for int array
template <>
inline void cuda_gb_get_channel<int>(cudaChannelFormatDesc &channel) {
  channel = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindSigned);
}

}

/// Page-locked Row Vector on Host
template <class numtyp>
class NVC_Host {
 public:
  NVC_Host() { _cols=0; }
  ~NVC_Host() { if (_cols>0) CUDA_SAFE_CALL(cudaFreeHost(_array)); }
  
  // Allocate page-locked memory with fast write/slow read on host
  inline void safe_alloc_w(const size_t cols) {
    _cols=cols;
    _row_bytes=cols*sizeof(numtyp);
    CUDA_SAFE_CALL(cudaHostAlloc((void **)&_array,_row_bytes,
                   cudaHostAllocWriteCombined));
    _end=_array+cols;
  }

  // Allocate page-locked memory with fast read/write on host
  inline void safe_alloc_rw(const size_t cols) {
    _cols=cols;
    _row_bytes=cols*sizeof(numtyp);
    CUDA_SAFE_CALL(cudaMallocHost((void **)&_array,_row_bytes));
    _end=_array+cols;
  }
  
  /// Free any memory associated with device
  inline void clear() 
    { if (_cols>0) { _cols=0; CUDA_SAFE_CALL(cudaFreeHost(_array)); } }

  /// Set each element to zero
  inline void zero() { memset(_array,0,row_bytes()); }
  
  /// Set first n elements to zero
  inline void zero(const int n) { memset(_array,0,n*sizeof(numtyp)); }

  inline numtyp * begin() { return _array; }
  inline const numtyp * begin() const { return _array; }
  inline numtyp * end() { return _end; }
  inline const numtyp * end() const { return _end; }

  inline size_t numel() const { return _cols; }
  inline size_t rows() const { return 1; }
  inline size_t cols() const { return _cols; }
  inline size_t row_size() const { return _cols; }
  inline size_t row_bytes() const { return _row_bytes; }
    
  inline numtyp & operator[](const int i) { return _array[i]; }
  inline const numtyp & operator[](const int i) const { return _array[i]; }
  
  /// Copy from device (numel is not bytes)
  inline void copy_from_device(const numtyp *device_p, size_t numel) {
    CUDA_SAFE_CALL(cudaMemcpy(_array,device_p,numel*sizeof(numtyp),
                              cudaMemcpyDeviceToHost));
  }
  
  /// Copy to device (numel is not bytes)
  inline void copy_to_device(numtyp *device_p, size_t numel) {
    CUDA_SAFE_CALL(cudaMemcpy(device_p,_array,numel*sizeof(numtyp),
                              cudaMemcpyHostToDevice));
  }
  
  /// Copy to 2D matrix on device (numel is not bytes)
  inline void copy_to_2Ddevice(numtyp *device_p, const size_t dev_row_size,
                               const size_t rows, const size_t cols) {
    CUDA_SAFE_CALL(cudaMemcpy2D(device_p,dev_row_size*sizeof(numtyp),
                                _array,cols*sizeof(numtyp),
                                cols*sizeof(numtyp),rows,
                                cudaMemcpyHostToDevice));
  }

  /// Asynchronous copy from device (numel is not bytes)
  inline void copy_from_device(const numtyp *device_p, size_t numel,
                               cudaStream_t &stream) {
    CUDA_SAFE_CALL_NO_SYNC(cudaMemcpyAsync(_array,device_p,numel*sizeof(numtyp),
                                           cudaMemcpyDeviceToHost,stream));
  }
  
  /// Asynchronous copy to device (numel is not bytes)
  inline void copy_to_device(numtyp *device_p, size_t numel,
                             cudaStream_t &stream) {
    CUDA_SAFE_CALL_NO_SYNC(cudaMemcpyAsync(device_p,_array,numel*sizeof(numtyp),
                                           cudaMemcpyHostToDevice,stream));
  }
  
  /// Asynchronous copy to 2D matrix on device (numel is not bytes)
  inline void copy_to_2Ddevice(numtyp *device_p, const size_t dev_row_size,
                               const size_t rows, const size_t cols,
                               cudaStream_t &stream) {
    CUDA_SAFE_CALL_NO_SYNC(cudaMemcpy2DAsync(device_p,
                                             dev_row_size*sizeof(numtyp),
                                             _array,cols*sizeof(numtyp),
                                             cols*sizeof(numtyp),rows,
                                             cudaMemcpyHostToDevice,stream));
  }
  
 private:
  numtyp *_array, *_end;
  size_t _row_bytes, _row_size, _rows, _cols;
};

/// Row vector on device 
template <class numtyp>
class NVC_Vec {
 public:
  NVC_Vec() { _cols=0; }
  ~NVC_Vec() { if (_cols>0) CUDA_SAFE_CALL(cudaFree(_array)); }
  
  // Row vector on device
  inline void safe_alloc(const size_t cols) {
    _cols=cols;
    _row_bytes=cols*sizeof(numtyp);
    CUDA_SAFE_CALL(cudaMalloc((void **)&_array,_row_bytes));
    _end=_array+cols;
  }
  
  // Row vector on device (allocate and assign texture and bind)
  inline void safe_alloc(const size_t cols, textureReference *t) 
    { safe_alloc(cols); assign_texture(t); bind(); }

  /// Free any memory associated with device
  inline void clear() 
    { if (_cols>0) { _cols=0; CUDA_SAFE_CALL(cudaFree(_array)); } }

  /// Set each element to zero
  inline void zero() { CUDA_SAFE_CALL(cudaMemset(_array,0,row_bytes())); }

  inline numtyp * begin() { return _array; }
  inline const numtyp * begin() const { return _array; }
  inline numtyp * end() { return _end; }
  inline const numtyp * end() const { return _end; }

  inline size_t numel() const { return _cols; }
  inline size_t rows() const { return 1; }
  inline size_t cols() const { return _cols; }
  inline size_t row_size() const { return _cols; }
  inline size_t row_bytes() const { return _row_bytes; }
    
  /// Copy from host
  inline void copy_from_host(const numtyp *host_p)
    { CUDA_SAFE_CALL(cudaMemcpy(_array,host_p,row_bytes(), 
                                cudaMemcpyHostToDevice)); }

  /// Asynchronous copy from host
  inline void copy_from_host(const numtyp *host_p, cudaStream_t &stream)
    { CUDA_SAFE_CALL_NO_SYNC(cudaMemcpyAsync(_array,host_p,row_bytes(), 
                             cudaMemcpyHostToDevice, stream)); }

  /// Copy to host
  inline void copy_to_host(numtyp *host_p)
    { CUDA_SAFE_CALL(cudaMemcpy(host_p,_array,row_bytes(),
                                cudaMemcpyDeviceToHost)); }

  /// Copy n elements to host
  inline void copy_to_host(numtyp *host_p, const int n)
    { CUDA_SAFE_CALL(cudaMemcpy(host_p,_array,n*sizeof(numtyp),
                                cudaMemcpyDeviceToHost)); }

  /// Cast and then copy to device
  template <class numtyp2>
  inline void cast_copy(const numtyp2 *buffer, NVC_HostT &host_write) {
    for (int i=0; i<numel(); i++)
      host_write[i]=static_cast<numtyp>(buffer[i]);
    copy_from_host(host_write.begin());
  }
  
  /// Assign a texture to matrix
  inline void assign_texture(textureReference *t) { _tex_ptr=t; }  

  /// Bind to texture
  inline void bind() {
    NVC::cuda_gb_get_channel<numtyp>(_channel);
    (*_tex_ptr).addressMode[0] = cudaAddressModeClamp;
    (*_tex_ptr).addressMode[1] = cudaAddressModeClamp;
    (*_tex_ptr).filterMode = cudaFilterModePoint;
    (*_tex_ptr).normalized = false;
    CUDA_SAFE_CALL(cudaBindTexture(NULL,_tex_ptr,_array,&_channel));
  }

  /// Unbind texture
  inline void unbind() { CUDA_SAFE_CALL(cudaUnbindTexture(_tex_ptr)); }

  /// Output the vector (debugging)
  inline void print(std::ostream &out) { print (out, numel()); }    

  // Output first n elements of vector
  inline void print(std::ostream &out, const int n) {
    numtyp *t=new numtyp[n];
    copy_to_host(t,n);
    for (int i=0; i<n; i++)
      out << t[i] << " ";
    delete []t;
  }
    
 private:
  numtyp *_array, *_end;
  size_t _row_bytes, _row_size, _rows, _cols;
  cudaChannelFormatDesc _channel;
  textureReference *_tex_ptr;
};

/// 2D Matrix on device (can have extra column storage to get correct alignment)
template <class numtyp>
class NVC_Mat {
 public:
  NVC_Mat() { _rows=0; }
  ~NVC_Mat() { if (_rows>0) CUDA_SAFE_CALL(cudaFree(_array)); }
  
  // Row major matrix on device
  // - Coalesced access using adjacent cols on same row
  // - NVC_Mat(row,col) given by array[row*row_size()+col]
  inline void safe_alloc(const size_t rows, const size_t cols) {
    _rows=rows;
    _cols=cols;
    CUDA_SAFE_CALL(cudaMallocPitch((void **)&_array,&_pitch, 
                   cols*sizeof(numtyp),rows));
   _row_size=_pitch/sizeof(numtyp);                
    _end=_array+_row_size*cols;
  }
  
  /// Free any memory associated with device
  inline void clear() 
    { if (_rows>0) { _rows=0; CUDA_SAFE_CALL(cudaFree(_array)); } }

  /// Set each element to zero
  inline void zero() { CUDA_SAFE_CALL(cudaMemset(_array,0, _pitch*_rows)); }

  inline numtyp * begin() { return _array; }
  inline const numtyp * begin() const { return _array; }
  inline numtyp * end() { return _end; }
  inline const numtyp * end() const { return _end; }


  inline size_t numel() const { return _cols*_rows; }
  inline size_t rows() const { return _rows; }
  inline size_t cols() const { return _cols; }
  inline size_t row_size() const { return _row_size; }
  inline size_t row_bytes() const { return _pitch; }
    
  /// Copy from host (elements not bytes)
  inline void copy_from_host(const numtyp *host_p, const size_t numel)
    { CUDA_SAFE_CALL(cudaMemcpy(_array,host_p,numel*sizeof(numtyp), 
                     cudaMemcpyHostToDevice)); }

  /// Asynchronous copy from host (elements not bytes)
  inline void copy_from_host(const numtyp *host_p, const size_t numel,
                             cudaStream_t &stream)
    { CUDA_SAFE_CALL_NO_SYNC(cudaMemcpyAsync(_array,host_p,numel*sizeof(numtyp), 
                             cudaMemcpyHostToDevice, stream)); }

  /// Asynchronous Copy from Host
  /** \note Used when the number of columns/rows allocated on host smaller than
    *       on device **/
  inline void copy_2Dfrom_host(const numtyp *host_p, const size_t rows,
                               const size_t cols, cudaStream_t &stream) {
    CUDA_SAFE_CALL_NO_SYNC(cudaMemcpy2DAsync(_array, _pitch, host_p,
                           cols*sizeof(numtyp), cols*sizeof(numtyp), rows,
                           cudaMemcpyHostToDevice,stream));
  }

 private:
  numtyp *_array, *_end;
  size_t _pitch, _row_size, _rows, _cols;
};

/// Const 2D Matrix on device (requires texture binding)
template <class numtyp>
class NVC_ConstMat {
 public:
  NVC_ConstMat() { _rows=0; }
  ~NVC_ConstMat() { if (_rows>0) CUDA_SAFE_CALL(cudaFreeArray(_array)); }

  /// Assign a texture to matrix
  inline void assign_texture(textureReference *t) { _tex_ptr=t; }  
      
  /// Row major matrix on device
  inline void safe_alloc(const size_t rows, const size_t cols) {
    _rows=rows;
    _cols=cols;

    NVC::cuda_gb_get_channel<numtyp>(_channel);
    CUDA_SAFE_CALL(cudaMallocArray(&_array, &_channel, cols, rows));
  }
  
  /// Row major matrix on device (Allocate and bind texture)
  inline void safe_alloc(const size_t rows, const size_t cols, 
                         textureReference *t) 
    { safe_alloc(rows,cols); assign_texture(t); bind(); }

  /// Bind to texture
  inline void bind() {
    (*_tex_ptr).addressMode[0] = cudaAddressModeClamp;
    (*_tex_ptr).addressMode[1] = cudaAddressModeClamp;
    (*_tex_ptr).filterMode = cudaFilterModePoint;
    (*_tex_ptr).normalized = false;
    CUDA_SAFE_CALL(cudaBindTextureToArray(_tex_ptr,_array,&_channel));
  }
  
  /// Unbind texture
  inline void unbind() { CUDA_SAFE_CALL(cudaUnbindTexture(_tex_ptr)); }
  
  /// Free any memory associated with device and unbind
  inline void clear() {
    if (_rows>0) { 
      _rows=0; 
      CUDA_SAFE_CALL(cudaUnbindTexture(_tex_ptr)); 
      CUDA_SAFE_CALL(cudaFreeArray(_array)); 
    } 
  }

  inline size_t numel() const { return _cols*_rows; }
  inline size_t rows() const { return _rows; }
  inline size_t cols() const { return _cols; }
  inline size_t row_size() const { return _cols; }
  inline size_t row_bytes() const { return _cols*sizeof(numtyp); }

  /// Copy from Host
  inline void copy_from_host(const numtyp *host_p) {
    CUDA_SAFE_CALL(cudaMemcpyToArray(_array, 0, 0, host_p,
                                     numel()*sizeof(numtyp),
                                     cudaMemcpyHostToDevice));
  }

  /// Copy from Host
  /** \note Used when the number of columns/rows allocated on host smaller than
    *       on device **/
  inline void copy_2Dfrom_host(const numtyp *host_p, const size_t rows,
                               const size_t cols) {
    CUDA_SAFE_CALL(cudaMemcpy2DToArray(_array, 0, 0, host_p, 
                         cols*sizeof(numtyp), cols*sizeof(numtyp), rows,
                         cudaMemcpyHostToDevice));
  }                               

  /// Asynchronous Copy from Host
  inline void copy_from_host(const numtyp *host_p, cudaStream_t &stream) {
    CUDA_SAFE_CALL_NO_SYNC(cudaMemcpyToArrayAsync(_array, 0, 0, host_p, 
                                                  numel()*sizeof(numtyp),
                                                  cudaMemcpyHostToDevice,
                                                  stream));
  }
      
  /// Asynchronous Copy from Host
  /** \note Used when the number of columns/rows allocated on host smaller than
    *       on device **/
  inline void copy_2Dfrom_host(const numtyp *host_p, const size_t rows,
                               const size_t cols, cudaStream_t &stream) {
    CUDA_SAFE_CALL_NO_SYNC(cudaMemcpy2DToArrayAsync(_array, 0, 0, host_p, 
                           cols*sizeof(numtyp), cols*sizeof(numtyp), rows,
                           cudaMemcpyHostToDevice,stream));
  }                               

  /// Cast buffer to numtyp in host_write and copy to array
  template <class numtyp2>
  inline void cast_copy(const numtyp2 *buffer, NVC_HostT &host_write) {
    int n=numel();                        
    for (int i=0; i<n; i++) {
      host_write[i]=static_cast<numtyp>(*buffer); buffer++;
    }
    copy_from_host(host_write.begin());
  }

  /// Cast buffer to numtyp in host_write and copy to array
  /** \note Used when the number of columns/rows allocated on host smaller than
    *       on device **/
  template <class numtyp2>
  inline void cast_copy2D(const numtyp2 *buffer, NVC_HostT &host_write,
                          const size_t rows, const size_t cols) {
    int n=rows*cols;                        
    for (int i=0; i<n; i++) {
      host_write[i]=static_cast<numtyp>(*buffer); buffer++;
    }
    copy_2Dfrom_host(host_write.begin(),rows,cols);
  }

  /// Cast buffer to numtyp in host_write and copy to array asynchronously
  template <class numtyp2>
  inline void cast_copy(const numtyp2 *buffer, NVC_HostT &host_write,
                        cudaStream_t &stream) {
    int n=numel();                        
    for (int i=0; i<n; i++) {
      host_write[i]=static_cast<numtyp>(*buffer); buffer++;
    }
    copy_from_host(host_write.begin(),stream);
  }

 private:
  size_t _rows, _cols;
  cudaArray *_array;
  cudaChannelFormatDesc _channel;
  textureReference *_tex_ptr;
};

#endif
