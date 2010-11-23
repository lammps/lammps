/***************************************************************************
                                 ucl_copy.h
                             -------------------
                               W. Michael Brown

  Routines for copying matrix/vector data onto and off coprocessor device

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Mon Jan 4 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */
   
/***************************************************************************
   The ucl_copy and ucl_cast_copy routines provide a general prototype for
   copying data between host and device memory (including texture memory)
   for the matrix and vector types in nvc_memory.
   
   For host/host and host/device transfers, typecasting is performed 
   automatically as necessary. 
   
   The routines are written so that all branches can be removed by the 
   compiler during template instantiation.
   
   The routines currently assume row-major ordering for all types.
   
   For asynchronous copy in the default command queue, async is boolean true;
   For asynchronous copy in a specified command queue, async is command queue
   Otherwise, set async to boolean false;
   
   When performing frequent data copies that require casting, it is more
   efficient to allocate a casting buffer once and then pass that buffer
   to the copy routine. This can be accomplished with the ucl_cast_copy
   routines.
   
   Examples 
      (x's represent alignment padding - to maintain alignment)
      (o's represent a larger matrix in memory)
      (vectors represented as single row)
   ----------------------------------------------------------------
       dst           src            command
   ----------------------------------------------------------------
    0 1 2 3 4 <-- 0 1 2 3 4          ucl_copy(dst,src,async)
    
    0 1 2 3   <-- 0 1 2 3 4          ucl_copy(dst,src,4,async)
    
    0 1 2     <-- 0 1 2 3 4 5        ucl_copy(dst,src,async)
    3 4 5 
   
    0 1 2 3 4 5 <-- 0 1 2            ucl_copy(dst,src,async)
                    3 4 5
                    
    0 1 2      <--  0 1 2            ucl_copy(dst,src,async)
    3 4 5           3 4 5
    
    0 1 2      <--  0 1 2            ucl_copy(dst,src,6,async)
    3 4 5           3 4 5
                    5 6 7

    0 1 2      <--  0  1  2  3       ucl_copy(dst,src,2,3,async)
    4 5 6           4  5  6  7
                    8  9  10 11
    
    0 1 2 x x  <--  0 1 2            ucl_copy(dst,src,async)
    3 4 5 x x       3 4 5
    
    0 1 2      <--  0 1 2 x x        ucl_copy(dst,src,async)
    3 4 5           3 4 5 x x
    
    0 1 2 o o  <--  0 1 2            ucl_copy(dst,src,2,3,async)
    3 4 5 o o       3 4 5
    o o o o o       

    0 1 2 o o  <--  0 1 2 3 4 5      ucl_copy(dst,src,2,3,async)
    3 4 5 o o       
    o o o o o       

    0 1 o o o  <--  0 1 2 3 4 5      ucl_copy(dst,src,2,2,async)
    2 3 o o o       
    o o o o o       

    0 1 2 o o  <--  0  1  2  3  4    ucl_copy(dst,src,2,3,async)
    5 6 7 o o       5  6  7  8  9
    o o o o o       10 11 12 13 14
    
    0 1 2 5 6 7  <--  0  1  2  3  4  ucl_copy(dst,src,2,3,async)
                      5  6  7  8  9
                      10 11 12 13 14
    
 ***************************************************************************/

// Only allow this file to be included by nvc_memory.h and ocl_memory.h
#ifdef UCL_COPY_ALLOW

// --------------------------------------------------------------------------
// - HOST-HOST COPY ROUTINES
// --------------------------------------------------------------------------

// Have to use specialization because some types don't have operator[]
template <int host_t1, int host_t2> struct _host_host_copy;

// Both on host
template <> struct _host_host_copy<1,1> {
  template <class mat1, class mat2>
  static inline void hhc(mat1 &dst, const mat2 &src, const size_t numel) {
    #ifdef UCL_DEBUG
    assert(mat1::PADDED==0 && mat2::PADDED==0);
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    #endif
    if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE && mat1::DATA_TYPE!=0)
      memcpy(dst.begin(),src.begin(),numel*sizeof(typename mat1::data_type));
    else
      for (size_t i=0; i<numel; i++)
        dst[i]=static_cast<typename mat1::data_type>(src[i]);
  }
  template <class mat1, class mat2>
  static inline void hhc(mat1 &dst, const mat2 &src, const size_t rows,
                         const size_t cols) {
    #ifdef UCL_DEBUG
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    #endif
    size_t dst_row_size, src_row_size;
    if (mat1::VECTOR)
      dst_row_size=cols;
    else
      dst_row_size=dst.row_size();
    if (mat2::VECTOR)
      src_row_size=cols;
    else
      src_row_size=src.row_size();
    if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE && mat1::DATA_TYPE!=0)
      for (size_t i=0; i<rows; i++)
        memcpy(dst.begin()+i*dst_row_size,src.begin()+i*src_row_size,
               cols*sizeof(typename mat1::data_type));
    else
      for (size_t j=0; j<rows; j++) {
        int dst_i=j*dst_row_size;
        int d_end=dst_i+cols;
        int src_i=j*src_row_size;
        for (; dst_i<d_end; dst_i++) {
          dst[dst_i]=static_cast<typename mat1::data_type>(src[src_i]);
          src_i++;
        }
      }
  }
};

// Should never be here
template <int host_t1, int host_t2> struct _host_host_copy {
  template <class mat1, class mat2>
  static inline void hhc(mat1 &dst, const mat2 &src, const size_t numel) {
    assert(0==1);
  }
  template <class mat1, class mat2>
  static inline void hhc(mat1 &dst, const mat2 &src, const size_t rows,
                         const size_t cols) {
    assert(0==1);
  }                         
};

// --------------------------------------------------------------------------
// - TEMPLATE HELPER FUNCTIONS FOR SPECIALIZED CASTING
// --------------------------------------------------------------------------

// Helper functions for ucl_cast_copy
template <int host_type1, int host_type2> struct _ucl_cast_copy;

// Destination is on host
template <int host_type2> struct _ucl_cast_copy<1,host_type2> {
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer) {
    ucl_mv_cpy(cast_buffer,src,numel*sizeof(typename mat2::data_type));
    for (size_t i=0; i<numel; i++)
      dst[i]=static_cast<typename mat1::data_type>(cast_buffer[i]);
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer,command_queue &cq) {
    ucl_mv_cpy(cast_buffer,src,numel*sizeof(typename mat2::data_type),cq);
    cast_buffer.sync();
    for (size_t i=0; i<numel; i++)
      dst[i]=static_cast<typename mat1::data_type>(cast_buffer[i]);
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer) {
    // Asynchronous currently pointless here 
    #ifdef UCL_DEBUG
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    assert(dst.numel()>=rows*cols && cast_buffer.numel()>=rows*cols);
    if (mat1::VECTOR==0) assert(dst.rows()>=rows && dst.cols()>=cols);    
    if (mat2::VECTOR==0) assert(src.rows()>=rows && src.cols()>=cols);    
    #endif    
    if (mat1::VECTOR) {
      ucl_mv_cpy(cast_buffer,cols*sizeof(typename mat2::data_type),src,
                 src.row_bytes(),cols*sizeof(typename mat2::data_type),rows);
      for (size_t i=0; i<rows*cols; i++)
        dst[i]=static_cast<typename mat1::data_type>(cast_buffer[i]);
    } else {
      if (mat2::VECTOR) 
        ucl_mv_cpy(cast_buffer,cols*sizeof(typename mat2::data_type),src,
                   cols*sizeof(typename mat2::data_type),
                   cols*sizeof(typename mat2::data_type),rows);
      else
        ucl_mv_cpy(cast_buffer,cols*sizeof(typename mat2::data_type),src,
                   src.row_bytes(),cols*sizeof(typename mat2::data_type),
                   rows);
      int dst_i=0;
      int buff_i=0;
      for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
          dst[dst_i]=static_cast<typename mat1::data_type>(cast_buffer[buff_i]);
          buff_i++;
          dst_i++;
        }
        dst_i+=dst.cols()-cols;
      }
    }
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer, 
                        command_queue &cq) {
    // Asynchronous currently pointless here 
    #ifdef UCL_DEBUG
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    assert(dst.numel()>=rows*cols && cast_buffer.numel()>=rows*cols);
    if (mat1::VECTOR==0) assert(dst.rows()>=rows && dst.cols()>=cols);    
    if (mat2::VECTOR==0) assert(src.rows()>=rows && src.cols()>=cols);    
    #endif    
    if (mat1::VECTOR) {
      ucl_mv_cpy(cast_buffer,cols*sizeof(typename mat2::data_type),src,
                 src.row_bytes(),cols*sizeof(typename mat2::data_type),rows,cq);
      cast_buffer.sync();           
      for (size_t i=0; i<rows*cols; i++)
        dst[i]=static_cast<typename mat1::data_type>(cast_buffer[i]);
    } else {
      if (mat2::VECTOR) 
        ucl_mv_cpy(cast_buffer,cols*sizeof(typename mat2::data_type),src,
                   cols*sizeof(typename mat2::data_type),
                   cols*sizeof(typename mat2::data_type),rows,cq);
      else
        ucl_mv_cpy(cast_buffer,cols*sizeof(typename mat2::data_type),src,
                   src.row_bytes(),cols*sizeof(typename mat2::data_type),
                   rows,cq);
      cast_buffer.sync();
      int dst_i=0;
      int buff_i=0;
      for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
          dst[dst_i]=static_cast<typename mat1::data_type>(cast_buffer[buff_i]);
          buff_i++;
          dst_i++;
        }
        dst_i+=dst.cols()-cols;
      }
    }
  }
};

// Source is on host
template <int host_type1> struct _ucl_cast_copy<host_type1,1> {
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer) {
    for (size_t i=0; i<numel; i++)
      cast_buffer[i]=static_cast<typename mat3::data_type>(src[i]);
    ucl_mv_cpy(dst,cast_buffer,numel*sizeof(typename mat1::data_type));
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer, command_queue &cq) {
    for (size_t i=0; i<numel; i++)
      cast_buffer[i]=static_cast<typename mat3::data_type>(src[i]);
    ucl_mv_cpy(dst,cast_buffer,numel*sizeof(typename mat1::data_type),cq);
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer) {
    #ifdef UCL_DEBUG
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    assert(src.numel()>=rows*cols && cast_buffer.numel()>=rows*cols);
    if (mat1::VECTOR==0) assert(dst.rows()>=rows && dst.cols()>=cols);
    if (mat2::VECTOR==0) assert(src.rows()>=rows && src.cols()>=cols);
    #endif
    if (mat2::VECTOR) {
      for (size_t i=0; i<rows*cols; i++)
        cast_buffer[i]=static_cast<typename mat3::data_type>(src[i]);
      ucl_mv_cpy(dst,dst.row_bytes(),cast_buffer,
                 cols*sizeof(typename mat1::data_type),
                 cols*sizeof(typename mat1::data_type),rows);
    } else if (mat1::VECTOR) {
      int src_i=0;
      int buf_i=0;
      for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
          cast_buffer[buf_i]=static_cast<typename mat3::data_type>(src[src_i]);
          buf_i++;
          src_i++;
        }
        src_i+=src.cols()-cols;
      }
      ucl_mv_cpy(dst,cast_buffer,cols*sizeof(typename mat1::data_type)*rows);
    } else {
      int src_i=0;
      int buf_i=0;
      for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
          cast_buffer[buf_i]=static_cast<typename mat3::data_type>(src[src_i]);
          buf_i++;
          src_i++;
        }
        src_i+=src.cols()-cols;
      }
      ucl_mv_cpy(dst,dst.row_bytes(),cast_buffer,
                 cols*sizeof(typename mat1::data_type),
                 cols*sizeof(typename mat1::data_type),rows);
    }
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer,
                        command_queue &cq) {
    #ifdef UCL_DEBUG
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    assert(src.numel()>=rows*cols && cast_buffer.numel()>=rows*cols);
    if (mat1::VECTOR==0) assert(dst.rows()>=rows && dst.cols()>=cols);    
    if (mat2::VECTOR==0) assert(src.rows()>=rows && src.cols()>=cols);    
    #endif
    if (mat2::VECTOR) {
      for (size_t i=0; i<rows*cols; i++)
        cast_buffer[i]=static_cast<typename mat3::data_type>(src[i]);
      ucl_mv_cpy(dst,dst.row_bytes(),
                 cast_buffer,cols*sizeof(typename mat1::data_type),
                 cols*sizeof(typename mat1::data_type),rows,cq);
    } else if (mat1::VECTOR) {
      int src_i=0;
      int buf_i=0;
      for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
          cast_buffer[buf_i]=static_cast<typename mat3::data_type>(src[src_i]);
          buf_i++;
          src_i++;
        }
        src_i+=src.cols()-cols;
      }
      ucl_mv_cpy(dst,cast_buffer,cols*sizeof(typename mat1::data_type)*rows,cq);
    } else {
      int src_i=0;
      int buf_i=0;
      for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
          cast_buffer[buf_i]=static_cast<typename mat3::data_type>(src[src_i]);
          buf_i++;
          src_i++;
        }
        src_i+=src.cols()-cols;
      }
      ucl_mv_cpy(dst,dst.row_bytes(),cast_buffer,
                 cols*sizeof(typename mat1::data_type),
                 cols*sizeof(typename mat1::data_type),rows,cq);
    }
  }
};

// Neither on host or both on host
template <> struct _ucl_cast_copy<1,1> {
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer, command_queue &cq) {
    assert(0==1);                        
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer) {
    assert(0==1);                        
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer) {
    assert(0==1);                        
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer,
                        command_queue &cq) {
    assert(0==1);                        
  }
};

// Neither on host or both on host
template <> struct _ucl_cast_copy<0,0> {
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer, command_queue &cq) {
    assert(0==1);                        
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t numel,
                        mat3 &cast_buffer) {
    assert(0==1);                        
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer) {
    assert(0==1);                        
  }
  template <class mat1, class mat2, class mat3>
  static inline void cc(mat1 &dst, const mat2 &src, const size_t rows,
                        const size_t cols, mat3 &cast_buffer,
                        command_queue &cq) {
    assert(0==1);                        
  }
};

// --------------------------------------------------------------------------
// - 1D COPY - SPECIFIED NUMBER OF BYTES
// --------------------------------------------------------------------------

/// Asynchronous copy of matrix/vector with cast (Device/Host transfer)
/** \param numel Number of elements (not bytes) to copy
  * \param cast_buffer Buffer on host with enough storage for casting
  * - If the data types for the two matrices are same, no cast performed
  * - Padding for 2D matrices is not considered in this routine. 
  * - Currently does not handle textures **/
template <class mat1, class mat2, class mat3>
inline void ucl_cast_copy(mat1 &dst, const mat2 &src, const size_t numel,
                          mat3 &cast_buffer, command_queue &cq) {
  #ifdef UCL_DEBUG
  assert(dst.numel()>=numel && src.numel()>=numel);
  assert(cast_buffer.numel()>=numel);
  assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
  #endif
  if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,numel,cq);
  else
    _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,numel,
                                                      cast_buffer,cq);
}

/// Asynchronous copy of matrix/vector with cast (Device/Host transfer)
/** \param numel Number of elements (not bytes) to copy
  * \param async Perform non-blocking copy on default stream
  * \param cast_buffer Buffer on host with enough storage for casting
  * - If the data types for the two matrices are same, no cast performed
  * - Padding for 2D matrices is not considered in this routine. 
  * - Currently does not handle textures **/
template <class mat1, class mat2, class mat3>
inline void ucl_cast_copy(mat1 &dst, const mat2 &src, const size_t numel,
                          mat3 &cast_buffer, const bool async) {
  #ifdef UCL_DEBUG
  assert(dst.numel()>=numel && src.numel()>=numel);
  assert(cast_buffer.numel()>=numel);
  assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
  #endif
  if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,numel,async);
  else if (async)
    _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,numel,
                                                      cast_buffer,dst.cq());
  else
    _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,numel,
                                                      cast_buffer);
}

/// Asynchronous copy of matrix/vector (memory already allocated)
/** \param numel Number of elements (not bytes) to copy
  * - If the data types of the two matrices are not the same,
  *   casting will be performed automatically as long as the copy is
  *   not device to device. For host/device transfers, a temporary
  *   buffer is created for copy. When multiple casts occur, it is
  *   more efficient to create a permanent casting buffer that can
  *   be passed to an alternative  copy routine.
  * - Padding for 2D matrices is not considered in this routine. 
  * - Currently does not handle textures **/
template <class mat1, class mat2>
inline void ucl_copy(mat1 &dst, const mat2 &src, const size_t numel,
                     command_queue &cq) {
  #ifdef UCL_DEBUG
  assert(dst.row_size()*dst.rows()>=numel && src.row_size()*src.rows()>=numel);
  assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
  assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
  #endif
  if (mat1::MEM_TYPE==1 && mat2::MEM_TYPE==1)
    _host_host_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::hhc(dst,src,numel);
  else if ((int)mat1::DATA_TYPE!=(int)mat2::DATA_TYPE && 
      (mat1::MEM_TYPE==1 || mat2::MEM_TYPE==1)) {
    if (mat1::MEM_TYPE==1) {
      UCL_H_Vec<typename mat2::data_type> cast_buffer;
      cast_buffer.alloc(numel,dst,UCL_RW_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,numel,
                                                        cast_buffer,cq);
    } else {
      UCL_H_Vec<typename mat1::data_type> cast_buffer;
      cast_buffer.alloc(numel,dst,UCL_WRITE_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,numel,
                                                        cast_buffer,cq);
    }
  } else 
    ucl_mv_cpy(dst,src,numel*sizeof(typename mat2::data_type),cq); 
}

/// Copy matrix/vector (memory already allocated)
/** \param numel Number of elements (not bytes) to copy
  * \param async Perform non-blocking copy (ignored for host to host copy)
  * - If the data types of the two matrices are not the same,
  *   casting will be performed automatically as long as the copy is
  *   not device to device. For host/device transfers, a temporary
  *   buffer is created for copy. When multiple casts occur, it is
  *   more efficient to create a permanent casting buffer that can
  *   be passed to an alternative  copy routine.
  * - Padding for 2D matrices is not considered in this routine. 
  * - The default stream is used for asynchronous copy
  * - Currently does not handle textures **/
template <class mat1, class mat2>
inline void ucl_copy(mat1 &dst, const mat2 &src, const size_t numel,
                     const bool async) {
  #ifdef UCL_DEBUG
  assert(dst.row_size()*dst.rows()>=numel && src.row_size()*src.rows()>=numel);
  assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
  #endif
  if (mat1::MEM_TYPE==1 && mat2::MEM_TYPE==1)
    _host_host_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::hhc(dst,src,numel);
  else if (async)
    ucl_copy(dst,src,numel,dst.cq());
  else if ((int)mat1::DATA_TYPE!=(int)mat2::DATA_TYPE &&
           (mat1::MEM_TYPE==1 || mat2::MEM_TYPE==1)) {
    if (mat1::MEM_TYPE==1) {
      UCL_H_Vec<typename mat2::data_type> cast_buffer;
      cast_buffer.alloc(numel,dst,UCL_RW_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,numel,
                                                        cast_buffer);
    } else {
      UCL_H_Vec<typename mat1::data_type> cast_buffer;
      cast_buffer.alloc(numel,dst,UCL_WRITE_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,numel,
                                                        cast_buffer);
    }
  } else
    ucl_mv_cpy(dst,src,numel*sizeof(typename mat2::data_type)); 
}

// --------------------------------------------------------------------------
// - 2D COPY - SPECIFIED NUMBER OF ROWS/COLS
// --------------------------------------------------------------------------

/// Asynchronous copy subset matrix rows/cols with cast (Device/Host transfer)
/** \param async Perform non-blocking copy on default stream
  * \param cast_buffer Buffer on host with enough storage for casting
  * - If src is a vector, routine assumes row-major rows by cols copy
  * - If src is a matrix, routine will copy upper left tile of matrix 
  * - If dst is a vector, routine assumes row-major rows by cols copy
  * - If dst is a matrix, routine will copy into left tile of matrix 
  * - If the data types for the two matrices are same, no cast performed
  * - Padding for 2D matrices is not considered in this routine. 
  * - Copy from vector to matrix and vice versa allowed
  * - Currently does not handle textures **/
template <class mat1, class mat2, class mat3>
inline void ucl_cast_copy(mat1 &dst, const mat2 &src, const size_t rows,
                          const size_t cols, mat3 &cast_buffer,
                          const bool async) {
  if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,rows,cols,async);
  else if (async)
    ucl_copy(dst,src,rows,cols,dst.cq());
  else
    _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,rows,cols,
                                                      cast_buffer);
}

/// Asynchronous copy subset matrix rows,cols with cast (Device/Host transfer)
/** \param cast_buffer Buffer on host with enough storage for casting
  * - If src is a vector, routine assumes row-major rows by cols copy
  * - If src is a matrix, routine will copy upper left tile of matrix 
  * - If dst is a vector, routine assumes row-major rows by cols copy
  * - If dst is a matrix, routine will copy into upper left tile of matrix 
  * - If the data types for the two matrices are same, no cast performed
  * - Padding for 2D matrices is not considered in this routine. 
  * - Copy from vector to matrix and vice versa allowed
  * - Currently does not handle textures **/
template <class mat1, class mat2, class mat3>
inline void ucl_cast_copy(mat1 &dst, const mat2 &src, const size_t rows,
                          const size_t cols, mat3 &cast_buffer, 
                          command_queue &cq) {
  if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,rows,cols,cq);
  else 
    _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,rows,cols,
                                                      cast_buffer,cq);
}

/// Asynchronous copy of subset matrix rows,cols (memory already allocated)
/** - If src is a vector, routine assumes row-major rows by cols copy
  * - If src is a matrix, routine will copy upper left tile of matrix 
  * - If dst is a vector, routine assumes row-major rows by cols copy
  * - If dst is a matrix, routine will copy into left tile of matrix 
  * - If the data types of the two matrices are not the same,
  *   casting will be performed automatically as long as the copy is 
  *   not device to device. For host/device transfers, a temporary
  *   buffer is created for copy. When multiple casts occur, it is
  *   more efficient to create a permanent casting buffer that can
  *   be passed to an alternative copy routine.
  * - The copy should handle padding for 2D alignment correctly
  * - Copy from vector to matrix and vice versa allowed
  * - Currently does not handle textures **/
template <class mat1, class mat2>
inline void ucl_copy(mat1 &dst, const mat2 &src, const size_t rows,
                     const size_t cols, command_queue &cq) {
  if (mat1::MEM_TYPE==1 && mat2::MEM_TYPE==1)
    _host_host_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::hhc(dst,src,rows,cols);
  else if ((int)mat1::DATA_TYPE!=(int)mat2::DATA_TYPE && 
           (mat1::MEM_TYPE==1 || mat2::MEM_TYPE==1)) {
    if (mat1::MEM_TYPE==1) {
      UCL_H_Vec<typename mat2::data_type> cast_buffer;
      cast_buffer.alloc(rows*cols,dst,UCL_RW_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,rows,cols,
                                                        cast_buffer,cq);
    } else {
      UCL_H_Vec<typename mat1::data_type> cast_buffer;
      cast_buffer.alloc(rows*cols,dst,UCL_WRITE_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,rows,cols,
                                                        cast_buffer,cq);
    }
  // If we are here, at least one of the matrices must have VECTOR=0
  } else if (mat1::VECTOR) {
    #ifdef UCL_DEBUG
    assert(dst.numel()>=rows*cols && src.rows()>=rows && src.cols()>=cols);
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    #endif
    ucl_mv_cpy(dst,cols*sizeof(typename mat1::data_type),src,src.row_bytes(),
                               cols*sizeof(typename mat1::data_type),rows,
                               cq);
  } else if (mat2::VECTOR) {
    #ifdef UCL_DEBUG
    assert(src.numel()>=rows*cols && dst.rows()>=rows && dst.cols()>=cols);
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    #endif
    ucl_mv_cpy(dst,dst.row_bytes(),src,cols*sizeof(typename mat1::data_type),
               cols*sizeof(typename mat1::data_type),rows,cq);
  } else {
    #ifdef UCL_DEBUG
    assert(src.rows()>=rows && src.cols()>=cols);
    assert(dst.rows()>=rows && dst.cols()>=cols);
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    #endif
    ucl_mv_cpy(dst,dst.row_bytes(),src,src.row_bytes(),
               cols*sizeof(typename mat1::data_type),rows,cq);
  }
}

/// Copy subset of matrix rows,cols (memory already allocated)
/** \param async Perform non-blocking copy (ignored for host to host copy)
  * - If src is a vector, routine assumes row-major rows by cols copy
  * - If src is a matrix, routine will copy upper left tile of matrix 
  * - If dst is a vector, routine assumes row-major rows by cols copy
  * - If dst is a matrix, routine will copy into left tile of matrix 
  * - If the data types of the two matrices are not the same,
  *   casting will be performed automatically as long as the copy is
  *   not device to device. For host/device transfers, a temporary
  *   buffer is created for copy. When multiple casts occur, it is
  *   more efficient to create a permanent casting buffer that can
  *   be passed to an alternative  copy routine.
  * - The copy should handle padding for 2D alignment correctly
  * - Copy from vector to matrix and vice versa allowed
  * - The default stream is used for asynchronous copy
  * - Currently does not handle textures **/
template <class mat1, class mat2>
inline void ucl_copy(mat1 &dst, const mat2 &src, const size_t rows,
                     const size_t cols, const bool async) {
  if (async)
    ucl_copy(dst,src,rows,cols,dst.cq());
  else if (mat1::MEM_TYPE==1 && mat2::MEM_TYPE==1)
    _host_host_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::hhc(dst,src,rows,cols);
  else if ((int)mat1::DATA_TYPE!=(int)mat2::DATA_TYPE && 
           (mat1::MEM_TYPE==1 || mat2::MEM_TYPE==1)) {
    if (mat1::MEM_TYPE==1) {
      UCL_H_Vec<typename mat2::data_type> cast_buffer;
      cast_buffer.alloc(rows*cols,dst,UCL_RW_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,rows,cols,
                                                        cast_buffer);
    } else {
      UCL_H_Vec<typename mat1::data_type> cast_buffer;
      cast_buffer.alloc(rows*cols,dst,UCL_WRITE_OPTIMIZED);
      _ucl_cast_copy<mat1::MEM_TYPE,mat2::MEM_TYPE>::cc(dst,src,rows,cols,
                                                        cast_buffer);
    }
  // If we are here, at least one of the matrices must have VECTOR=0
  } else if (mat1::VECTOR) {
    #ifdef UCL_DEBUG
    assert(dst.numel()>=rows*cols && src.rows()>=rows && src.cols()>=cols);
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    assert(mat2::VECTOR==0);
    #endif
    ucl_mv_cpy(dst,cols*sizeof(typename mat1::data_type),src,src.row_bytes(),
                   cols*sizeof(typename mat1::data_type),rows);
  } else if (mat2::VECTOR) {
    #ifdef UCL_DEBUG
    assert(src.numel()>=rows*cols && dst.rows()>=rows && dst.cols()>=cols);
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    assert(mat1::VECTOR==0);
    #endif
    ucl_mv_cpy(dst,dst.row_bytes(),src,cols*sizeof(typename mat1::data_type),
               cols*sizeof(typename mat1::data_type),rows);
  } else {
    #ifdef UCL_DEBUG
    assert(src.rows()>=rows && src.cols()>=cols);
    assert(dst.rows()>=rows && dst.cols()>=cols);
    assert(mat1::ROW_MAJOR==1 && mat2::ROW_MAJOR==1);
    #endif
    ucl_mv_cpy(dst,dst.row_bytes(),src,src.row_bytes(),
               cols*sizeof(typename mat1::data_type),rows);
  }
}

// --------------------------------------------------------------------------
// - 1D/2D COPY
// --------------------------------------------------------------------------

/// Asynchronous copy of matrix/vector with cast (Device/Host transfer)
/** \param async Perform non-blocking copy on default stream
  * \param cast_buffer Buffer on host with enough storage for casting
  * - If the data types for the two matrices are same, no cast performed
  * - The number of bytes copied is determined by entire src data
  * - Padding for 2D matrices is not considered in this routine. 
  * - Copy from vector to matrix and vice versa allowed
  * - Currently does not handle textures **/
template <class mat1, class mat2, class mat3>
inline void ucl_cast_copy(mat1 &dst, const mat2 &src,
                          mat3 &cast_buffer, const bool async) {
  if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,async);
  else if (mat2::PADDED==1 || (mat1::PADDED==1 && mat2::VECTOR==0) )
    ucl_cast_copy(dst,src,src.rows(),src.cols(),cast_buffer,async);
  else if (mat1::PADDED==1)
    ucl_cast_copy(dst,src,dst.rows(),dst.cols(),cast_buffer,async);
  else
    ucl_cast_copy(dst,src,src.numel(),cast_buffer,async);
}

/// Asynchronous copy of matrix/vector with cast (Device/Host transfer)
/** \param cast_buffer Buffer on host with enough storage for casting
  * - If the data types for the two matrices are same, no cast performed
  * - The number of bytes copied is determined by entire src data
  * - Padding for 2D matrices is not considered in this routine. 
  * - Copy from vector to matrix and vice versa allowed
  * - Currently does not handle textures **/
template <class mat1, class mat2, class mat3>
inline void ucl_cast_copy(mat1 &dst, const mat2 &src,
                          mat3 &cast_buffer, command_queue &cq) {
  if ((int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,cq);
  else if (mat2::PADDED==1 || (mat1::PADDED==1 && mat2::VECTOR==0) )
    ucl_copy(dst,src,src.rows(),src.cols(),cast_buffer,cq);
  else if (mat1::PADDED==1)
    ucl_copy(dst,src,dst.rows(),dst.cols(),cast_buffer,cq);
  else
    ucl_copy(dst,src,src.numel(),cast_buffer,cq);
}

/// Asynchronous copy of matrix/vector (memory already allocated)
/** - The number of bytes copied is determined by entire src data
  * - If the data types of the two matrices are not the same,
  *   casting will be performed automatically as long as the copy is 
  *   not device to device. For host/device transfers, a temporary
  *   buffer is created for copy. When multiple casts occur, it is
  *   more efficient to create a permanent casting buffer that can
  *   be passed to an alternative copy routine.
  * - The copy should handle padding for 2D alignment correctly
  * - Copy from vector to matrix and vice versa allowed
  * - Currently does not handle textures **/
template <class mat1, class mat2>
inline void ucl_copy(mat1 &dst, const mat2 &src, command_queue &cq) {
  if (dst.row_bytes()==src.row_bytes() &&
      src.kind()!=UCL_VIEW && dst.kind()!=UCL_VIEW &&
      (int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,src.row_size()*src.rows(),cq);
  else if (mat2::PADDED==1 || (mat1::PADDED==1 && mat2::VECTOR==0) )
    ucl_copy(dst,src,src.rows(),src.cols(),cq);
  else if (mat1::PADDED==1)
    ucl_copy(dst,src,dst.rows(),dst.cols(),cq);
  else
    ucl_copy(dst,src,src.numel(),cq);
}

/// Copy matrix/vector (memory already allocated)
/** \param async Perform non-blocking copy (ignored for host to host copy)
  * - The number of bytes copied is determined by entire src data
  * - If the data types of the two matrices are not the same,
  *   casting will be performed automatically as long as the copy is
  *   not device to device. For host/device transfers, a temporary
  *   buffer is created for copy. When multiple casts occur, it is
  *   more efficient to create a permanent casting buffer that can
  *   be passed to an alternative  copy routine.
  * - The copy should handle padding for 2D alignment correctly
  * - Copy from vector to matrix and vice versa allowed
  * - The default stream is used for asynchronous copy
  * - Currently does not handle textures **/
template <class mat1, class mat2>
inline void ucl_copy(mat1 &dst, const mat2 &src, const bool async) {
  if (async)
    ucl_copy(dst,src,dst.cq());
  else if (dst.row_bytes()==src.row_bytes() && 
           src.kind()!=UCL_VIEW && dst.kind()!=UCL_VIEW &&
           (int)mat1::DATA_TYPE==(int)mat2::DATA_TYPE)
    ucl_copy(dst,src,src.row_size()*src.rows(),async);
  else if (mat2::PADDED==1 || (mat1::PADDED==1 && mat2::VECTOR==0) )
    ucl_copy(dst,src,src.rows(),src.cols(),async);
  else if (mat1::PADDED==1)
    ucl_copy(dst,src,dst.rows(),dst.cols(),async);
  else
    ucl_copy(dst,src,src.numel(),async);
}

#endif

