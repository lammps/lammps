/***************************************************************************
                                 ucl_h_vec.h
                             -------------------
                               W. Michael Brown

  Vector Container on Host

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Thu Jun 25 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

// Only allow this file to be included by CUDA and OpenCL specific headers
#ifdef _UCL_MAT_ALLOW

/// Row Vector on Host with options for pinning (page locked)
template <class numtyp>
class UCL_H_Vec : public UCL_BaseMat {
 public:
   // Traits for copying data
   // MEM_TYPE is 0 for device, 1 for host, and 2 for image
   enum traits {
     DATA_TYPE = _UCL_DATA_ID<numtyp>::id,
     MEM_TYPE = 1,
     PADDED = 0,
     ROW_MAJOR = 1,
     VECTOR = 1
   };
   typedef numtyp data_type; 
   
  UCL_H_Vec() : _kind(UCL_VIEW), _cols(0) {
    #ifdef _OCL_MAT
    _carray=(cl_mem)(0);
    #endif
  }
  ~UCL_H_Vec() { if (_kind!=UCL_VIEW) _host_free(*this,_kind); }
  
  /// Construct with n columns
  /** \sa alloc() **/
  UCL_H_Vec(const size_t n, UCL_Device &device, 
            const enum UCL_MEMOPT kind=UCL_RW_OPTIMIZED) 
    { _cols=0; _kind=UCL_VIEW; alloc(n,device,kind); }
  
  /// Set up host vector with 'cols' columns and reserve memory
  /** The kind parameter controls memory pinning as follows:
    * - UCL_NOT_PINNED      - Memory is not pinned
    * - UCL_WRITE_OPTIMIZED - Memory can be pinned (write-combined)
    * - UCL_RW_OPTIMIZED    - Memory can be pinned 
    * \param cq Default command queue for operations copied from another mat
    * \return UCL_SUCCESS if the memory allocation is successful **/
  template <class mat_type>
  inline int alloc(const size_t cols, mat_type &cq,
                   const enum UCL_MEMOPT kind=UCL_RW_OPTIMIZED) {
    clear();

    _row_bytes=cols*sizeof(numtyp);
    int err=_host_alloc(*this,cq,_row_bytes,kind);

    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes
                << " bytes on host.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _cols=cols;
    _kind=kind;
    _end=_array+cols;
    return err;
  }    

  /// Set up host vector with 'cols' columns and reserve memory
  /** The kind parameter controls memory pinning as follows:
    * - UCL_NOT_PINNED      - Memory is not pinned
    * - UCL_WRITE_OPTIMIZED - Memory can be pinned (write-combined)
    * - UCL_RW_OPTIMIZED    - Memory can be pinned 
    * \param device Used to get the default command queue for operations
    * \return UCL_SUCCESS if the memory allocation is successful **/
  inline int alloc(const size_t cols, UCL_Device &device,
                   const enum UCL_MEMOPT kind=UCL_RW_OPTIMIZED) {
    clear();

    _row_bytes=cols*sizeof(numtyp);
    int err=_host_alloc(*this,device,_row_bytes,kind);

    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes
                << " bytes on host.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif 
      _row_bytes=0;
      return err;
    }

    _cols=cols;
    _kind=kind;
    _end=_array+cols;
    return err;
  }
  
  /// Return the type of memory allocation
  /** Returns UCL_READ_WRITE, UCL_WRITE_ONLY, UCL_READ_ONLY, or UCL_VIEW **/ 
  inline enum UCL_MEMOPT kind() const { return _kind; }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device container on the host is not supported **/ 
  template <class ucl_type>
  inline void view(ucl_type &input, const size_t rows, const size_t cols) {
    #ifdef UCL_DEBUG
    assert(rows==1);
    #endif
    clear();
    _kind=UCL_VIEW;
    _cols=cols;
    _row_bytes=_cols*sizeof(numtyp);
    this->_cq=input.cq();
    _array=input.begin();
    _end=_array+_cols;
    #ifdef _OCL_MAT
    _carray=input.cbegin();
    #endif
  }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device container on the host is not supported 
    * \param stride Number of _elements_ between the start of each row **/ 
  template <class ucl_type>
  inline void view(ucl_type &input, const size_t rows, const size_t cols,
                   const size_t stride) { view(input,rows,cols); }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view 
    * - Viewing a device container on the host is not supported **/ 
  template <class ucl_type>
  inline void view(ucl_type &input, const size_t cols)
    { view(input,1,cols); }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view 
    * - Viewing a device container on the host is not supported **/ 
  template <class ucl_type>
  inline void view(ucl_type &input) 
    { view(input,input.rows()*input.row_size()); }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device pointer on the host is not supported **/ 
  template <class ptr_type>
  inline void view(ptr_type *input, const size_t rows, const size_t cols,
                   UCL_Device &dev) {
    #ifdef UCL_DEBUG
    assert(rows==1);
    #endif
    clear();
    _kind=UCL_VIEW;
    _cols=cols;
    _row_bytes=_cols*sizeof(numtyp);
    this->_cq=dev.cq();
    _array=input;
    _end=_array+_cols;
    
    #ifdef _OCL_MAT
    _host_alloc(*this,dev,_row_bytes,UCL_VIEW);
    #endif 
  }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device pointer on the host is not supported 
    * \param stride Number of _elements_ between the start of each row **/ 
  template <class ptr_type>
  inline void view(ptr_type *input, const size_t rows, const size_t cols,
                   const size_t stride, UCL_Device &dev) 
    { view(input,rows,cols,stride); }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device pointer on the host is not supported **/ 
  template <class ptr_type>
  inline void view(ptr_type *input, const size_t cols, UCL_Device &dev)
    { view(input,1,cols,dev); }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device container on the host is not supported **/ 
  template <class ucl_type>
  inline void view_offset(const size_t offset,ucl_type &input,const size_t rows,
                          const size_t cols) {
    #ifdef UCL_DEBUG
    assert(rows==1);
    #endif
    clear();
    _kind=UCL_VIEW;
    _cols=cols;
    _row_bytes=_cols*sizeof(numtyp);
    this->_cq=input.cq();
    _array=input.begin()+offset;
    _end=_array+_cols;
    #ifdef _OCL_MAT
    _host_alloc(*this,input,_row_bytes,UCL_VIEW);
    #endif
  }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device container on the host is not supported 
    * \param stride Number of _elements_ between the start of each row **/ 
  template <class ucl_type>
  inline void view_offset(const size_t offset,ucl_type &input,const size_t rows,
                          const size_t cols, const size_t stride) 
    { view_offset(offset,input,rows,cols); }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view 
    * - Viewing a device container on the host is not supported **/ 
  template <class ucl_type>
  inline void view_offset(const size_t offset,ucl_type &input,const size_t cols)
    { view_offset(offset,input,1,cols); }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view 
    * - Viewing a device container on the host is not supported **/ 
  template <class ucl_type>
  inline void view_offset(const size_t offset, ucl_type &input) 
    { view_offset(offset,input,input.rows()*input.row_size()-offset); }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device pointer on the host is not supported **/ 
  template <class ptr_type>
  inline void view_offset(const size_t offset,ptr_type *input,const size_t rows,
                          const size_t cols, UCL_Device &dev)
    { view(input+offset,rows,cols,dev); }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device pointer on the host is not supported 
    * \param stride Number of _elements_ between the start of each row **/ 
  template <class ptr_type>
  inline void view_offset(const size_t offset,ptr_type *input,const size_t rows,
                          const size_t cols,const size_t stride,UCL_Device &dev) 
    { view(input+offset,rows,cols,stride,dev); }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - Viewing a device pointer on the host is not supported **/ 
  template <class ptr_type>
  inline void view_offset(const size_t offset, ptr_type *input, 
                          const size_t cols, UCL_Device &dev)
    { view(input+offset,1,cols,dev); }
  
  /// Free memory and set size to 0
  inline void clear() 
    { if (_kind!=UCL_VIEW) {_kind=UCL_VIEW; _cols=0; _host_free(*this,_kind);}}

  /// Resize the allocation to contain cols elements
  /** \note Cannot be used on views **/
  inline int resize(const int cols) {
    assert(_kind!=UCL_VIEW);
    _row_bytes=cols*sizeof(numtyp);
    int err=_host_resize(*this,_row_bytes);
    
    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes
                << " bytes on host.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _cols=cols;
    _end=_array+cols;
    return err;
  }
    
  /// Resize (only if bigger) the allocation to contain cols elements
  /** \note Cannot be used on views **/
  inline int resize_ib(const int cols)
    { if (cols>_cols) return resize(cols); else return UCL_SUCCESS; }

  /// Set each element to zero
  inline void zero() { _host_zero(_array,row_bytes()); }
  
  /// Set first n elements to zero
  inline void zero(const int n) { _host_zero(_array,n*sizeof(numtyp)); }

  /// Get host pointer to first element
  inline numtyp * begin() { return _array; }
  /// Get host pointer to first element
  inline const numtyp * begin() const { return _array; }
  /// Get host pointer to one past last element
  inline numtyp * end() { return _end; }
  /// Get host pointer to one past last element
  inline const numtyp * end() const { return _end; }

  /// Get the number of elements
  inline size_t numel() const { return _cols; }
  /// Get the number of rows
  inline size_t rows() const { return 1; }
  /// Get the number of columns
  inline size_t cols() const { return _cols; }
  ///Get the size of a row (including any padding) in elements
  inline size_t row_size() const { return _cols; }
  /// Get the size of a row (including any padding) in bytes
  inline size_t row_bytes() const { return _row_bytes; }
  /// Get the size in bytes of 1 element
  inline int element_size() const { return sizeof(numtyp); }
    
  /// Get element at index i
  inline numtyp & operator[](const int i) { return _array[i]; }
  /// Get element at index i
  inline const numtyp & operator[](const int i) const { return _array[i]; }
  /// 2D access (row should always be 0) 
  inline numtyp & operator()(const int row, const int col) 
    { return _array[col]; }
  /// 2D access (row should always be 0) 
  inline const numtyp & operator()(const int row, const int col) const
    { return _array[col]; }
  
  /// Returns pointer to memory pointer for allocation on host
  inline numtyp ** host_ptr() { return &_array; }
  
  /// Return the offset (in elements) from begin() pointer where data starts
  /** \note Always 0 for host matrices and CUDA APIs **/
  inline size_t offset() const { return 0; }
  /// Return the offset (in bytes) from begin() pointer where data starts
  /** \note Always 0 for host matrices and CUDA APIs **/
  inline size_t byteoff() const { return 0; }
  
  #ifdef _OCL_MAT
  /// For OpenCL, returns a reference to the cl_mem object
  inline device_ptr & cbegin() { return _carray; }
  /// For OpenCL, returns a reference to the cl_mem object
  inline const device_ptr & cbegin() const { return _carray; }
  #endif
  
 private:
  enum UCL_MEMOPT _kind;
  numtyp *_array, *_end;
  size_t _row_bytes, _cols;

  #ifdef _OCL_MAT
  device_ptr _carray;
  #endif
};

#endif

