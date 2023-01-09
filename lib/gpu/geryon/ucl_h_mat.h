/***************************************************************************
                                 ucl_h_mat.h
                             -------------------
                               W. Michael Brown

  Matrix Container on Host

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

/// Matrix on Host with options for pinning (page locked)
template <class numtyp>
class UCL_H_Mat : public UCL_BaseMat {
 public:
   // Traits for copying data
   // MEM_TYPE is 0 for device, 1 for host, and 2 for image
   enum traits {
     DATA_TYPE = _UCL_DATA_ID<numtyp>::id,
     MEM_TYPE = 1,
     PADDED = 0,
     ROW_MAJOR = 1,
     VECTOR = 0
   };
   typedef numtyp data_type;

  UCL_H_Mat() : _cols(0) {
    #ifdef _OCL_MAT
    _carray=(cl_mem)(0);
    #endif
  }
  ~UCL_H_Mat() { _host_free(*this); }

  /// Construct with specied number of rows and columns
  /** \sa alloc() **/
  UCL_H_Mat(const size_t rows, const size_t cols, UCL_Device &device,
            const enum UCL_MEMOPT kind=UCL_READ_WRITE)
    { _cols=0; _kind=UCL_VIEW; alloc(rows,cols,device,kind); }

  /// Set up host matrix with specied # of rows/cols and reserve memory
  /** The kind parameter controls memory pinning as follows:
    * - UCL_READ_WRITE - Specify that you will read and write from host
    * - UCL_WRITE_ONLY - Specify that you will only write from host
    * - UCL_READ_ONLY  - Specify that you will only read from host
    * - UCL_NOT_PINNED - Memory is not pinned/page-locked on host
    * \param cq Default command queue for operations copied from another mat
    * \return UCL_SUCCESS if the memory allocation is successful **/
  template <class mat_type>
  inline int alloc(const size_t rows, const size_t cols, mat_type &cq,
                   const enum UCL_MEMOPT kind=UCL_READ_WRITE,
                   const enum UCL_MEMOPT kind2=UCL_NOT_SPECIFIED) {
    clear();

    _row_bytes=cols*sizeof(numtyp);
    int err=_host_alloc(*this,cq,_row_bytes*rows,kind,kind2);
    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes*_rows
                << " bytes on host.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _cols=cols;
    _rows=rows;
    _kind=kind;
    _end=_array+rows*cols;
    return err;
  }

  /// Set up host matrix with specied # of rows/cols and reserve memory
  /** The kind parameter controls memory pinning as follows:
    * - UCL_READ_WRITE - Specify that you will read and write from host
    * - UCL_WRITE_ONLY - Specify that you will only write from host
    * - UCL_READ_ONLY  - Specify that you will only read from host
    * - UCL_NOT_PINNED - Memory is not pinned/page-locked on host
    * \param device Used to get the default command queue for operations
    * \return UCL_SUCCESS if the memory allocation is successful **/
  inline int alloc(const size_t rows, const size_t cols, UCL_Device &device,
                   const enum UCL_MEMOPT kind=UCL_READ_WRITE,
                   const enum UCL_MEMOPT kind2=UCL_NOT_SPECIFIED) {
    clear();

    _row_bytes=cols*sizeof(numtyp);
    int err=_host_alloc(*this,device,_row_bytes*rows,kind,kind2);
    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes*_rows
                << " bytes on host.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _cols=cols;
    _rows=rows;
    _kind=kind;
    _end=_array+rows*cols;
    return err;
  }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device container on the host is not supported
    * \param stride Number of _elements_ between the start of each row **/
  template <class ucl_type>
  inline void view(ucl_type &input, const size_t rows, const size_t cols,
                   const size_t stride) {
    assert(rows==1 || stride==cols);
    clear();
    _kind=UCL_VIEW;
    _cols=cols;
    _rows=rows;
    _row_bytes=stride*sizeof(numtyp);
    this->_cq=input.cq();
    _array=input.begin();
    _end=_array+_cols;
    #ifdef _OCL_MAT
    _carray=input.cbegin();
    // When viewing outside host allocation with discrete main memory on accelerator,
    // no cl_buffer object is created to avoid unnecessary creation of device allocs
    if (_carray!=(cl_mem)(0))
      CL_SAFE_CALL(clRetainMemObject(input.cbegin()));
    CL_SAFE_CALL(clRetainCommandQueue(input.cq()));
    #endif
  }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device container on the host is not supported **/
  template <class ucl_type>
  inline void view(ucl_type &input, const size_t rows, const size_t cols)
    { view(input,rows,cols,input.row_size()); }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
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
    *   allocating container when using CUDA APIs
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view when using CUDA APIs
    * - Viewing a device container on the host is not supported **/
  template <class ucl_type>
  inline void view(ucl_type &input)
    { view(input,input.rows(),input.cols()); }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device pointer on the host is not supported
    * \param stride Number of _elements_ between the start of each row **/
  template <class ptr_type>
  inline void view(ptr_type *input, const size_t rows, const size_t cols,
                   const size_t stride, UCL_Device &dev) {
    assert(rows==1 || stride==cols);
    clear();
    _kind=UCL_VIEW;
    _cols=cols;
    _rows=rows;
    _row_bytes=stride*sizeof(numtyp);
    this->_cq=dev.cq();
    _array=input;
    _end=_array+_cols;

    #ifdef _OCL_MAT
    _host_view(*this,dev,_row_bytes*rows);
    #endif
  }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device pointer on the host is not supported **/
  template <class ptr_type>
  inline void view(ptr_type *input, const size_t rows, const size_t cols,
                   UCL_Device &dev) { view(input,rows,cols,cols,dev); }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device pointer on the host is not supported **/
  template <class ptr_type>
  inline void view(ptr_type *input, const size_t cols, UCL_Device &dev)
    { view(input,1,cols,dev); }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device container on the host is not supported
    * \param stride Number of _elements_ between the start of each row **/
  template <class ucl_type>
  inline void view_offset(const size_t offset,ucl_type &input,const size_t rows,
                          const size_t cols, const size_t stride) {
    assert(rows==1 || stride==cols);
    clear();
    _kind=UCL_VIEW;
    _cols=cols;
    _rows=rows;
    _row_bytes=stride*sizeof(numtyp);
    this->_cq=input.cq();
    _array=input.begin()+offset;
    _end=_array+_cols;
    #ifdef _OCL_MAT
    _host_view(*this,input,offset*sizeof(numtyp),_row_bytes*_rows);
    #endif
  }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device container on the host is not supported **/
  template <class ucl_type>
  inline void view_offset(const size_t offset,ucl_type &input,const size_t rows,
                          const size_t cols)
    { view_offset(offset,input,rows,cols,input.row_size()); }

  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
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
    *   allocating container when using CUDA APIs
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view
    * - Viewing a device container on the host is not supported **/
  template <class ucl_type>
  inline void view_offset(const size_t offset, ucl_type &input) {
    if (input.rows()==1)
      view_offset(offset,input,1,input.cols()-offset);
    else
      view_offset(offset,input,input.rows()-offset/input.row_size(),
                  input.cols());
  }

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
    *   allocating container when using CUDA APIs
    * - Viewing a device pointer on the host is not supported
    * \param stride Number of _elements_ between the start of each row **/
  template <class ptr_type>
  inline void view_offset(const size_t offset,ptr_type *input,const size_t rows,
                          const size_t cols,const size_t stride,UCL_Device &dev)
    { view(input+offset,rows,cols,stride,dev); }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container when using CUDA APIs
    * - Viewing a device pointer on the host is not supported **/
  template <class ptr_type>
  inline void view_offset(const size_t offset, ptr_type *input,
                          const size_t cols, UCL_Device &dev)
    { view(input+offset,1,cols,dev); }

  /// Free memory and set size to 0
  inline void clear()
    { _host_free(*this); _cols=0; _kind=UCL_VIEW; }

  /// Resize the allocation to rows x cols elements
  /** \note Cannot be used on views **/
  inline int resize(const int rows, const int cols) {
    assert(_kind!=UCL_VIEW);

    _row_bytes=cols*sizeof(numtyp);
    int err=_host_resize(*this,_row_bytes*rows);
    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes*_rows
                << " bytes on host.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _cols=cols;
    _rows=rows;
    _end=_array+rows*cols;
    return err;
  }

  /// Resize (only if bigger) the allocation to contain rows x cols elements
  /** \note Cannot be used on views **/
  inline int resize_ib(const int rows, const int cols)
    { if (cols>_cols || rows>_rows) return resize(rows,cols);
      else return UCL_SUCCESS; }

  /// Set each element to zero
  inline void zero() { _host_zero(_array,_rows*row_bytes()); }
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
  inline size_t numel() const { return _rows*_cols; }
  /// Get the number of rows
  inline size_t rows() const { return _rows; }
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
    { return _array[row*_cols+col]; }
  /// 2D access (row should always be 0)
  inline const numtyp & operator()(const int row, const int col) const
    { return _array[row*_cols+col]; }

  /// Returns pointer to memory pointer for allocation on host
  inline numtyp ** host_ptr() { return &_array; }

  /// Return the offset (in elements) from begin() pointer where data starts
  /** \note Always 0 for host matrices and CUDA APIs **/
  inline size_t offset() const { return 0; }
  /// Return the offset (in bytes) from begin() pointer where data starts
  /** \note Always 0 for host matrices and CUDA APIs **/
  inline size_t byteoff() const { return 0; }

  #ifdef _OCL_MAT
  /// Returns an API specific device pointer (cl_mem& for OpenCL, void ** for CUDA)
  inline device_ptr & cbegin() { return _carray; }
  /// Returns an API specific device pointer (cl_mem& for OpenCL, void ** for CUDA)
  inline const device_ptr & cbegin() const { return _carray; }
  #else
  /// Returns an API specific device pointer (cl_mem& for OpenCL, void ** for CUDA)
  inline void ** cbegin() { return (void **)&_array; }
  /// Returns an API specific device pointer (cl_mem& for OpenCL, void ** for CUDA)
  inline const void ** cbegin() const { return (const void **)&_array; }
  #endif

 private:
  numtyp *_array, *_end;
  size_t _row_bytes, _rows, _cols;

  #ifdef _OCL_MAT
  device_ptr _carray;
  #endif
};

#endif

