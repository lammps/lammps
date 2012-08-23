/***************************************************************************
                                 ucl_d_vec.h
                             -------------------
                               W. Michael Brown

  Vector Container on Device

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

/// Row vector on device 
template <class numtyp>
class UCL_D_Vec : public UCL_BaseMat {
 public:
  // Traits for copying data
  // MEM_TYPE is 0 for device, 1 for host, and 2 for image
  enum traits {
    DATA_TYPE = _UCL_DATA_ID<numtyp>::id,
    MEM_TYPE = 0,
    PADDED = 0,
    ROW_MAJOR = 1,
    VECTOR = 1
  };
  typedef numtyp data_type; 

  UCL_D_Vec() : _cols(0), _kind(UCL_VIEW) {}
  ~UCL_D_Vec() { if (_kind!=UCL_VIEW) _device_free(*this); }

  /// Construct with n columns
  /** \sa alloc() **/
  UCL_D_Vec(const size_t n, UCL_Device &device,
            const enum UCL_MEMOPT kind=UCL_READ_WRITE) : 
    _cols(0), _kind(UCL_VIEW) { alloc(n,device,kind); }

  /// Set up host vector with 'cols' columns and reserve memory
  /** The kind parameter controls memory optimizations as follows:
    * - UCL_READ_WRITE - Specify that you will read and write in kernels
    * - UCL_WRITE_ONLY - Specify that you will only write in kernels
    * - UCL_READ_ONLY  - Specify that you will only read in kernels
    * \param cq Default command queue for operations copied from another mat
    * \return UCL_SUCCESS if the memory allocation is successful **/
  template <class mat_type>
  inline int alloc(const size_t cols, mat_type &cq,
                   const enum UCL_MEMOPT kind=UCL_READ_WRITE) {
                        
    clear();

    _row_bytes=cols*sizeof(numtyp);
    int err=_device_alloc(*this,cq,_row_bytes,kind);
    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes
                << " bytes on device.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _kind=kind;
    _cols=cols;
    #ifndef _UCL_DEVICE_PTR_MAT
    _end=_array+cols;
    #endif
    #ifdef _OCL_MAT
    _offset=0;
    #endif
    return err; 
  }    

  /// Set up host vector with 'cols' columns and reserve memory
  /** The kind parameter controls memory optimizations as follows:
    * - UCL_READ_WRITE - Specify that you will read and write in kernels
    * - UCL_WRITE_ONLY - Specify that you will only write in kernels
    * - UCL_READ_ONLY  - Specify that you will only read in kernels
    * \param device Used to get the default command queue for operations
    * \return UCL_SUCCESS if the memory allocation is successful **/
  inline int alloc(const size_t cols, UCL_Device &device,
                   const enum UCL_MEMOPT kind=UCL_READ_WRITE) {
    clear();
    _row_bytes=cols*sizeof(numtyp);
    int err=_device_alloc(*this,device,_row_bytes,kind);
    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes
                << " bytes on device.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _kind=kind;
    _cols=cols;
    #ifndef _UCL_DEVICE_PTR_MAT
    _end=_array+cols;
    #endif
    #ifdef _OCL_MAT
    _offset=0;
    #endif
    return err; 
  }

  /// Return the type of memory allocation
  /** Returns UCL_READ_WRITE, UCL_WRITE_ONLY, UCL_READ_ONLY, or UCL_VIEW **/ 
  inline enum UCL_MEMOPT kind() const { return _kind; }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container **/
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
    #ifdef _OCL_MAT
    _offset=0;
    _array=input.cbegin();
    #else
    _device_view(&_array,input.begin());
    #endif
    
    #ifndef _UCL_DEVICE_PTR_MAT
    _end=_array+_cols;
    #endif
  }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
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
    *   will be used for view **/
  template <class ucl_type>
  inline void view(ucl_type &input, const size_t cols)
    { view(input,1,cols); }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view **/
  template <class ucl_type>
  inline void view(ucl_type &input) 
    { view(input,input.rows()*input.row_size()); }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container **/
  template <class ptr_type>
  inline void view(ptr_type input, const size_t rows, const size_t cols,
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
    #ifndef _UCL_DEVICE_PTR_MAT
    _end=_array+_cols;
    #endif
    #ifdef _OCL_MAT
    _offset=0;
    #endif
  }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * \param stride Number of _elements_ between the start of each row **/ 
  template <class ptr_type>
  inline void view(ptr_type input, const size_t rows, const size_t cols,
                   const size_t stride, UCL_Device &dev) 
    { view(input,rows,cols,stride); }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container **/
  template <class ptr_type>
  inline void view(ptr_type input, const size_t cols, UCL_Device &dev)
    { view(input,1,cols,dev); }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container **/
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
    #ifdef _OCL_MAT
    _array=input.begin();
    _offset=offset;
    #else
    _device_view(&_array,input.begin(),offset,sizeof(numtyp));
    #endif
    
    #ifndef _UCL_DEVICE_PTR_MAT
    _end=_array+_cols;
    #endif
  }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
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
    *   will be used for view **/
  template <class ucl_type>
  inline void view_offset(const size_t offset,ucl_type &input,const size_t cols)
    { view_offset(offset,input,1,cols); }
  
  /// Do not allocate memory, instead use an existing allocation from Geryon
  /** This function must be passed a Geryon vector or matrix container.
    * No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * - If a matrix is used a input, all elements (including padding)
    *   will be used for view **/
  template <class ucl_type>
  inline void view_offset(const size_t offset, ucl_type &input) 
    { view_offset(offset,input,input.rows()*input.row_size()-offset); }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container **/
  template <class ptr_type>
  inline void view_offset(const size_t offset,ptr_type input,const size_t rows,
                          const size_t cols, UCL_Device &dev) {
    #ifdef UCL_DEBUG
    assert(rows==1);
    #endif
    clear();
    _kind=UCL_VIEW;
    _cols=cols;
    _row_bytes=_cols*sizeof(numtyp);
    this->_cq=dev.cq();
    
    #ifdef _OCL_MAT
    _array=input;
    _offset=offset;
    #else
    #ifdef _UCL_DEVICE_PTR_MAT
    _array=input+offset*sizeof(numtyp);
    #else
    _array=input+offset;
    #endif
    #endif
    
    #ifndef _UCL_DEVICE_PTR_MAT
    _end=_array+_cols;
    #endif
  }
  
  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container 
    * \param stride Number of _elements_ between the start of each row **/ 
  template <class ptr_type>
  inline void view_offset(const size_t offset,ptr_type input,const size_t rows,
                          const size_t cols,const size_t stride,UCL_Device &dev) 
    { view_offset(offset,input,rows,cols,stride); }

  /// Do not allocate memory, instead use an existing allocation
  /** - No memory is freed when the object is destructed.
    * - The view does not prevent the memory from being freed by the
    *   allocating container **/
  template <class ptr_type>
  inline void view_offset(const size_t offset, ptr_type input, 
                          const size_t cols, UCL_Device &dev)
    { view_offset(offset,input,1,cols,dev); }
  
  /// Free memory and set size to 0
  inline void clear() 
    { if (_kind!=UCL_VIEW) { _cols=0; _kind=UCL_VIEW; _device_free(*this); } }

  /// Resize the allocation to contain cols elements
  /** \note Cannot be used on views **/
  inline int resize(const int cols) {
    assert(_kind!=UCL_VIEW);

    _row_bytes=cols*sizeof(numtyp);
    int err=_device_resize(*this,_row_bytes);
    if (err!=UCL_SUCCESS) {
      #ifndef UCL_NO_EXIT
      std::cerr << "UCL Error: Could not allocate " << _row_bytes
                << " bytes on device.\n";
      _row_bytes=0;
      UCL_GERYON_EXIT;
      #endif
      _row_bytes=0;
      return err;
    }

    _cols=cols;
    #ifndef _UCL_DEVICE_PTR_MAT
    _end=_array+cols;
    #endif
    #ifdef _OCL_MAT
    _offset=0;
    #endif
    return err; 
  }
    
  /// Resize (only if bigger) the allocation to contain cols elements
  /** \note Cannot be used on views **/
  inline int resize_ib(const int cols)
    { if (cols>_cols) return resize(cols); else return UCL_SUCCESS; }

  /// Set each element to zero
  inline void zero() { _device_zero(*this,row_bytes()); }

  /// Set first n elements to zero
  inline void zero(const int n) { _device_zero(*this,n*sizeof(numtyp)); }

  #ifdef _UCL_DEVICE_PTR_MAT
  /// For OpenCL, returns a (void *) device pointer to memory allocation
  inline device_ptr & begin() { return _array; }
  /// For OpenCL, returns a (void *) device pointer to memory allocation
  inline const device_ptr & begin() const { return _array; }
  #else
  /// For CUDA-RT, get device pointer to first element
  inline numtyp * & begin() { return _array; }
  /// For CUDA-RT, get device pointer to first element
  inline numtyp * const & begin() const { return _array; }
  /// For CUDA-RT, get device pointer to one past last element
  inline numtyp * end() { return _end; }
  /// For CUDA-RT, get device pointer to one past last element
  inline numtyp * end() const { return _end; }
  #endif
  
  #ifdef _UCL_DEVICE_PTR_MAT
  /// Returns an API specific device pointer
  /** - For OpenCL, returns a &cl_mem object
    * - For CUDA Driver, returns a &CUdeviceptr
    * - For CUDA-RT, returns void** **/
  inline device_ptr & cbegin() { return _array; }
  /// Returns an API specific device pointer
  /** - For OpenCL, returns a &cl_mem object
    * - For CUDA Driver, returns a &CUdeviceptr
    * - For CUDA-RT, returns void** **/
  inline const device_ptr & cbegin() const { return _array; }
  #else
  /// Returns an API specific device pointer
  /** - For OpenCL, returns a &cl_mem object
    * - For CUDA Driver, returns a &CUdeviceptr
    * - For CUDA-RT, returns numtyp** **/
  inline numtyp ** cbegin() { return &_array; }
  /// Returns an API specific device pointer
  /** - For OpenCL, returns a &cl_mem object
    * - For CUDA Driver, returns a &CUdeviceptr
    * - For CUDA-RT, returns numtyp** **/
  inline const numtyp ** cbegin() const { return &_array; }
  /// For CUDA-RT, allocate row vector and bind texture
  inline void safe_alloc(const size_t cols, UCL_Device &dev,
                         textureReference *t) 
    { alloc(cols,dev); assign_texture(t); bind(); }
  /// For CUDA-RT, assign a texture to matrix
  inline void assign_texture(textureReference *t) { _tex_ptr=t; }  
  /// For CUDA-RT, bind to texture
  inline void bind() {
    cuda_gb_get_channel<numtyp>(_channel);
    (*_tex_ptr).addressMode[0] = cudaAddressModeClamp;
    (*_tex_ptr).addressMode[1] = cudaAddressModeClamp;
    (*_tex_ptr).filterMode = cudaFilterModePoint;
    (*_tex_ptr).normalized = false;
    CUDA_SAFE_CALL(cudaBindTexture(NULL,_tex_ptr,_array,&_channel));
  }
  /// For CUDA-RT, unbind texture
  inline void unbind() { CUDA_SAFE_CALL(cudaUnbindTexture(_tex_ptr)); }
  #endif

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
  
  #ifdef _OCL_MAT
  /// Return the offset (in elements) from begin() pointer where data starts
  /** \note Always 0 for host matrices and CUDA APIs **/
  inline size_t offset() const { return _offset; }
  #else
  /// Return the offset (in elements) from begin() pointer where data starts
  /** \note Always 0 for host matrices and CUDA APIs **/
  inline size_t offset() const { return 0; }
  #endif

  /// Return the offset (in bytes) from begin() pointer where data starts
  /** \note Always 0 for host matrices and CUDA APIs **/
  inline size_t byteoff() const { return offset()*sizeof(numtyp); }

 private:
  size_t _row_bytes, _row_size, _rows, _cols;
  enum UCL_MEMOPT _kind;
  
  #ifdef _UCL_DEVICE_PTR_MAT
  device_ptr _array;
  #else
  numtyp *_array,*_end;
  cudaChannelFormatDesc _channel;
  textureReference *_tex_ptr;
  #endif

  #ifdef _OCL_MAT
  size_t _offset;
  #endif
};

#endif

