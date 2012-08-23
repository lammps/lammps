/***************************************************************************
                                 ucl_matrix.h
                             -------------------
                               W. Michael Brown

  Matrix Container on Host

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Thu May 10 2012
    copyright            : (C) 2012 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   This software is distributed under the Simplified BSD License.
   ----------------------------------------------------------------------- */

// Only allow this file to be included by CUDA and OpenCL specific headers
#ifdef _UCL_MAT_ALLOW

/// Matrix S-Object
template <class hosttype, class devtype>
class UCL_Matrix {
 public:
  // Traits for copying data
  // MEM_TYPE is 0 for device, 1 for host, and 2 for image
  enum traits {
    DATA_TYPE = _UCL_DATA_ID<hosttype>::id,
    MEM_TYPE = 1,
    PADDED = 0,
    ROW_MAJOR = 1,
    VECTOR = 0
  };
  typedef hosttype data_type; 

  /// Host Allocation
  UCL_H_Mat<hosttype> host;
  
  /// Device Allocation
  UCL_D_Mat<devtype> device;

  UCL_Matrix() { }
  ~UCL_Matrix() { }
  
  /// Construct with specied number of rows and columns
  /** \sa alloc() **/
  UCL_Matrix(const size_t rows, const size_t cols, UCL_Device &acc, 
             const enum UCL_MEMOPT kind1=UCL_RW_OPTIMIZED,
             const enum UCL_MEMOPT kind2=UCL_READ_WRITE)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        alloc(host,device,_buffer,rows,cols,acc,kind1,kind2); }
  
  /// Set up host matrix with specied # of rows/cols and reserve memory
  /** The kind1 parameter controls memory pinning as follows:
    * - UCL_NOT_PINNED      - Memory is not pinned
    * - UCL_WRITE_OPTIMIZED - Memory can be pinned (write-combined)
    * - UCL_RW_OPTIMIZED    - Memory can be pinned 
    * The kind2 parameter controls memory optimizations as follows:
    * - UCL_READ_WRITE - Specify that you will read and write in kernels
    * - UCL_WRITE_ONLY - Specify that you will only write in kernels
    * - UCL_READ_ONLY  - Specify that you will only read in kernels
    * \note When passing a command queue instead of a device, the device
    *       allocation is always performed. Even if the device shares memory
    *       with the host.
    * \param cq Default command queue for operations copied from another mat
    * \return UCL_SUCCESS if the memory allocation is successful **/
  template <class mat_type>
  inline int alloc(const size_t rows, const size_t cols, mat_type &cq,
                   const enum UCL_MEMOPT kind1=UCL_RW_OPTIMIZED,
                   const enum UCL_MEMOPT kind2=UCL_READ_WRITE)
    { return _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        alloc(host,device,_buffer,rows,cols,cq,kind1,kind2); }
  
  /// Set up host matrix with specied # of rows/cols and reserve memory
  /** The kind1 parameter controls memory pinning as follows:
    * - UCL_NOT_PINNED      - Memory is not pinned
    * - UCL_WRITE_OPTIMIZED - Memory can be pinned (write-combined)
    * - UCL_RW_OPTIMIZED    - Memory can be pinned 
    * The kind2 parameter controls memory optimizations as follows:
    * - UCL_READ_WRITE - Specify that you will read and write in kernels
    * - UCL_WRITE_ONLY - Specify that you will only write in kernels
    * - UCL_READ_ONLY  - Specify that you will only read in kernels
    * \param device Used to get the default command queue for operations
    * \return UCL_SUCCESS if the memory allocation is successful **/
  inline int alloc(const size_t rows, const size_t cols, UCL_Device &acc,
                   const enum UCL_MEMOPT kind1=UCL_RW_OPTIMIZED,
                   const enum UCL_MEMOPT kind2=UCL_READ_WRITE)
    { return _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        alloc(host,device,_buffer,rows,cols,acc,kind1,kind2); }
  
  /// Free memory and set size to 0
  inline void clear() 
    { host.clear(); device.clear(); }

  /// Resize the allocation to contain cols elements
  inline int resize(const int rows, const int cols) {
    assert(host.kind()!=UCL_VIEW);
    int err=host.resize(rows,cols);
    if (err!=UCL_SUCCESS)
      return err;
    return _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
      dev_resize(device,host,_buffer,rows,cols);
  }
    
  /// Resize (only if bigger) the allocation to contain cols elements
  inline int resize_ib(const int new_rows, const int new_cols)
    { if (new_rows>rows() || new_cols>cols()) return resize(new_rows,new_cols); 
      else return UCL_SUCCESS; }

  /// Set each element to zero
  inline void zero() { host.zero(); device.zero(); }
  
  /// Set first n elements to zero
  inline void zero(const int n) { host.zero(n); device.zero(n); }

  /// Get the number of elements
  inline size_t numel() const { return host.numel(); }
  /// Get the number of rows
  inline size_t rows() const { return host.rows(); }
  /// Get the number of columns
  inline size_t cols() const { return host.cols(); }
  /// Get the memory usage (bytes) of the s-object (including any buffers)
  inline size_t host_mem_usage() 
    { return host.row_bytes()*host.rows()+_buffer.row_bytes()*_buffer.rows(); }
  /// Get the memory usage (bytes) of the s-object (including any buffers)
  inline size_t device_mem_usage() 
    { return device.row_bytes()*device.rows(); }
    
  /// Get element at index i
  inline hosttype & operator[](const int i) { return host[i]; }
  /// Get element at index i
  inline const hosttype & operator[](const int i) const { return host[i]; }
  /// 2D access (row should always be 0) 
  inline hosttype & operator()(const int row, const int col) 
    { return host(row,col); }
  /// 2D access (row should always be 0) 
  inline const hosttype & operator()(const int row, const int col) const
    { return host(row,col); }
  
  /// Returns pointer to memory pointer for allocation on host
  inline hosttype ** host_ptr() { return host.host_ptr(); }
  
  /// Return the default command queue/stream associated with this data
  inline command_queue & cq() { return host.cq(); }
  /// Block until command_queue associated with matrix is complete
  inline void sync() { host.sync(); }

  ///Get the size of a row on the host (including any padding) in elements
  inline size_t row_size() const { return host.row_size(); }
  /// Get the size of a row on the host(including any padding) in bytes
  inline size_t row_bytes() const { return host.row_bytes(); }
  /// Get the size on the host in bytes of 1 element
  inline int element_size() const { return sizeof(hosttype); }


  /// Update the allocation on the host asynchronously
  inline void update_host() 
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(host,device,_buffer,true); }
  /// Update the allocation on the host (true for asynchronous copy)
  inline void update_host(const bool async)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(host,device,_buffer,async); }
  /// Update the allocation on the host (using command queue)
  inline void update_host(command_queue &cq)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(host,device,_buffer,cq); }
  /// Update the first n elements on the host (true for asynchronous copy)
  inline void update_host(const int n, const bool async)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(host,device,n,_buffer,async); }
  /// Update the first n elements on the host (using command queue)
  inline void update_host(const int n, command_queue &cq)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(host,device,n,_buffer,cq); }
  /// Update slice on the host (true for asynchronous copy)
  inline void update_host(const int rows, const int cols, const bool async)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(host,device,rows,cols,_buffer,async); }
  /// Update slice on the host (using command queue)
  inline void update_host(const int rows, const int cols, command_queue &cq)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(host,device,rows,cols,_buffer,cq); }


  /// Update the allocation on the device asynchronously
  inline void update_device() 
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(device,host,_buffer,true); }
  /// Update the allocation on the device (true for asynchronous copy)
  inline void update_device(const bool async)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(device,host,_buffer,async); }
  /// Update the allocation on the device (using command queue)
  inline void update_device(command_queue &cq)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(device,host,_buffer,cq); }
  /// Update the first n elements on the device (true for asynchronous copy)
  inline void update_device(const int n, const bool async)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(device,host,n,_buffer,async); }
  /// Update the first n elements on the device (using command queue)
  inline void update_device(const int n, command_queue &cq)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(device,host,n,_buffer,cq); }
  /// Update slice on the device (true for asynchronous copy)
  inline void update_device(const int rows, const int cols, const bool async)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(device,host,rows,cols,_buffer,async); }
  /// Update slice on the device (using command queue)
  inline void update_device(const int rows, const int cols, command_queue &cq)
    { _ucl_s_obj_help< ucl_same_type<hosttype,devtype>::ans >::
        copy(device,host,rows,cols,_buffer,cq); }


 private:
  UCL_H_Mat<devtype> _buffer;
};

#endif

