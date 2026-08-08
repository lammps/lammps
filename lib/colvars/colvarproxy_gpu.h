#ifndef COLVARPROXY_GPU_H
#define COLVARPROXY_GPU_H

#include "colvar_gpu_support.h"
#include "colvarmodule.h"

/**
 * @file colvarproxy_gpu.h
 * @brief Declaration of the class for GPU memory management
 */

/**
 * @brief Class for managing GPU memory allocation and data transfer
 */
class colvarproxy_gpu {
public:
  /// \brief Constructor
  colvarproxy_gpu(): support_gpu(false) {}
  /// \brief Whether the proxy supports GPU
  bool has_gpu_support() const {
    return support_gpu;
  }
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP) || defined (COLVARS_SYCL)
  /// \brief Get the default CUDA stream from the proxy
  virtual cudaStream_t get_default_stream() {return (cudaStream_t)0;}
  /**
   * @brief Template function to allocate host-pinned memory
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated host-pinned memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int allocate_host(T **pp, const size_t len) {
    return allocate_host_T((void **)pp, len, sizeof(T));
  }
  /**
   * @brief Template function to deallocate host-pinned memory
   *
   * @tparam T The type of elements to deallocate
   * @param[in,out] pp Pointer to the pointer that holds the allocated host-pinned memory
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int deallocate_host(T **pp) {
    return deallocate_host_T((void **)pp);
  }
  /**
   * @brief Template function to allocate device memory
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated device memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int allocate_device(T **pp, const size_t len) {
    return allocate_device_T((void **)pp, len, sizeof(T));
  }
  /**
   * @brief Template function to reallocate device memory
   *
   * This function first deallocates any existing memory pointed to by `*pp`,
   * then allocates new device memory for `len` elements of type `T`.
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated device memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int reallocate_device(T **pp, const size_t len) {
    int error_code = COLVARS_OK;
    error_code |= deallocate_device(pp);
    error_code |= allocate_device_T((void **)pp, len, sizeof(T));
    return error_code;
  }
  /**
   * @brief Template function to reallocate host-pinned memory
   *
   * This function first deallocates any existing memory pointed to by `*pp`,
   * then allocates new host-pinned memory for `len` elements of type `T`.
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated host-pinned memory
   * @param[in] len Number of elements to allocate
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int reallocate_host(T **pp, const size_t len) {
    int error_code = COLVARS_OK;
    error_code |= deallocate_host(pp);
    error_code |= allocate_host_T((void **)pp, len, sizeof(T));
    return error_code;
  }
  /**
   * @brief Template function to allocate device memory asynchronously
   *
   * @tparam T The type of elements to allocate
   * @param[out] pp Pointer to the pointer that will hold the allocated device memory
   * @param[in] len Number of elements to allocate
   * @param[in] stream The CUDA stream to use for the allocation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int allocate_device_async(T **pp, const size_t len, cudaStream_t stream) {
    return allocate_device_T_async((void **)pp, len, sizeof(T), stream);
  }
  /**
   * @brief Template function to deallocate device memory
   *
   * @tparam T The type of elements to deallocate
   * @param[in,out] pp Pointer to the pointer that holds the allocated device memory
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int deallocate_device(T **pp) {
    return deallocate_device_T((void **)pp);
  }
  /**
   * @brief Template function to deallocate device memory asynchronously
   *
   * @tparam T The type of elements to deallocate
   * @param[in,out] pp Pointer to the pointer that holds the allocated device memory
   * @param[in] stream The CUDA stream to use for the deallocation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int deallocate_device_async(T **pp, cudaStream_t stream) {
    return deallocate_device_T_async((void **)pp, stream);
  }
  /**
   * @brief Template function to clear a device array to zero
   *
   * @tparam T The type of elements in the array
   * @param[in] data Pointer to the device array to clear
   * @param[in] ndata Number of elements in the array
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int clear_device_array(T *data, const size_t ndata) {
    return clear_device_array_T(data, ndata, sizeof(T));
  }
  /**
   * @brief Template function to clear a device array to zero asynchronously
   *
   * @tparam T The type of elements in the array
   * @param[in] data Pointer to the device array to clear
   * @param[in] ndata Number of elements in the array
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int clear_device_array_async(T *data, const size_t ndata, cudaStream_t stream) {
    return clear_device_array_T_async(data, ndata, sizeof(T), stream);
  }
  /**
   * @brief Template function to copy data from host to device
   *
   * @tparam T The type of elements to copy
   * @param[in] h_array Pointer to the host array
   * @param[in] d_array Pointer to the device array
   * @param[in] array_len Number of elements to copy
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_HtoD(const T *h_array, T *d_array, size_t array_len) {
    return copy_HtoD_T(h_array, d_array, array_len, sizeof(T));
  }
  /**
   * @brief Template function to copy data from host to device asynchronously
   *
   * @tparam T The type of elements to copy
   * @param[in] h_array Pointer to the host array
   * @param[out] d_array Pointer to the device array
   * @param[in] array_len Number of elements to copy
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_HtoD_async(const T *h_array, T *d_array, size_t array_len, cudaStream_t stream) {
    return copy_HtoD_T_async(h_array, d_array, array_len, sizeof(T), stream);
  }
  /**
   * @brief Template function to copy data from device to host
   *
   * @tparam T The type of elements to copy
   * @param[in] d_array Pointer to the device array
   * @param[out] h_array Pointer to the host array
   * @param[in] array_len Number of elements to copy
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoH(const T *d_array, T *h_array, size_t array_len) {
    return copy_DtoH_T(d_array, h_array, array_len, sizeof(T));
  }
  /**
   * @brief Template function to copy data from device to host asynchronously
   *
   * @tparam T The type of elements to copy
   * @param[in] d_array Pointer to the device array
   * @param[out] h_array Pointer to the host array
   * @param[in] array_len Number of elements to copy
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoH_async(const T *d_array, T *h_array, size_t array_len, cudaStream_t stream) {
    return copy_DtoH_T_async(d_array, h_array, array_len, sizeof(T), stream);
  }
  /**
   * @brief Template function to copy data from device to device
   *
   * @tparam T The type of elements to copy
   * @param[in] d_src Pointer to the source device array
   * @param[out] d_dst Pointer to the destination device array
   * @param[in] array_len Number of elements to copy
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoD(const T *d_src, T *d_dst, size_t array_len) {
    return copy_DtoD_T(d_src, d_dst, array_len, sizeof(T));
  }
  /**
   * @brief Template function to copy data from device to device asynchronously
   *
   * @tparam T The type of elements to copy
   * @param[in] d_src Pointer to the source device array
   * @param[out] d_dst Pointer to the destination device array
   * @param[in] array_len Number of elements to copy
   * @param[in] stream The CUDA stream to use for the operation
   * @return COLVARS_OK if successful, otherwise COLVARS_ERROR
   */
  template <typename T>
  int copy_DtoD_async(const T *d_src, T *d_dst, size_t array_len, cudaStream_t stream) {
    return copy_DtoD_T_async(d_src, d_dst, array_len, sizeof(T), stream);
  }
  /// @brief Memory management and data transfer implementations
  /// @{
  virtual int allocate_host_T(void **pp, const size_t len, const size_t sizeofT);
  virtual int deallocate_host_T(void **pp);
  virtual int allocate_device_T(void **pp, const size_t len, const size_t sizeofT);
  virtual int deallocate_device_T(void **pp);
  virtual int clear_device_array_T(void *data, const size_t ndata, const size_t sizeofT);
  virtual int allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, cudaStream_t stream);
  virtual int deallocate_device_T_async(void **pp, cudaStream_t stream);
  virtual int clear_device_array_T_async(void *data, const size_t ndata, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_HtoD_T(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT);
  virtual int copy_HtoD_T_async(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_DtoH_T(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT);
  virtual int copy_DtoH_T_async(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  virtual int copy_DtoD_T(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT);
  virtual int copy_DtoD_T_async(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT, cudaStream_t stream);
  /// @}
  /// @brief Functions to get device pointers for atom properties
  /// This functions should be overridden in derived proxy classes that manage actual GPU memory
  /// @{
  virtual float* proxy_atoms_masses_gpu_float() {return nullptr;}
  virtual float* proxy_atoms_charges_gpu_float() {return nullptr;}
  virtual cvm::real* proxy_atoms_masses_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_charges_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_positions_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_total_forces_gpu() {return nullptr;}
  virtual cvm::real* proxy_atoms_new_colvar_forces_gpu() {return nullptr;}
  /// @}
  /**
   * @brief This function will be called after atom groups are calculated on GPU.
   *
   * This function is useful when additional information is needed to transfer 
   * from the proxy. For example, the proxy can copy the lattice vectors in a 
   * separate stream, and this function can wait for that stream to complete.
   */
  virtual int wait_for_extra_info_ready();
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)
  /// \brief Destructor
  virtual ~colvarproxy_gpu();
protected:
  /// \brief Whether the proxy supports GPU
  bool support_gpu;
};

#endif // COLVARPROXY_GPU_H
