#include "colvarproxy_gpu.h"
#include "colvarmodule.h"

using namespace colvars_gpu;

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP) || defined (COLVARS_SYCL)
int colvarproxy_gpu::allocate_host_T(void **pp, const size_t len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaHostAlloc(pp, sizeofT*len, cudaHostAllocDefault));
  return error_code;
}

int colvarproxy_gpu::deallocate_host_T(void **pp) {
  int error_code = COLVARS_OK;
  if (*pp != nullptr) {
    error_code |= checkGPUError(cudaFreeHost((void *)(*pp)));
    *pp = nullptr;
  }
  return error_code;
}

int colvarproxy_gpu::allocate_device_T(void **pp, const size_t len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMalloc(pp, sizeofT*len));
  return error_code;
}

int colvarproxy_gpu::deallocate_device_T(void **pp) {
  int error_code = COLVARS_OK;
  if (*pp != nullptr) {
    error_code |= checkGPUError(cudaFree((void *)(*pp)));
    *pp = nullptr;
  }
  return error_code;
}

int colvarproxy_gpu::allocate_device_T_async(void **pp, const size_t len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMallocAsync(pp, sizeofT*len, stream));
  return error_code;
}

int colvarproxy_gpu::deallocate_device_T_async(void **pp, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  if (*pp != nullptr) {
    error_code |= checkGPUError(cudaFreeAsync((void *)(*pp), stream));
    *pp = nullptr;
  }
  return error_code;
}

int colvarproxy_gpu::clear_device_array_T(void *data, const size_t ndata, const size_t sizeofT) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemset(data, 0, sizeofT*ndata));
  return error_code;
}

int colvarproxy_gpu::clear_device_array_T_async(void *data, const size_t ndata, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemsetAsync(data, 0, sizeofT*ndata, stream));
  return error_code;
}

int colvarproxy_gpu::copy_HtoD_T(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemcpy(d_array, h_array, sizeofT*array_len, cudaMemcpyHostToDevice));
  return error_code;
}

int colvarproxy_gpu::copy_HtoD_T_async(const void *h_array, void *d_array, size_t array_len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemcpyAsync(d_array, h_array, sizeofT*array_len, cudaMemcpyHostToDevice, stream));
  return error_code;
}

int colvarproxy_gpu::copy_DtoH_T(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemcpy(h_array, d_array, sizeofT*array_len, cudaMemcpyDeviceToHost));
  return error_code;
}

int colvarproxy_gpu::copy_DtoH_T_async(const void *d_array, void *h_array, size_t array_len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemcpyAsync(h_array, d_array, sizeofT*array_len, cudaMemcpyDeviceToHost, stream));
  return error_code;
}

int colvarproxy_gpu::copy_DtoD_T(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemcpy(d_dst, d_src, sizeofT*array_len, cudaMemcpyDeviceToDevice));
  return error_code;
}

int colvarproxy_gpu::copy_DtoD_T_async(const void *d_src, void *d_dst, size_t array_len, const size_t sizeofT, cudaStream_t stream) {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaMemcpyAsync(d_dst, d_src, sizeofT*array_len, cudaMemcpyDeviceToDevice, stream));
  return error_code;
}

int colvarproxy_gpu::wait_for_extra_info_ready() {
  return COLVARS_OK;
}

#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)

colvarproxy_gpu::~colvarproxy_gpu() {
}
