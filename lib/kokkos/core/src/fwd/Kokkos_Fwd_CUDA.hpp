// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CUDA_FWD_HPP_
#define KOKKOS_CUDA_FWD_HPP_
#if defined(KOKKOS_ENABLE_CUDA)
namespace Kokkos {

class CudaSpace;            ///< Memory space on Cuda GPU
class CudaUVMSpace;         ///< Memory space on Cuda GPU with UVM
class CudaHostPinnedSpace;  ///< Memory space on Host accessible to Cuda GPU
class Cuda;                 ///< Execution space for Cuda GPU

namespace Impl {

template <class ExecSpace>
void cuda_prefetch_pointer(const ExecSpace& /*space*/, const void* /*ptr*/,
                           size_t /*bytes*/, bool /*to_device*/) {}

void cuda_prefetch_pointer(const Cuda& space, const void* ptr, size_t bytes,
                           bool to_device);

}  // namespace Impl
}  // namespace Kokkos
#endif
#endif
