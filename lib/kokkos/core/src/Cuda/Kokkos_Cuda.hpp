//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_CUDA_HPP
#define KOKKOS_CUDA_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_CUDA)

#include <Kokkos_Core_fwd.hpp>

#include <iosfwd>
#include <vector>

#include <impl/Kokkos_AnalyzePolicy.hpp>
#include <Cuda/Kokkos_CudaSpace.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>  // CUDA_SAFE_CALL

#include <Kokkos_Parallel.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class CudaInternal;
}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace Impl {
namespace Experimental {
enum class CudaLaunchMechanism : unsigned {
  Default        = 0,
  ConstantMemory = 1,
  GlobalMemory   = 2,
  LocalMemory    = 4
};

constexpr inline CudaLaunchMechanism operator|(CudaLaunchMechanism p1,
                                               CudaLaunchMechanism p2) {
  return static_cast<CudaLaunchMechanism>(static_cast<unsigned>(p1) |
                                          static_cast<unsigned>(p2));
}
constexpr inline CudaLaunchMechanism operator&(CudaLaunchMechanism p1,
                                               CudaLaunchMechanism p2) {
  return static_cast<CudaLaunchMechanism>(static_cast<unsigned>(p1) &
                                          static_cast<unsigned>(p2));
}

template <CudaLaunchMechanism l>
struct CudaDispatchProperties {
  CudaLaunchMechanism launch_mechanism = l;
};
}  // namespace Experimental

enum class ManageStream : bool { no, yes };

}  // namespace Impl
/// \class Cuda
/// \brief Kokkos Execution Space that uses CUDA to run on GPUs.
///
/// An "execution space" represents a parallel execution model.  It tells Kokkos
/// how to parallelize the execution of kernels in a parallel_for or
/// parallel_reduce.  For example, the Threads execution space uses
/// C++11 threads on a CPU, the OpenMP execution space uses the OpenMP language
/// extensions, and the Serial execution space executes "parallel" kernels
/// sequentially.  The Cuda execution space uses NVIDIA's CUDA programming
/// model to execute kernels in parallel on GPUs.
class Cuda {
 public:
  //! \name Type declarations that all Kokkos execution spaces must provide.
  //@{

  //! Tag this class as a kokkos execution space
  using execution_space = Cuda;

#if defined(KOKKOS_ENABLE_CUDA_UVM)
  //! This execution space's preferred memory space.
  using memory_space = CudaUVMSpace;
#else
  //! This execution space's preferred memory space.
  using memory_space = CudaSpace;
#endif

  //! This execution space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  //! The size_type best suited for this execution space.
  using size_type = memory_space::size_type;

  //! This execution space's preferred array layout.
  using array_layout = LayoutLeft;

  //!
  using scratch_memory_space = ScratchMemorySpace<Cuda>;

  //@}
  //--------------------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

  /// \brief True if and only if this method is being called in a
  ///   thread-parallel function.

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined(__CUDA_ARCH__)
    return true;
#else
    return false;
#endif
  }
#endif

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void impl_static_fence(const std::string& name);

  void fence(const std::string& name =
                 "Kokkos::Cuda::fence(): Unnamed Instance Fence") const;

  /** \brief  Return the maximum amount of concurrency.  */
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static int concurrency();
#else
  int concurrency() const;
#endif

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  //@}
  //--------------------------------------------------
  //! \name  Cuda space instances

  Cuda();

  explicit Cuda(cudaStream_t stream) : Cuda(stream, Impl::ManageStream::no) {}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  template <typename T = void>
  KOKKOS_DEPRECATED_WITH_COMMENT(
      "Cuda execution space should be constructed explicitly.")
  Cuda(cudaStream_t stream)
      : Cuda(stream) {}
#endif

  Cuda(cudaStream_t stream, Impl::ManageStream manage_stream);

  KOKKOS_DEPRECATED Cuda(cudaStream_t stream, bool manage_stream);

  //--------------------------------------------------------------------------
  //! Free any resources being consumed by the device.
  static void impl_finalize();

  //! Has been initialized
  static int impl_is_initialized();

  //! Initialize, telling the CUDA run-time library which device to use.
  static void impl_initialize(InitializationSettings const&);

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  /// \brief Cuda device architecture of the selected device.
  ///
  /// This matches the __CUDA_ARCH__ specification.
  KOKKOS_DEPRECATED static size_type device_arch() {
    const cudaDeviceProp cudaProp = Cuda().cuda_device_prop();
    return cudaProp.major * 100 + cudaProp.minor;
  }

  //! Query device count.
  KOKKOS_DEPRECATED static size_type detect_device_count() {
    int count;
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&count));
    return count;
  }

  /** \brief  Detect the available devices and their architecture
   *          as defined by the __CUDA_ARCH__ specification.
   */
  KOKKOS_DEPRECATED static std::vector<unsigned> detect_device_arch() {
    int count;
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&count));
    std::vector<unsigned> out;
    for (int i = 0; i < count; ++i) {
      cudaDeviceProp prop;
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceProperties(&prop, i));
      out.push_back(prop.major * 100 + prop.minor);
    }
    return out;
  }
#endif

  cudaStream_t cuda_stream() const;
  int cuda_device() const;
  const cudaDeviceProp& cuda_device_prop() const;

  //@}
  //--------------------------------------------------------------------------

  static const char* name();

  inline Impl::CudaInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }
  uint32_t impl_instance_id() const noexcept;

 private:
  friend bool operator==(Cuda const& lhs, Cuda const& rhs) {
    return lhs.impl_internal_space_instance() ==
           rhs.impl_internal_space_instance();
  }
  friend bool operator!=(Cuda const& lhs, Cuda const& rhs) {
    return !(lhs == rhs);
  }
  Kokkos::Impl::HostSharedPtr<Impl::CudaInternal> m_space_instance;
};

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Cuda> {
  /// \brief An ID to differentiate (for example) Serial from OpenMP in Tooling
  static constexpr DeviceType id = DeviceType::Cuda;
  static int device_id(const Cuda& exec) { return exec.cuda_device(); }
};
}  // namespace Experimental
}  // namespace Tools
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::CudaSpace,
                         Kokkos::Cuda::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

#if defined(KOKKOS_ENABLE_CUDA_UVM)

// If forcing use of UVM everywhere
// then must assume that CudaUVMSpace
// can be a stand-in for CudaSpace.
// This will fail when a strange host-side execution space
// that defines CudaUVMSpace as its preferredmemory space.

template <>
struct MemorySpaceAccess<Kokkos::CudaUVMSpace,
                         Kokkos::Cuda::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

#endif

}  // namespace Impl
}  // namespace Kokkos

#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDA_HPP */
