/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_3
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#else
KOKKOS_IMPL_WARNING("Including non-public Kokkos header files is not allowed.")
#endif
#endif
#ifndef KOKKOS_CUDA_HPP
#define KOKKOS_CUDA_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_CUDA)

#include <Kokkos_Core_fwd.hpp>

#include <iosfwd>
#include <vector>

#include <impl/Kokkos_AnalyzePolicy.hpp>
#include <Kokkos_CudaSpace.hpp>
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
class CudaExec;
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
  KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined(__CUDA_ARCH__)
    return true;
#else
    return false;
#endif
  }

  /** \brief  Set the device in a "sleep" state.
   *
   * This function sets the device in a "sleep" state in which it is
   * not ready for work.  This may consume less resources than if the
   * device were in an "awake" state, but it may also take time to
   * bring the device from a sleep state to be ready for work.
   *
   * \return True if the device is in the "sleep" state, else false if
   *   the device is actively working and could not enter the "sleep"
   *   state.
   */
  static bool sleep();

  /// \brief Wake the device from the 'sleep' state so it is ready for work.
  ///
  /// \return True if the device is in the "ready" state, else "false"
  ///  if the device is actively working (which also means that it's
  ///  awake).
  static bool wake();

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
  static int concurrency();

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  //@}
  //--------------------------------------------------
  //! \name  Cuda space instances

  Cuda();

  Cuda(cudaStream_t stream, bool manage_stream = false);

  //--------------------------------------------------------------------------
  //! Free any resources being consumed by the device.
  static void impl_finalize();

  //! Has been initialized
  static int impl_is_initialized();

  //! Initialize, telling the CUDA run-time library which device to use.
  static void impl_initialize(InitializationSettings const&);

  /// \brief Cuda device architecture of the selected device.
  ///
  /// This matches the __CUDA_ARCH__ specification.
  static size_type device_arch();

  //! Query device count.
  static size_type detect_device_count();

  /** \brief  Detect the available devices and their architecture
   *          as defined by the __CUDA_ARCH__ specification.
   */
  static std::vector<unsigned> detect_device_arch();

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

namespace Impl {

template <class DT, class... DP>
struct ZeroMemset<Kokkos::Cuda, DT, DP...> {
  ZeroMemset(const Kokkos::Cuda& exec_space_instance,
             const View<DT, DP...>& dst,
             typename View<DT, DP...>::const_value_type&) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaMemsetAsync(
        dst.data(), 0,
        dst.size() * sizeof(typename View<DT, DP...>::value_type),
        exec_space_instance.cuda_stream()));
  }

  ZeroMemset(const View<DT, DP...>& dst,
             typename View<DT, DP...>::const_value_type&) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        cudaMemset(dst.data(), 0,
                   dst.size() * sizeof(typename View<DT, DP...>::value_type)));
  }
};
}  // namespace Impl
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
