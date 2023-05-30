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
#ifndef KOKKOS_OPENMP_HPP
#define KOKKOS_OPENMP_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_OPENMP)

#include <Kokkos_Core_fwd.hpp>

#include <cstddef>
#include <iosfwd>
#include <Kokkos_HostSpace.hpp>

#ifdef KOKKOS_ENABLE_HBWSPACE
#include <Kokkos_HBWSpace.hpp>
#endif

#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>

#include <vector>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace Impl {
class OpenMPInternal;
}

/// \class OpenMP
/// \brief Kokkos device for multicore processors in the host memory space.
class OpenMP {
 public:
  //! Tag this class as a kokkos execution space
  using execution_space = OpenMP;

  using memory_space =
#ifdef KOKKOS_ENABLE_HBWSPACE
      Experimental::HBWSpace;
#else
      HostSpace;
#endif

  //! This execution space preferred device_type
  using device_type          = Kokkos::Device<execution_space, memory_space>;
  using array_layout         = LayoutRight;
  using size_type            = memory_space::size_type;
  using scratch_memory_space = ScratchMemorySpace<OpenMP>;

  OpenMP();

  /// \brief Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  /// \brief is the instance running a parallel algorithm
  inline static bool in_parallel(OpenMP const& = OpenMP()) noexcept;

  /// \brief Wait until all dispatched functors complete on the given instance
  ///
  ///  This is a no-op on OpenMP
  static void impl_static_fence(std::string const& name);

  void fence(std::string const& name =
                 "Kokkos::OpenMP::fence: Unnamed Instance Fence") const;

  /// \brief Does the given instance return immediately after launching
  /// a parallel algorithm
  ///
  /// This always returns false on OpenMP
  inline static bool is_asynchronous(OpenMP const& = OpenMP()) noexcept;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
  /// \brief Partition the default instance into new instances without creating
  ///  new masters
  ///
  /// This is a no-op on OpenMP since the default instance cannot be partitioned
  /// without promoting other threads to 'master'
  static std::vector<OpenMP> partition(...);

  /// Non-default instances should be ref-counted so that when the last
  /// is destroyed the instance resources are released
  ///
  /// This is a no-op on OpenMP since a non default instance cannot be created
  static OpenMP create_instance(...);

  /// \brief Partition the default instance and call 'f' on each new 'master'
  /// thread
  ///
  /// Func is a functor with the following signiture
  ///   void( int partition_id, int num_partitions )
  template <typename F>
  KOKKOS_DEPRECATED static void partition_master(
      F const& f, int requested_num_partitions = 0,
      int requested_partition_size = 0);
#endif

  // use UniqueToken
  static int concurrency();

  static void impl_initialize(InitializationSettings const&);

  /// \brief is the default execution space initialized for current 'master'
  /// thread
  static bool impl_is_initialized() noexcept;

  /// \brief Free any resources being consumed by the default execution space
  static void impl_finalize();

  inline static int impl_thread_pool_size() noexcept;

  /** \brief  The rank of the executing thread in this thread pool */
  KOKKOS_INLINE_FUNCTION
  static int impl_thread_pool_rank() noexcept;

  inline static int impl_thread_pool_size(int depth);

  // use UniqueToken
  inline static int impl_max_hardware_threads() noexcept;

  // use UniqueToken
  KOKKOS_INLINE_FUNCTION
  static int impl_hardware_thread_id() noexcept;

  static int impl_get_current_max_threads() noexcept;

  Impl::OpenMPInternal* impl_internal_space_instance() const {
#ifdef KOKKOS_IMPL_WORKAROUND_ICE_IN_TRILINOS_WITH_OLD_INTEL_COMPILERS
    return m_space_instance;
#else
    return m_space_instance.get();
#endif
  }

  static constexpr const char* name() noexcept { return "OpenMP"; }
  uint32_t impl_instance_id() const noexcept { return 1; }

 private:
#ifdef KOKKOS_IMPL_WORKAROUND_ICE_IN_TRILINOS_WITH_OLD_INTEL_COMPILERS
  Impl::OpenMPInternal* m_space_instance;
#else
  Kokkos::Impl::HostSharedPtr<Impl::OpenMPInternal> m_space_instance;
#endif
};

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<OpenMP> {
  static constexpr DeviceType id = DeviceType::OpenMP;
  static int device_id(const OpenMP&) { return 0; }
};
}  // namespace Experimental
}  // namespace Tools
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::OpenMP::memory_space,
                         Kokkos::OpenMP::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#include <OpenMP/Kokkos_OpenMP_Instance.hpp>
#include <OpenMP/Kokkos_OpenMP_Team.hpp>
#include <OpenMP/Kokkos_OpenMP_Parallel.hpp>
#include <OpenMP/Kokkos_OpenMP_Task.hpp>

#include <KokkosExp_MDRangePolicy.hpp>
/*--------------------------------------------------------------------------*/

#endif /* #if defined( KOKKOS_ENABLE_OPENMP ) && defined( _OPENMP ) */
#endif /* #ifndef KOKKOS_OPENMP_HPP */
