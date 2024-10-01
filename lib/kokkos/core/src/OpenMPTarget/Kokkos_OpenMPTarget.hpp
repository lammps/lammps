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
#ifndef KOKKOS_OPENMPTARGET_HPP
#define KOKKOS_OPENMPTARGET_HPP

#include <Kokkos_Core_fwd.hpp>

#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(_OPENMP)

#include <omp.h>

#include <cstddef>
#include <iosfwd>
#include <OpenMPTarget/Kokkos_OpenMPTargetSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>
#include <KokkosExp_MDRangePolicy.hpp>
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
namespace Impl {
class OpenMPTargetInternal;
}

/// \class OpenMPTarget
/// \brief Kokkos device for multicore processors in the host memory space.
class OpenMPTarget {
 public:
  //------------------------------------
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! Tag this class as a kokkos execution space
  using execution_space = OpenMPTarget;
  using memory_space    = OpenMPTargetSpace;
  //! This execution space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  using array_layout = LayoutLeft;
  using size_type    = memory_space::size_type;

  using scratch_memory_space = ScratchMemorySpace<OpenMPTarget>;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED inline static bool in_parallel() {
    return omp_in_parallel();
  }
#endif

  static void fence(const std::string& name =
                        "Kokkos::OpenMPTarget::fence: Unnamed Instance Fence");

  static void impl_static_fence(const std::string& name);

  /** \brief  Return the maximum amount of concurrency.  */
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static int concurrency();
#else
  int concurrency() const;
#endif

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  static const char* name();

  //! Free any resources being consumed by the device.
  static void impl_finalize();

  //! Has been initialized
  static int impl_is_initialized();

  //! Initialize, telling the CUDA run-time library which device to use.
  static void impl_initialize(InitializationSettings const&);

  inline Impl::OpenMPTargetInternal* impl_internal_space_instance() const {
    return m_space_instance;
  }

  OpenMPTarget();
  uint32_t impl_instance_id() const noexcept;

 private:
  friend bool operator==(OpenMPTarget const& lhs, OpenMPTarget const& rhs) {
    return lhs.impl_internal_space_instance() ==
           rhs.impl_internal_space_instance();
  }
  friend bool operator!=(OpenMPTarget const& lhs, OpenMPTarget const& rhs) {
    return !(lhs == rhs);
  }
  Impl::OpenMPTargetInternal* m_space_instance;
};
}  // namespace Experimental

namespace Impl {
template <>
struct MemorySpaceAccess<
    Kokkos::Experimental::OpenMPTargetSpace,
    Kokkos::Experimental::OpenMPTarget::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};
}  // namespace Impl

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<::Kokkos::Experimental::OpenMPTarget> {
  static constexpr DeviceType id =
      ::Kokkos::Profiling::Experimental::DeviceType::OpenMPTarget;
  static int device_id(const Kokkos::Experimental::OpenMPTarget&) {
    return omp_get_default_device();
  }
};
}  // namespace Experimental
}  // namespace Tools

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelFor_MDRange.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_ParallelReduce_MDRange.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Task.hpp>

/*--------------------------------------------------------------------------*/

#endif /* #if defined( KOKKOS_ENABLE_OPENMPTARGET ) && defined( _OPENMP ) */
#endif /* #ifndef KOKKOS_OPENMPTARGET_HPP */
