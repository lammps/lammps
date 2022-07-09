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

#ifndef KOKKOS_OPENMPTARGET_HPP
#define KOKKOS_OPENMPTARGET_HPP

#include <Kokkos_Core_fwd.hpp>

#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(_OPENMP)

#include <omp.h>

#include <cstddef>
#include <iosfwd>
#include <Kokkos_OpenMPTargetSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <KokkosExp_MDRangePolicy.hpp>
#include <impl/Kokkos_ExecSpaceInitializer.hpp>
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

  inline static bool in_parallel() { return omp_in_parallel(); }

  static void fence();
  static void fence(const std::string&);

  static void impl_static_fence();
  static void impl_static_fence(const std::string&);
  /** \brief  Return the maximum amount of concurrency.  */
  static int concurrency();

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream&, const bool detail = false);

  static const char* name();

  //! Free any resources being consumed by the device.
  void impl_finalize();

  //! Has been initialized
  static int impl_is_initialized();

  //! Initialize, telling the CUDA run-time library which device to use.
  void impl_initialize();

  inline Impl::OpenMPTargetInternal* impl_internal_space_instance() const {
    return m_space_instance;
  }

  OpenMPTarget();
  uint32_t impl_instance_id() const noexcept;

 private:
  Impl::OpenMPTargetInternal* m_space_instance;
};
}  // namespace Experimental

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

namespace Impl {

class OpenMPTargetSpaceInitializer : public ExecSpaceInitializerBase {
 public:
  OpenMPTargetSpaceInitializer()  = default;
  ~OpenMPTargetSpaceInitializer() = default;
  void initialize(const InitArguments& args) final;
  void finalize(const bool) final;
  void fence() final;
  void fence(const std::string&) final;
  void print_configuration(std::ostream& msg, const bool detail) final;
};

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#include <OpenMPTarget/Kokkos_OpenMPTarget_Exec.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel_MDRange.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Task.hpp>

/*--------------------------------------------------------------------------*/

#endif /* #if defined( KOKKOS_ENABLE_OPENMPTARGET ) && defined( _OPENMP ) */
#endif /* #ifndef KOKKOS_OPENMPTARGET_HPP */
