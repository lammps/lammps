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

#ifndef KOKKOS_ROCM_HPP
#define KOKKOS_ROCM_HPP

#include <Kokkos_Core_fwd.hpp>

#if defined(KOKKOS_ENABLE_ROCM)

class dim3 {
 public:
  int x, y, z;
  dim3(int _x, int _y, int _z) : x(_x), y(_y), z(_z){};
};

#include <ROCm/hc_math_std.hpp>
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <cstddef>
#include <iosfwd>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ROCmSpace.hpp>
#include <ROCm/Kokkos_ROCm_Exec.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

#include <hc.hpp>
#include <hc_am.hpp>
#include <amp_math.h>

#if defined(__HCC_ACCELERATOR__)

using namespace ::Concurrency::precise_math;

#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class ROCmExec;
}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
/// \class ROCm
/// \brief Kokkos device for multicore processors in the host memory space.
class ROCm {
 public:
  //------------------------------------
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! Tag this class as a kokkos execution space
  typedef ROCm execution_space;
  typedef ROCmSpace memory_space;
  typedef Kokkos::Device<execution_space, memory_space> device_type;

  typedef LayoutLeft array_layout;
  typedef HostSpace::size_type size_type;

  typedef ScratchMemorySpace<ROCm> scratch_memory_space;

  ~ROCm() {}
  ROCm();
  //  explicit ROCm( const int instance_id );

  ROCm(ROCm&&)      = default;
  ROCm(const ROCm&) = default;
  ROCm& operator=(ROCm&&) = default;
  ROCm& operator=(const ROCm&) = default;

  //@}
  //------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

  KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined(__HCC_ACCELERATOR__)
    return true;
#else
    return false;
#endif
  }

  /** \brief  Set the device in a "sleep" state. */
  static bool sleep();

  /** \brief Wake the device from the 'sleep' state. A noop for OpenMP. */
  static bool wake();

  /** \brief Wait until all dispatched functors complete. A noop for OpenMP. */
  static void impl_static_fence();

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  static void fence();
#else
  void fence() const;
#endif

  /// \brief Print configuration information to the given output stream.
  static void print_configuration(std::ostream&, const bool detail = false);

  /// \brief Free any resources being consumed by the device.
  static void finalize();

  /** \brief  Initialize the device.
   *
   */
  struct SelectDevice {
    int rocm_device_id;
    SelectDevice() : rocm_device_id(1) {}
    explicit SelectDevice(int id) : rocm_device_id(id + 1) {}
  };

  int rocm_device() const { return m_device; }
  bool isAPU();
  bool isAPU(int device);

  static void initialize(const SelectDevice = SelectDevice());

  static int is_initialized();

  //  static size_type device_arch();

  //  static size_type detect_device_count();

  static int concurrency();
  static const char* name();

 private:
  int m_device;
};
}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::Experimental::ROCmSpace,
                         Kokkos::Experimental::ROCm::scratch_memory_space> {
  enum { assignable = false };
  enum { accessible = true };
  enum { deepcopy = false };
};

template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::Experimental::ROCm::memory_space,
    Kokkos::Experimental::ROCm::scratch_memory_space> {
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void*) {}
};

template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::HostSpace, Kokkos::Experimental::ROCm::scratch_memory_space> {
  enum { value = false };
  inline static void verify(void) {
    Kokkos::Experimental::ROCmSpace::access_error();
  }
  inline static void verify(const void* p) {
    Kokkos::Experimental::ROCmSpace::access_error(p);
  }
};

}  // namespace Impl
}  // namespace Kokkos

#define threadIdx_x (hc_get_workitem_id(0))
#define threadIdx_y (hc_get_workitem_id(1))
#define threadIdx_z (hc_get_workitem_id(2))

#define blockIdx_x (hc_get_group_id(0))
#define blockIdx_y (hc_get_group_id(1))
#define blockIdx_z (hc_get_group_id(2))

#define blockDim_x (hc_get_group_size(0))
#define blockDim_y (hc_get_group_size(1))
#define blockDim_z (hc_get_group_size(2))

#define gridDim_x (hc_get_num_groups(0))
#define gridDim_y (hc_get_num_groups(1))
#define gridDim_z (hc_get_num_groups(2))

#include <ROCm/Kokkos_ROCm_Parallel.hpp>
#include <ROCm/Kokkos_ROCm_Task.hpp>

#endif
#endif
