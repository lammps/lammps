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

#ifndef KOKKOS_OPENMPTARGET_INSTANCE_HPP
#define KOKKOS_OPENMPTARGET_INSTANCE_HPP

#include <Kokkos_Core.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

enum class openmp_fence_is_static { yes, no };

class OpenMPTargetInternal {
 private:
  OpenMPTargetInternal()                            = default;
  OpenMPTargetInternal(const OpenMPTargetInternal&) = default;
  OpenMPTargetInternal& operator=(const OpenMPTargetInternal&) = default;

 public:
  void fence(openmp_fence_is_static is_static = openmp_fence_is_static::no);
  void fence(const std::string& name,
             openmp_fence_is_static is_static = openmp_fence_is_static::no);

  /** \brief  Return the maximum amount of concurrency.  */
  int concurrency();

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose) const;

  static const char* name();

  //! Free any resources being consumed by the device.
  void impl_finalize();

  //! Has been initialized
  int impl_is_initialized();
  uint32_t impl_get_instance_id() const noexcept;
  //! Initialize, telling the CUDA run-time library which device to use.
  void impl_initialize();

  static OpenMPTargetInternal* impl_singleton();

 private:
  bool m_is_initialized  = false;
  uint32_t m_instance_id = Kokkos::Tools::Experimental::Impl::idForInstance<
      Kokkos::Experimental::OpenMPTarget>(reinterpret_cast<uintptr_t>(this));
};
}  // Namespace Impl
}  // Namespace Experimental
}  // Namespace Kokkos

#endif  // KOKKOS_OPENMPTARGET_INSTANCE_HPP
