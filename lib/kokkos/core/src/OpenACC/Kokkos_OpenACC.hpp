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

#ifndef KOKKOS_OPENACC_HPP
#define KOKKOS_OPENACC_HPP

#include <OpenACC/Kokkos_OpenACCSpace.hpp>
#include <Kokkos_Concepts.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <OpenACC/Kokkos_OpenACC_Traits.hpp>

#include <openacc.h>

#include <iosfwd>
#include <string>

namespace Kokkos::Experimental::Impl {
class OpenACCInternal;
}

namespace Kokkos::Experimental {

class OpenACC {
  Impl::OpenACCInternal* m_space_instance = nullptr;

 public:
  using execution_space = OpenACC;
  using memory_space    = OpenACCSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using array_layout = LayoutLeft;
  using size_type    = memory_space::size_type;

  using scratch_memory_space = ScratchMemorySpace<OpenACC>;

  OpenACC();

  static void impl_initialize(InitializationSettings const& settings);
  static void impl_finalize();
  static bool impl_is_initialized();

  void print_configuration(std::ostream& os, bool verbose = false) const;

  void fence(std::string const& name =
                 "Kokkos::OpenACC::fence(): Unnamed Instance Fence") const;
  static void impl_static_fence(std::string const& name);

  static char const* name() { return "OpenACC"; }
  static int concurrency() { return 256000; }  // FIXME_OPENACC
  static bool in_parallel() { return acc_on_device(acc_device_not_host); }
  uint32_t impl_instance_id() const noexcept;
};

}  // namespace Kokkos::Experimental

template <>
struct Kokkos::Tools::Experimental::DeviceTypeTraits<
    ::Kokkos::Experimental::OpenACC> {
  static constexpr DeviceType id =
      ::Kokkos::Profiling::Experimental::DeviceType::OpenACC;
  // FIXME_OPENACC: Need to return the device id from the execution space
  // instance. In fact, acc_get_device_num() will return the same value as the
  // device id from the execution space instance except for the host fallback
  // case, where the device id may need to be updated with the value of
  // acc_get_device_num().
  static int device_id(const Kokkos::Experimental::OpenACC&) {
    using Kokkos::Experimental::Impl::OpenACC_Traits;
    return acc_get_device_num(OpenACC_Traits::dev_type);
  }
};

#endif
