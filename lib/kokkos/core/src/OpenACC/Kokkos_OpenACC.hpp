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

#ifndef KOKKOS_OPENACC_HPP
#define KOKKOS_OPENACC_HPP

#include <OpenACC/Kokkos_OpenACCSpace.hpp>
#include <Kokkos_Concepts.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <OpenACC/Kokkos_OpenACC_Traits.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>

#include <openacc.h>

#include <iosfwd>
#include <string>

// FIXME_OPENACC: Below macro is temporarily enabled to avoid issues on existing
// OpenACC compilers not supporting lambda with parallel loops.
// LLVM/Clacc compiler does not need this.
#ifndef KOKKOS_COMPILER_CLANG
#define KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS
#endif

namespace Kokkos::Experimental::Impl {
class OpenACCInternal;
}

namespace Kokkos::Experimental {

class OpenACC {
  Kokkos::Impl::HostSharedPtr<Impl::OpenACCInternal> m_space_instance;

  friend bool operator==(OpenACC const& lhs, OpenACC const& rhs) {
    return lhs.impl_internal_space_instance() ==
           rhs.impl_internal_space_instance();
  }
  friend bool operator!=(OpenACC const& lhs, OpenACC const& rhs) {
    return !(lhs == rhs);
  }

 public:
  using execution_space = OpenACC;
  using memory_space    = OpenACCSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using array_layout = LayoutLeft;
  using size_type    = memory_space::size_type;

  using scratch_memory_space = ScratchMemorySpace<OpenACC>;

  OpenACC();

  explicit OpenACC(int async_arg);

  static void impl_initialize(InitializationSettings const& settings);
  static void impl_finalize();
  static bool impl_is_initialized();

  void print_configuration(std::ostream& os, bool verbose = false) const;

  void fence(std::string const& name =
                 "Kokkos::OpenACC::fence(): Unnamed Instance Fence") const;
  static void impl_static_fence(std::string const& name);

  static char const* name() { return "OpenACC"; }
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static int concurrency() { return 256000; }  // FIXME_OPENACC
#else
  int concurrency() const { return 256000; }  // FIXME_OPENACC
#endif
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED static bool in_parallel() {
    return acc_on_device(acc_device_not_host);
  }
#endif
  uint32_t impl_instance_id() const noexcept;
  Impl::OpenACCInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }

  int acc_async_queue() const;
  int acc_device_number() const;
};

}  // namespace Kokkos::Experimental

template <>
struct Kokkos::Tools::Experimental::DeviceTypeTraits<
    ::Kokkos::Experimental::OpenACC> {
  static constexpr DeviceType id =
      ::Kokkos::Profiling::Experimental::DeviceType::OpenACC;
  static int device_id(const Kokkos::Experimental::OpenACC& accInstance) {
    return accInstance.acc_device_number();
  }
};

#endif
