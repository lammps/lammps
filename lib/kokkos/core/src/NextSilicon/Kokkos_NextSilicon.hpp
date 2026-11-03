// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_HPP
#define KOKKOS_NEXTSILICON_HPP

#include <nextapi/intrinsics.h>

#include <NextSilicon/Kokkos_NextSilicon_Instance.hpp>
#include <NextSilicon/Kokkos_NextSiliconSpace.hpp>
#include <Kokkos_Concepts.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>

namespace Kokkos {
namespace Experimental::Impl {
class NextSiliconInternal;
}  // namespace Experimental::Impl

namespace Experimental {

class NextSilicon {
  Kokkos::Impl::HostSharedPtr<Impl::NextSiliconInternal> m_space_instance;

  friend bool operator==(NextSilicon const& lhs, NextSilicon const& rhs) {
    return lhs.impl_internal_space_instance() ==
           rhs.impl_internal_space_instance();
  }
  friend bool operator!=(NextSilicon const& lhs, NextSilicon const& rhs) {
    return !(lhs == rhs);
  }

 public:
  using execution_space = NextSilicon;
  using memory_space    = NextSiliconSharedSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using array_layout = LayoutLeft;
  using size_type    = memory_space::size_type;
  using index_type   = memory_space::index_type;

  using scratch_memory_space = ScratchMemorySpace<NextSilicon>;

  NextSilicon(const NextSilicon&) = default;
  NextSilicon(NextSilicon&& other) noexcept
      : NextSilicon(static_cast<const NextSilicon&>(other)) {}
  NextSilicon& operator=(const NextSilicon&) = default;
  NextSilicon& operator=(NextSilicon&& other) noexcept {
    return *this = static_cast<const NextSilicon&>(other);
  }
  NextSilicon();
  ~NextSilicon();

  static void impl_initialize(InitializationSettings const& settings);
  static void impl_finalize();

  void print_configuration(std::ostream& os, bool verbose = false) const;

  void fence(std::string const& name =
                 "Kokkos::NextSilicon::fence(): Unnamed Instance Fence") const;
  static void impl_static_fence(std::string const& name);

  static char const* name() { return "NextSilicon"; }
  static int concurrency() {
    return 64 * 1024; /* FIXME_NEXTSILICON - move to nextapi call */
  }
  uint32_t impl_instance_id() const noexcept;
  Impl::NextSiliconInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }
};

}  // namespace Experimental

namespace Impl {
template <>
struct MemorySpaceAccess<Experimental::NextSiliconSharedSpace,
                         Experimental::NextSilicon::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
};
}  // namespace Impl
}  // namespace Kokkos

template <>
struct Kokkos::Tools::Experimental::DeviceTypeTraits<
    ::Kokkos::Experimental::NextSilicon> {
  static constexpr DeviceType id =
      ::Kokkos::Profiling::Experimental::DeviceType::NextSilicon;

  static int device_id(const Kokkos::Experimental::NextSilicon& /*exec*/) {
    // FIXME_NEXTSILICON - return something meaningful here when supported by
    // NextAPI -- see NextSilicon ticket CS-611
    return 0;
  }
};

#endif  // KOKKOS_NEXTSILICON_HPP
