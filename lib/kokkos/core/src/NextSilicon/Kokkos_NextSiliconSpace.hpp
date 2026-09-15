// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_NEXTSILICON_SPACE_HPP
#define KOKKOS_NEXTSILICON_SPACE_HPP

#include <Kokkos_Concepts.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <iosfwd>
#include <type_traits>

namespace Kokkos::Experimental {

class NextSilicon;

class NextSiliconSharedSpace {
 public:
  using memory_space    = NextSiliconSharedSpace;
  using execution_space = NextSilicon;
  using device_type     = Kokkos::Device<execution_space, memory_space>;
  using size_type       = size_t;
  using index_type      = std::make_signed_t<size_type>;

  NextSiliconSharedSpace() = default;

  /**\brief  Allocate untracked memory in the space */
  template <typename ExecutionSpace>
  void* allocate(const ExecutionSpace& exec_space,
                 const size_t arg_alloc_size) const {
    return allocate(exec_space, "[unlabeled]", arg_alloc_size);
  }
  template <typename ExecutionSpace>
  void* allocate(const ExecutionSpace& exec_space, const char* arg_label,
                 const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const {
    return impl_allocate(exec_space, arg_label, arg_alloc_size,
                         arg_logical_size);
  }
  void* allocate(const size_t arg_alloc_size) const;
  void* allocate(const char* arg_label, const size_t arg_alloc_size,
                 const size_t arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* const arg_alloc_ptr, const size_t arg_alloc_size) const;
  void deallocate(const char* arg_label, void* const arg_alloc_ptr,
                  const size_t arg_alloc_size,
                  const size_t arg_logical_size = 0) const;

  static constexpr char const* name() { return "NextSiliconSharedSpace"; }

 private:
  template <typename ExecutionSpace>
  void* impl_allocate(const ExecutionSpace& exec_space, const char* arg_label,
                      const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      const Kokkos::Tools::SpaceHandle arg_handle =
                          Kokkos::Tools::make_space_handle(name())) const {
    if (!std::is_same_v<ExecutionSpace, Kokkos::Experimental::NextSilicon>) {
      exec_space.fence();
    }
    return impl_allocate(arg_label, arg_alloc_size, arg_logical_size,
                         arg_handle);
  }
  void* impl_allocate(const char* arg_label, const size_t arg_alloc_size,
                      const size_t arg_logical_size = 0,
                      const Kokkos::Tools::SpaceHandle =
                          Kokkos::Tools::make_space_handle(name())) const;
  void impl_deallocate(const char* arg_label, void* const arg_alloc_ptr,
                       const size_t arg_alloc_size,
                       const size_t arg_logical_size = 0,
                       const Kokkos::Tools::SpaceHandle =
                           Kokkos::Tools::make_space_handle(name())) const;
};

}  // namespace Kokkos::Experimental

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

static_assert(Kokkos::Impl::MemorySpaceAccess<
              Experimental::NextSiliconSharedSpace,
              Experimental::NextSiliconSharedSpace>::assignable);

template <>
struct MemorySpaceAccess<Kokkos::Experimental::NextSiliconSharedSpace,
                         Kokkos::HostSpace> {
  enum : bool { assignable = false };
  enum : bool {
    accessible = false
  };  // This is a lie, because NextSilicon can access HostSpace. Set this way
      // so the HostMirror of NextSiliconSharedSpace uses the host execution
      // space.
};

template <>
struct MemorySpaceAccess<Kokkos::HostSpace,
                         Kokkos::Experimental::NextSiliconSharedSpace> {
  enum : bool { assignable = false };
  enum : bool {
    accessible = true
  };  // Unlike <NextSiliconSharedSpace, HostSpace>, this is true so that
      // HostMirror uses NextSiliconSharedSpace as the memory space
};

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

#endif
