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

#ifndef KOKKOS_VIEW_ALLOC_HPP
#define KOKKOS_VIEW_ALLOC_HPP

#include <cstring>
#include <type_traits>
#include <string>

#include <impl/Kokkos_Tools.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

namespace Kokkos::Impl {

template <typename T>
bool is_zero_byte(const T& x) {
  constexpr std::byte all_zeroes[sizeof(T)] = {};
  return std::memcmp(&x, all_zeroes, sizeof(T)) == 0;
}

//----------------------------------------------------------------------------

/*
 *  The construction, assignment to default, and destruction
 *  are merged into a single functor.
 *  Primarily to work around an unresolved CUDA back-end bug
 *  that would lose the destruction cuda device function when
 *  called from the shared memory tracking destruction.
 *  Secondarily to have two fewer partial specializations.
 */
template <class DeviceType, class ValueType,
          bool IsScalar = std::is_scalar<ValueType>::value>
struct ViewValueFunctor;

template <class DeviceType, class ValueType>
struct ViewValueFunctor<DeviceType, ValueType, false /* is_scalar */> {
  using ExecSpace = typename DeviceType::execution_space;

  struct DestroyTag {};
  struct ConstructTag {};

  ExecSpace space;
  ValueType* ptr;
  size_t n;
  std::string name;
  bool default_exec_space;

  template <class _ValueType = ValueType>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<std::is_default_constructible<_ValueType>::value>
      operator()(ConstructTag const&, const size_t i) const {
    new (ptr + i) ValueType();
  }

  KOKKOS_INLINE_FUNCTION void operator()(DestroyTag const&,
                                         const size_t i) const {
    (ptr + i)->~ValueType();
  }

  ViewValueFunctor()                        = default;
  ViewValueFunctor(const ViewValueFunctor&) = default;
  ViewValueFunctor& operator=(const ViewValueFunctor&) = default;

  ViewValueFunctor(ExecSpace const& arg_space, ValueType* const arg_ptr,
                   size_t const arg_n, std::string arg_name)
      : space(arg_space),
        ptr(arg_ptr),
        n(arg_n),
        name(std::move(arg_name)),
        default_exec_space(false) {
    functor_instantiate_workaround();
  }

  ViewValueFunctor(ValueType* const arg_ptr, size_t const arg_n,
                   std::string arg_name)
      : space(ExecSpace{}),
        ptr(arg_ptr),
        n(arg_n),
        name(std::move(arg_name)),
        default_exec_space(true) {
    functor_instantiate_workaround();
  }

  template <typename Dummy = ValueType>
  std::enable_if_t<std::is_trivial<Dummy>::value &&
                   std::is_trivially_copy_assignable<ValueType>::value>
  construct_dispatch() {
    ValueType value{};
// On A64FX memset seems to do the wrong thing with regards to first touch
// leading to the significant performance issues
#ifndef KOKKOS_ARCH_A64FX
    if (Impl::is_zero_byte(value)) {
      uint64_t kpID = 0;
      if (Kokkos::Profiling::profileLibraryLoaded()) {
        // We are not really using parallel_for here but using beginParallelFor
        // instead of begin_parallel_for (and adding "via memset") is the best
        // we can do to indicate that this is not supposed to be tunable (and
        // doesn't really execute a parallel_for).
        Kokkos::Profiling::beginParallelFor(
            "Kokkos::View::initialization [" + name + "] via memset",
            Kokkos::Profiling::Experimental::device_id(space), &kpID);
      }
      (void)ZeroMemset(
          space, Kokkos::View<ValueType*, typename DeviceType::memory_space,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(ptr, n));

      if (Kokkos::Profiling::profileLibraryLoaded()) {
        Kokkos::Profiling::endParallelFor(kpID);
      }
      if (default_exec_space)
        space.fence("Kokkos::Impl::ViewValueFunctor: View init/destroy fence");
    } else {
#endif
      parallel_for_implementation<ConstructTag>();
#ifndef KOKKOS_ARCH_A64FX
    }
#endif
  }

  template <typename Dummy = ValueType>
  std::enable_if_t<!(std::is_trivial<Dummy>::value &&
                     std::is_trivially_copy_assignable<ValueType>::value)>
  construct_dispatch() {
    parallel_for_implementation<ConstructTag>();
  }

  template <typename Tag>
  void parallel_for_implementation() {
    using PolicyType =
        Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int64_t>, Tag>;
    PolicyType policy(space, 0, n);
    uint64_t kpID = 0;
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      const std::string functor_name =
          (std::is_same_v<Tag, DestroyTag>
               ? "Kokkos::View::destruction [" + name + "]"
               : "Kokkos::View::initialization [" + name + "]");
      Kokkos::Profiling::beginParallelFor(
          functor_name, Kokkos::Profiling::Experimental::device_id(space),
          &kpID);
    }

#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same<ExecSpace, Kokkos::Cuda>::value) {
      Kokkos::Impl::cuda_prefetch_pointer(space, ptr, sizeof(ValueType) * n,
                                          true);
    }
#endif
    const Kokkos::Impl::ParallelFor<ViewValueFunctor, PolicyType> closure(
        *this, policy);
    closure.execute();
    if (default_exec_space || std::is_same_v<Tag, DestroyTag>)
      space.fence("Kokkos::Impl::ViewValueFunctor: View init/destroy fence");
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endParallelFor(kpID);
    }
  }

  void construct_shared_allocation() { construct_dispatch(); }

  void destroy_shared_allocation() {
#ifdef KOKKOS_ENABLE_IMPL_VIEW_OF_VIEWS_DESTRUCTOR_PRECONDITION_VIOLATION_WORKAROUND
    if constexpr (std::is_same_v<typename ExecSpace::memory_space,
                                 Kokkos::HostSpace>)
      for (size_t i = 0; i < n; ++i) (ptr + i)->~ValueType();
    else
#endif
    {
      parallel_for_implementation<DestroyTag>();
    }
  }

  // This function is to ensure that the functor with DestroyTag is instantiated
  // This is a workaround to avoid "cudaErrorInvalidDeviceFunction" error later
  // when the function is queried with cudaFuncGetAttributes
  void functor_instantiate_workaround() {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET)
    if (false) {
      parallel_for_implementation<DestroyTag>();
    }
#endif
  }
};

template <class DeviceType, class ValueType>
struct ViewValueFunctor<DeviceType, ValueType, true /* is_scalar */> {
  using ExecSpace  = typename DeviceType::execution_space;
  using PolicyType = Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int64_t>>;

  ExecSpace space;
  ValueType* ptr;
  size_t n;
  std::string name;
  bool default_exec_space;

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t i) const { ptr[i] = ValueType(); }

  ViewValueFunctor()                        = default;
  ViewValueFunctor(const ViewValueFunctor&) = default;
  ViewValueFunctor& operator=(const ViewValueFunctor&) = default;

  ViewValueFunctor(ExecSpace const& arg_space, ValueType* const arg_ptr,
                   size_t const arg_n, std::string arg_name)
      : space(arg_space),
        ptr(arg_ptr),
        n(arg_n),
        name(std::move(arg_name)),
        default_exec_space(false) {}

  ViewValueFunctor(ValueType* const arg_ptr, size_t const arg_n,
                   std::string arg_name)
      : space(ExecSpace{}),
        ptr(arg_ptr),
        n(arg_n),
        name(std::move(arg_name)),
        default_exec_space(true) {}

  template <typename Dummy = ValueType>
  std::enable_if_t<std::is_trivial<Dummy>::value &&
                   std::is_trivially_copy_assignable<Dummy>::value>
  construct_shared_allocation() {
    // Shortcut for zero initialization
// On A64FX memset seems to do the wrong thing with regards to first touch
// leading to the significant performance issues
#ifndef KOKKOS_ARCH_A64FX
    ValueType value{};
    if (Impl::is_zero_byte(value)) {
      uint64_t kpID = 0;
      if (Kokkos::Profiling::profileLibraryLoaded()) {
        // We are not really using parallel_for here but using beginParallelFor
        // instead of begin_parallel_for (and adding "via memset") is the best
        // we can do to indicate that this is not supposed to be tunable (and
        // doesn't really execute a parallel_for).
        Kokkos::Profiling::beginParallelFor(
            "Kokkos::View::initialization [" + name + "] via memset",
            Kokkos::Profiling::Experimental::device_id(space), &kpID);
      }

      (void)ZeroMemset(
          space, Kokkos::View<ValueType*, typename DeviceType::memory_space,
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>>(ptr, n));

      if (Kokkos::Profiling::profileLibraryLoaded()) {
        Kokkos::Profiling::endParallelFor(kpID);
      }
      if (default_exec_space)
        space.fence("Kokkos::Impl::ViewValueFunctor: View init/destroy fence");
    } else {
#endif
      parallel_for_implementation();
#ifndef KOKKOS_ARCH_A64FX
    }
#endif
  }

  template <typename Dummy = ValueType>
  std::enable_if_t<!(std::is_trivial<Dummy>::value &&
                     std::is_trivially_copy_assignable<Dummy>::value)>
  construct_shared_allocation() {
    parallel_for_implementation();
  }

  void parallel_for_implementation() {
    PolicyType policy(space, 0, n);
    uint64_t kpID = 0;
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::beginParallelFor(
          "Kokkos::View::initialization [" + name + "]",
          Kokkos::Profiling::Experimental::device_id(space), &kpID);
    }
#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same<ExecSpace, Kokkos::Cuda>::value) {
      Kokkos::Impl::cuda_prefetch_pointer(space, ptr, sizeof(ValueType) * n,
                                          true);
    }
#endif
    const Kokkos::Impl::ParallelFor<ViewValueFunctor, PolicyType> closure(
        *this, policy);
    closure.execute();
    if (default_exec_space)
      space.fence(
          "Kokkos::Impl::ViewValueFunctor: Fence after setting values in "
          "view");
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endParallelFor(kpID);
    }
  }

  void destroy_shared_allocation() {}
};

template <class DeviceType, class ValueType>
struct ViewValueFunctorSequentialHostInit {
  using ExecSpace = typename DeviceType::execution_space;
  using MemSpace  = typename DeviceType::memory_space;
  static_assert(SpaceAccessibility<HostSpace, MemSpace>::accessible);

  ValueType* ptr;
  size_t n;

  ViewValueFunctorSequentialHostInit() = default;

  ViewValueFunctorSequentialHostInit(ExecSpace const& /*arg_space*/,
                                     ValueType* const arg_ptr,
                                     size_t const arg_n,
                                     std::string /*arg_name*/)
      : ptr(arg_ptr), n(arg_n) {}

  ViewValueFunctorSequentialHostInit(ValueType* const arg_ptr,
                                     size_t const arg_n,
                                     std::string /*arg_name*/)
      : ptr(arg_ptr), n(arg_n) {}

  void construct_shared_allocation() {
    if constexpr (std::is_trivial_v<ValueType>) {
      // value-initialization is equivalent to filling with zeros
      std::memset(static_cast<void*>(ptr), 0, n * sizeof(ValueType));
    } else {
      for (size_t i = 0; i < n; ++i) {
        new (ptr + i) ValueType();
      }
    }
  }

  void destroy_shared_allocation() {
    if constexpr (std::is_trivially_destructible_v<ValueType>) {
      // do nothing, don't bother calling the destructor
    } else {
      for (size_t i = 0; i < n; ++i) {
        (ptr + i)->~ValueType();
      }
    }
  }
};

}  // namespace Kokkos::Impl

#endif  // KOKKOS_VIEW_ALLOC_HPP
