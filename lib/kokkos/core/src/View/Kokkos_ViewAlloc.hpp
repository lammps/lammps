// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
#include <optional>

#include <impl/Kokkos_Tools.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

namespace Kokkos::Impl {

template <class DeviceType, class ValueType>
struct ViewValueFunctor {
  using ExecSpace = typename DeviceType::execution_space;

  struct DestroyTag {};
  struct ConstructTag {};

  ExecSpace space;
  ValueType* ptr;
  size_t n;
  std::string name;
  bool default_exec_space;

  template <class SameValueType = ValueType>
  KOKKOS_FUNCTION
      std::enable_if_t<std::is_default_constructible_v<SameValueType>>
      operator()(ConstructTag, const size_t i) const {
    new (ptr + i) ValueType();
  }

  KOKKOS_FUNCTION void operator()(DestroyTag, const size_t i) const {
    // When instantiating a View on host execution space with a host only
    // destructor the workaround for CUDA device symbol instantiation tries to
    // still compile a destruction kernel for the device, and issues a warning
    // for host from host-device
#ifdef KOKKOS_ENABLE_CUDA
    if constexpr (std::is_same_v<ExecSpace, Cuda>) {
      KOKKOS_IF_ON_DEVICE(((ptr + i)->~ValueType();))
    } else {
      KOKKOS_IF_ON_HOST(((ptr + i)->~ValueType();))
    }
#else
    (ptr + i)->~ValueType();
#endif
  }

  ViewValueFunctor()                                   = default;
  ViewValueFunctor(const ViewValueFunctor&)            = default;
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
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endParallelFor(kpID);
    }
    if (default_exec_space || std::is_same_v<Tag, DestroyTag>) {
      space.fence(std::is_same_v<Tag, DestroyTag>
                      ? "Kokkos::View::destruction before deallocate"
                      : "Kokkos::View::initialization");
    }
  }

  // Shortcut for zero initialization
  void zero_memset_implementation() {
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

    (void)ZeroMemset(space, ptr, n * sizeof(ValueType));

    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endParallelFor(kpID);
    }
    if (default_exec_space) {
      space.fence("Kokkos::View::initialization via memset");
    }
  }

  void construct_shared_allocation() {
    if constexpr (std::is_trivially_default_constructible_v<ValueType>) {
      // value-initialization is equivalent to filling with zeros
      zero_memset_implementation();
    } else {
      parallel_for_implementation<ConstructTag>();
    }
  }

  void destroy_shared_allocation() {
    if constexpr (std::is_trivially_destructible_v<ValueType>) {
      // do nothing, don't bother calling the destructor
    } else {
#ifdef KOKKOS_ENABLE_IMPL_VIEW_OF_VIEWS_DESTRUCTOR_PRECONDITION_VIOLATION_WORKAROUND
      if constexpr (std::is_same_v<typename ExecSpace::memory_space,
                                   Kokkos::HostSpace>)
        for (size_t i = 0; i < n; ++i) (ptr + i)->~ValueType();
      else
#endif
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
    if constexpr (std::is_trivially_default_constructible_v<ValueType>) {
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

template <class ElementType, class MemorySpace, class ExecutionSpace,
          bool Initialize, bool SequentialInit>
Kokkos::Impl::SharedAllocationRecord<void, void>* make_shared_allocation_record(
    const size_t& required_span_size, std::string_view label,
    const MemorySpace& memory_space,
    const std::optional<ExecutionSpace> exec_space,
    std::bool_constant<Initialize>, std::bool_constant<SequentialInit>) {
  static_assert(SpaceAccessibility<ExecutionSpace, MemorySpace>::accessible);

  // Use this for constructing and destroying the view
  using device_type  = Kokkos::Device<ExecutionSpace, MemorySpace>;
  using functor_type = std::conditional_t<
      SequentialInit,
      ViewValueFunctorSequentialHostInit<device_type, ElementType>,
      ViewValueFunctor<device_type, ElementType>>;
  using record_type =
      Kokkos::Impl::SharedAllocationRecord<MemorySpace, functor_type>;

  /* Force alignment of allocations on on 8 byte boundaries even for
   * element types smaller than 8 bytes */
  static constexpr std::size_t align_mask = 0x7;

  // Calculate the total size of the memory, in bytes, and make sure it is
  // byte-aligned
  const std::size_t alloc_size =
      (required_span_size * sizeof(ElementType) + align_mask) & ~align_mask;

  auto* record =
      exec_space
          ? record_type::allocate(*exec_space, memory_space, std::string{label},
                                  alloc_size)
          : record_type::allocate(memory_space, std::string{label}, alloc_size);

  auto ptr = static_cast<ElementType*>(record->data());

  auto functor =
      exec_space ? functor_type(*exec_space, ptr, required_span_size,
                                std::string{label})
                 : functor_type(ptr, required_span_size, std::string{label});

  //  Only initialize if the allocation is non-zero.
  //  May be zero if one of the dimensions is zero.
  if constexpr (Initialize) {
    if (alloc_size) {
      // Assume destruction is only required when construction is requested.
      // The ViewValueFunctor has both value construction and destruction
      // operators.
      record->m_destroy = std::move(functor);

      // Construct values
      record->m_destroy.construct_shared_allocation();
    }
  }

  return record;
}

}  // namespace Kokkos::Impl

#endif  // KOKKOS_VIEW_ALLOC_HPP
