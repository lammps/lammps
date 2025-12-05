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
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_MDSPAN_ACCESSOR_HPP
#define KOKKOS_MDSPAN_ACCESSOR_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>
#include <Kokkos_Core_fwd.hpp>
#include <desul/atomics.hpp>

namespace Kokkos {

// For now use the accessors in Impl namespace, as an
// implementation detail for rebasing View on mdspan
namespace Impl {

template <class MemorySpace, class NestedAccessor>
struct SpaceAwareAccessor {
  // Part of Accessor Requirements
  using element_type     = typename NestedAccessor::element_type;
  using reference        = typename NestedAccessor::reference;
  using data_handle_type = typename NestedAccessor::data_handle_type;
  using offset_policy =
      SpaceAwareAccessor<MemorySpace, typename NestedAccessor::offset_policy>;

  // Specific to SpaceAwareAccessor
  using memory_space         = MemorySpace;
  using nested_accessor_type = NestedAccessor;

  static_assert(is_memory_space_v<memory_space>);

  KOKKOS_DEFAULTED_FUNCTION
  constexpr SpaceAwareAccessor() = default;

  template <
      class OtherMemorySpace, class OtherNestedAccessorType,
      std::enable_if_t<
          MemorySpaceAccess<MemorySpace, OtherMemorySpace>::assignable &&
              std::is_constructible_v<NestedAccessor, OtherNestedAccessorType>,
          int> = 0>
  KOKKOS_FUNCTION constexpr SpaceAwareAccessor(
      const SpaceAwareAccessor<OtherMemorySpace, OtherNestedAccessorType>&
          other) noexcept
      : nested_acc(other.nested_acc) {}

  KOKKOS_FUNCTION
  SpaceAwareAccessor(const NestedAccessor& acc) : nested_acc(acc) {}

  KOKKOS_FUNCTION
  explicit operator NestedAccessor() const { return nested_acc; }

  KOKKOS_FUNCTION
  constexpr reference access(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(
        "Kokkos::SpaceAwareAccessor ERROR: attempt to access inaccessible "
        "memory space");
    return nested_acc.access(p, i);
  }

  KOKKOS_FUNCTION
  constexpr typename offset_policy::data_handle_type offset(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return nested_acc.offset(p, i);
  }

  // Canonical way for accessing nested accessor see ISO C++
  // [linalg.scaled.scaledaccessor]
  KOKKOS_FUNCTION
  constexpr const NestedAccessor& nested_accessor() const noexcept {
    return nested_acc;
  }

 private:
// We either compile with our custom mdspan impl
// in which case we discover inside it whether no_unique_address
// works, or we use C++23 in which case it better be available
#ifdef MDSPAN_IMPL_NO_UNIQUE_ADDRESS
  MDSPAN_IMPL_NO_UNIQUE_ADDRESS
#else
  [[no_unique_address]]
#endif
  NestedAccessor nested_acc;
  template <class, class>
  friend struct SpaceAwareAccessor;
};

template <class NestedAccessor>
struct SpaceAwareAccessor<AnonymousSpace, NestedAccessor> {
  // Part of Accessor Requirements
  using element_type     = typename NestedAccessor::element_type;
  using reference        = typename NestedAccessor::reference;
  using data_handle_type = typename NestedAccessor::data_handle_type;

  using offset_policy =
      SpaceAwareAccessor<AnonymousSpace,
                         typename NestedAccessor::offset_policy>;

  // Specific to SpaceAwareAccessor
  using memory_space         = AnonymousSpace;
  using nested_accessor_type = NestedAccessor;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr SpaceAwareAccessor() = default;

  template <class OtherMemorySpace, class OtherNestedAccessorType,
            std::enable_if_t<std::is_constructible_v<NestedAccessor,
                                                     OtherNestedAccessorType>,
                             int> = 0>
  KOKKOS_FUNCTION constexpr SpaceAwareAccessor(
      const SpaceAwareAccessor<OtherMemorySpace, OtherNestedAccessorType>&
          other) noexcept
      : nested_acc(other.nested_acc) {}

  KOKKOS_FUNCTION
  SpaceAwareAccessor(const NestedAccessor& acc) : nested_acc(acc) {}

  KOKKOS_FUNCTION
  explicit operator NestedAccessor() const { return nested_acc; }

  KOKKOS_FUNCTION
  constexpr reference access(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return nested_acc.access(p, i);
  }

  KOKKOS_FUNCTION
  constexpr typename offset_policy::data_handle_type offset(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return nested_acc.offset(p, i);
  }

  // Canonical way for accessing nested accessor see ISO C++
  // [linalg.scaled.scaledaccessor]
  KOKKOS_FUNCTION
  constexpr const NestedAccessor& nested_accessor() const noexcept {
    return nested_acc;
  }

 private:
// We either compile with our custom mdspan impl
// in which case we discover inside it whether no_unique_address
// works, or we use C++23 in which case it better be available
#ifdef MDSPAN_IMPL_NO_UNIQUE_ADDRESS
  MDSPAN_IMPL_NO_UNIQUE_ADDRESS
#else
  [[no_unique_address]]
#endif
  NestedAccessor nested_acc;
  template <class, class>
  friend struct SpaceAwareAccessor;
};

// Like atomic_accessor_relaxed proposed for ISO C++26 but with
// defaulted memory scope - similar to how desul's AtomicRef has a memory scope
template <class ElementType, class MemoryScope = desul::MemoryScopeDevice>
struct AtomicAccessorRelaxed {
  using element_type = ElementType;
  using reference =
      desul::AtomicRef<ElementType, desul::MemoryOrderRelaxed, MemoryScope>;
  using data_handle_type = ElementType*;
  using offset_policy    = AtomicAccessorRelaxed;

  KOKKOS_DEFAULTED_FUNCTION
  AtomicAccessorRelaxed() = default;

  // Conversions from non-const to const element type
  template <class OtherElementType,
            std::enable_if_t<std::is_convertible_v<
                OtherElementType (*)[], element_type (*)[]>>* = nullptr>
  KOKKOS_FUNCTION constexpr AtomicAccessorRelaxed(
      Kokkos::default_accessor<OtherElementType>) noexcept {}

  template <class OtherElementType,
            std::enable_if_t<std::is_convertible_v<
                OtherElementType (*)[], element_type (*)[]>>* = nullptr>
  KOKKOS_FUNCTION constexpr AtomicAccessorRelaxed(
      AtomicAccessorRelaxed<OtherElementType, MemoryScope>) noexcept {}

  template <class OtherElementType,
            std::enable_if_t<std::is_convertible_v<
                element_type (*)[], OtherElementType (*)[]>>* = nullptr>
  KOKKOS_FUNCTION explicit operator default_accessor<OtherElementType>() const {
    return default_accessor<OtherElementType>{};
  }

  KOKKOS_FUNCTION
  reference access(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return reference(p[i]);
  }

  KOKKOS_FUNCTION
  data_handle_type offset(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return p + i;
  }
};

//=====================================================================
//============= Reference Counted Accessor and DataHandle =============
//=====================================================================

template <class ElementType, class MemorySpace>
class ReferenceCountedDataHandle {
 public:
  using value_type   = ElementType;
  using pointer      = value_type*;
  using reference    = value_type&;
  using memory_space = MemorySpace;

  KOKKOS_DEFAULTED_FUNCTION
  ReferenceCountedDataHandle() = default;

  // this only ever works on host
  explicit ReferenceCountedDataHandle(SharedAllocationRecord<void, void>* rec) {
    m_tracker.assign_allocated_record_to_uninitialized(rec);
    m_handle = static_cast<pointer>(get_record()->data());
  }

  KOKKOS_FUNCTION
  ReferenceCountedDataHandle(const SharedAllocationTracker& tracker,
                             pointer data_handle)
      : m_tracker(tracker), m_handle(data_handle) {}

  // unmanaged ctor
  template <class OtherElementType,
            class = std::enable_if_t<std::is_convertible_v<
                OtherElementType (*)[], value_type (*)[]>>>
  KOKKOS_FUNCTION ReferenceCountedDataHandle(OtherElementType* ptr)
      : m_tracker(), m_handle(ptr) {}

  // subview ctor
  template <class OtherElementType,
            class = std::enable_if_t<std::is_convertible_v<
                OtherElementType (*)[], value_type (*)[]>>>
  KOKKOS_FUNCTION ReferenceCountedDataHandle(
      const ReferenceCountedDataHandle& other, OtherElementType* ptr)
      : m_tracker(other.m_tracker), m_handle(ptr) {}

  // converting ctor
  template <class OtherElementType,
            class = std::enable_if_t<std::is_convertible_v<
                OtherElementType (*)[], value_type (*)[]>>>
  KOKKOS_FUNCTION ReferenceCountedDataHandle(
      const ReferenceCountedDataHandle<OtherElementType, memory_space>& other)
      : m_tracker(other.m_tracker), m_handle(other.m_handle) {}

  template <
      class OtherElementType, class OtherSpace,
      class = std::enable_if_t<
          std::is_convertible_v<OtherElementType (*)[], value_type (*)[]> &&
          SpaceAccessibility<memory_space,
                             typename OtherSpace::memory_space>::assignable>>
  KOKKOS_FUNCTION ReferenceCountedDataHandle(
      const ReferenceCountedDataHandle<OtherElementType, OtherSpace>& other)
      : m_tracker(other.m_tracker), m_handle(other.m_handle) {}

  KOKKOS_FUNCTION
  pointer get() const noexcept { return m_handle; }
  KOKKOS_FUNCTION
  explicit operator pointer() const noexcept { return m_handle; }

  bool has_record() const { return m_tracker.has_record(); }
  auto* get_record() const { return m_tracker.get_record<memory_space>(); }
  int use_count() const noexcept { return m_tracker.use_count(); }

  std::string get_label() const { return m_tracker.get_label<memory_space>(); }
  KOKKOS_FUNCTION
  const SharedAllocationTracker& tracker() const noexcept { return m_tracker; }

  KOKKOS_FUNCTION
  friend bool operator==(const ReferenceCountedDataHandle& lhs,
                         const value_type* rhs) {
    return lhs.m_handle == rhs;
  }

  KOKKOS_FUNCTION
  friend bool operator==(const value_type* lhs,
                         const ReferenceCountedDataHandle& rhs) {
    return lhs == rhs.m_handle;
  }

 private:
  template <class OtherElementType, class OtherSpace>
  friend class ReferenceCountedDataHandle;

  template <class OtherElementType, class OtherSpace, class NestedAccessor>
  friend class ReferenceCountedAccessor;

  SharedAllocationTracker m_tracker;
  pointer m_handle = nullptr;
};

// Helper function used by View to extract raw pointer from data_handle
template <class ElementType, class MemorySpace>
KOKKOS_INLINE_FUNCTION constexpr auto ptr_from_data_handle(
    const ReferenceCountedDataHandle<ElementType, MemorySpace>& handle) {
  return handle.get();
}

template <class T>
struct IsReferenceCountedDataHandle : std::false_type {};

template <class ElementType, class MemorySpace>
struct IsReferenceCountedDataHandle<
    ReferenceCountedDataHandle<ElementType, MemorySpace>> : std::true_type {};

template <class T>
constexpr bool IsReferenceCountedDataHandleV =
    IsReferenceCountedDataHandle<T>::value;

template <class ElementType, class MemorySpace, class NestedAccessor>
class ReferenceCountedAccessor;

template <class Accessor>
struct IsReferenceCountedAccessor : std::false_type {};

template <class ElementType, class MemorySpace, class NestedAccessor>
struct IsReferenceCountedAccessor<
    ReferenceCountedAccessor<ElementType, MemorySpace, NestedAccessor>>
    : std::true_type {};

template <class T>
constexpr bool IsReferenceCountedAccessorV =
    IsReferenceCountedAccessor<T>::value;

template <class ElementType, class MemorySpace, class NestedAccessor>
class ReferenceCountedAccessor {
 public:
  using element_type     = ElementType;
  using data_handle_type = ReferenceCountedDataHandle<ElementType, MemorySpace>;
  using reference        = typename NestedAccessor::reference;
  using offset_policy =
      ReferenceCountedAccessor<ElementType, MemorySpace,
                               typename NestedAccessor::offset_policy>;
  using memory_space = MemorySpace;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr ReferenceCountedAccessor() noexcept = default;

  template <
      class OtherElementType, class OtherNestedAccessor,
      class = std::enable_if_t<
          std::is_convertible_v<OtherElementType (*)[], element_type (*)[]> &&
          std::is_constructible_v<NestedAccessor, OtherNestedAccessor>>>
  KOKKOS_FUNCTION constexpr ReferenceCountedAccessor(
      const ReferenceCountedAccessor<OtherElementType, MemorySpace,
                                     OtherNestedAccessor>&) {}

  template <
      class OtherElementType, class OtherSpace, class OtherNestedAccessor,
      class = std::enable_if_t<
          std::is_convertible_v<OtherElementType (*)[], element_type (*)[]> &&
          SpaceAccessibility<memory_space,
                             typename OtherSpace::memory_space>::assignable &&
          std::is_constructible_v<NestedAccessor, OtherNestedAccessor>>>
  KOKKOS_FUNCTION constexpr ReferenceCountedAccessor(
      const ReferenceCountedAccessor<OtherElementType, OtherSpace,
                                     OtherNestedAccessor>&) {}

  template <class OtherElementType,
            class = std::enable_if_t<std::is_convertible_v<
                OtherElementType (*)[], element_type (*)[]>>>
  KOKKOS_FUNCTION constexpr ReferenceCountedAccessor(
      const default_accessor<OtherElementType>&) {}

  template <class DstAccessor,
            typename = std::enable_if_t<
                !IsReferenceCountedAccessor<DstAccessor>::value &&
                std::is_convertible_v<NestedAccessor, DstAccessor>>>
  KOKKOS_FUNCTION operator DstAccessor() const {
    return m_nested_acc;
  }

  KOKKOS_FUNCTION
  constexpr reference access(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const {
    return m_nested_acc.access(p.get(), i);
  }

  KOKKOS_FUNCTION
  constexpr data_handle_type offset(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OpenACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const {
    return data_handle_type{p, m_nested_acc.offset(p.get(), i)};
  }

  KOKKOS_FUNCTION
  constexpr auto nested_accessor() const { return m_nested_acc; }

 private:
#ifdef MDSPAN_IMPL_NO_UNIQUE_ADDRESS
  MDSPAN_IMPL_NO_UNIQUE_ADDRESS
#else
  [[no_unique_address]]
#endif
  NestedAccessor m_nested_acc;
};

template <class ElementType, class MemorySpace>
using CheckedReferenceCountedAccessor =
    SpaceAwareAccessor<MemorySpace,
                       ReferenceCountedAccessor<ElementType, MemorySpace,
                                                default_accessor<ElementType>>>;

template <class ElementType, class MemorySpace,
          class MemoryScope = desul::MemoryScopeDevice>
using CheckedRelaxedAtomicAccessor =
    SpaceAwareAccessor<MemorySpace, AtomicAccessorRelaxed<ElementType>>;

template <class ElementType, class MemorySpace,
          class MemoryScope = desul::MemoryScopeDevice>
using CheckedReferenceCountedRelaxedAtomicAccessor = SpaceAwareAccessor<
    MemorySpace, ReferenceCountedAccessor<ElementType, MemorySpace,
                                          AtomicAccessorRelaxed<ElementType>>>;

}  // namespace Impl
}  // namespace Kokkos

#endif
