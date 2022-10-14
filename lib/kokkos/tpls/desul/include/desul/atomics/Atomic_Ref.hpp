/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMIC_REF_IMPL_HPP_
#define DESUL_ATOMIC_REF_IMPL_HPP_

#include <cstddef>
#include <memory>
#include <type_traits>

#include "desul/atomics/Common.hpp"
#include "desul/atomics/Generic.hpp"
#include "desul/atomics/Macros.hpp"

namespace desul {
namespace Impl {

// TODO current implementation is missing the following:
// * member functions
//   * wait
//   * notify_one
//   * notify_all

template <typename T,
          typename MemoryOrder,
          typename MemoryScope,
          bool = std::is_integral<T>{},
          bool = std::is_floating_point<T>{}>
struct basic_atomic_ref;

// base class for non-integral, non-floating-point, non-pointer types
template <typename T, typename MemoryOrder, typename MemoryScope>
struct basic_atomic_ref<T, MemoryOrder, MemoryScope, false, false> {
  static_assert(std::is_trivially_copyable<T>{}, "");

 private:
  T* _ptr;

  // 1/2/4/8/16-byte types must be aligned to at least their size
  static constexpr int _min_alignment = (sizeof(T) & (sizeof(T) - 1)) || sizeof(T) > 16
                                            ? 0
                                            : sizeof(T);

 public:
  using value_type = T;

  static constexpr bool is_always_lock_free = atomic_always_lock_free(sizeof(T));

  static constexpr std::size_t required_alignment = _min_alignment > alignof(T)
                                                        ? _min_alignment
                                                        : alignof(T);

  basic_atomic_ref() = delete;
  basic_atomic_ref& operator=(basic_atomic_ref const&) = delete;

  basic_atomic_ref(basic_atomic_ref const&) = default;

  explicit basic_atomic_ref(T& obj) : _ptr(std::addressof(obj)) {}

  T operator=(T desired) const noexcept {
    this->store(desired);
    return desired;
  }

  operator T() const noexcept { return this->load(); }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION void store(T desired,
                            _MemoryOrder order = _MemoryOrder()) const noexcept {
    atomic_store(_ptr, desired, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T load(_MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T exchange(T desired,
                            _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, desired, order, MemoryScope());
  }

  DESUL_FUNCTION bool is_lock_free() const noexcept {
    return atomic_is_lock_free<sizeof(T), required_alignment>();
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(T& expected,
                                            T desired,
                                            SuccessMemoryOrder success,
                                            FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_weak(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(
      T& expected, T desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_weak(expected,
                                 desired,
                                 order,
                                 cmpexch_failure_memory_order<_MemoryOrder>(),
                                 MemoryScope());
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T& expected,
      T desired,
      SuccessMemoryOrder success,
      FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_strong(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T& expected, T desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_strong(expected,
                                   desired,
                                   order,
                                   cmpexch_failure_memory_order<_MemoryOrder>(),
                                   MemoryScope());
  }
};

// base class for atomic_ref<integral-type>
template <typename T, typename MemoryOrder, typename MemoryScope>
struct basic_atomic_ref<T, MemoryOrder, MemoryScope, true, false> {
  static_assert(std::is_integral<T>{}, "");

 private:
  T* _ptr;

 public:
  using value_type = T;
  using difference_type = value_type;

  static constexpr bool is_always_lock_free = atomic_always_lock_free(sizeof(T));

  static constexpr std::size_t required_alignment = sizeof(T) > alignof(T) ? sizeof(T)
                                                                           : alignof(T);

  basic_atomic_ref() = delete;
  basic_atomic_ref& operator=(basic_atomic_ref const&) = delete;

  explicit basic_atomic_ref(T& obj) : _ptr(&obj) {}

  basic_atomic_ref(basic_atomic_ref const&) = default;

  T operator=(T desired) const noexcept {
    this->store(desired);
    return desired;
  }

  operator T() const noexcept { return this->load(); }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION void store(T desired,
                            _MemoryOrder order = _MemoryOrder()) const noexcept {
    atomic_store(_ptr, desired, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T load(_MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T exchange(T desired,
                            _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, desired, order, MemoryScope());
  }

  DESUL_FUNCTION bool is_lock_free() const noexcept {
    return atomic_is_lock_free<sizeof(T), required_alignment>();
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(T& expected,
                                            T desired,
                                            SuccessMemoryOrder success,
                                            FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_weak(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(
      T& expected, T desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_weak(expected,
                                 desired,
                                 order,
                                 cmpexch_failure_memory_order<_MemoryOrder>(),
                                 MemoryScope());
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T& expected,
      T desired,
      SuccessMemoryOrder success,
      FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_strong(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T& expected, T desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_strong(expected,
                                   desired,
                                   order,
                                   cmpexch_failure_memory_order<_MemoryOrder>(),
                                   MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_add(value_type arg, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_add(_ptr, arg, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_sub(value_type arg, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_sub(_ptr, arg, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_and(value_type arg, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_and(_ptr, arg, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_or(value_type arg, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_or(_ptr, arg, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_xor(value_type arg, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_xor(_ptr, arg, order, MemoryScope());
  }

  DESUL_FUNCTION value_type operator++() const noexcept {
    return atomic_add_fetch(_ptr, value_type(1), MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator++(int) const noexcept { return fetch_add(1); }

  DESUL_FUNCTION value_type operator--() const noexcept {
    return atomic_sub_fetch(_ptr, value_type(1), MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator--(int) const noexcept { return fetch_sub(1); }

  DESUL_FUNCTION value_type operator+=(value_type arg) const noexcept {
    atomic_add_fetch(_ptr, arg, MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator-=(value_type arg) const noexcept {
    atomic_sub_fetch(_ptr, arg, MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator&=(value_type arg) const noexcept {
    atomic_and_fetch(_ptr, arg, MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator|=(value_type arg) const noexcept {
    atomic_or_fetch(_ptr, arg, MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator^=(value_type arg) const noexcept {
    atomic_xor_fetch(_ptr, arg, MemoryOrder(), MemoryScope());
  }
};

// base class for atomic_ref<floating-point-type>
template <typename T, typename MemoryOrder, typename MemoryScope>
struct basic_atomic_ref<T, MemoryOrder, MemoryScope, false, true> {
  static_assert(std::is_floating_point<T>{}, "");

 private:
  T* _ptr;

 public:
  using value_type = T;
  using difference_type = value_type;

  static constexpr bool is_always_lock_free = atomic_always_lock_free(sizeof(T));

  static constexpr std::size_t required_alignment = alignof(T);

  basic_atomic_ref() = delete;
  basic_atomic_ref& operator=(basic_atomic_ref const&) = delete;

  explicit basic_atomic_ref(T& obj) : _ptr(&obj) {}

  basic_atomic_ref(basic_atomic_ref const&) = default;

  T operator=(T desired) const noexcept {
    this->store(desired);
    return desired;
  }

  operator T() const noexcept { return this->load(); }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION void store(T desired,
                            _MemoryOrder order = _MemoryOrder()) const noexcept {
    atomic_store(_ptr, desired, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T load(_MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T exchange(T desired,
                            _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, desired, order, MemoryScope());
  }

  DESUL_FUNCTION bool is_lock_free() const noexcept {
    return atomic_is_lock_free<sizeof(T), required_alignment>();
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(T& expected,
                                            T desired,
                                            SuccessMemoryOrder success,
                                            FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_weak(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(
      T& expected, T desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_weak(expected,
                                 desired,
                                 order,
                                 cmpexch_failure_memory_order<_MemoryOrder>(),
                                 MemoryScope());
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T& expected,
      T desired,
      SuccessMemoryOrder success,
      FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_strong(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T& expected, T desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_strong(expected,
                                   desired,
                                   order,
                                   cmpexch_failure_memory_order<_MemoryOrder>(),
                                   MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_add(value_type arg, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_add(_ptr, arg, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_sub(value_type arg, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_sub(_ptr, arg, order, MemoryScope());
  }

  DESUL_FUNCTION value_type operator+=(value_type arg) const noexcept {
    atomic_add_fetch(_ptr, arg, MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator-=(value_type arg) const noexcept {
    atomic_sub_fetch(_ptr, arg, MemoryOrder(), MemoryScope());
  }
};

// base class for atomic_ref<pointer-type>
template <typename T, typename MemoryOrder, typename MemoryScope>
struct basic_atomic_ref<T*, MemoryOrder, MemoryScope, false, false> {
 private:
  T** _ptr;

 public:
  using value_type = T*;
  using difference_type = std::ptrdiff_t;

  static constexpr bool is_always_lock_free = atomic_always_lock_free(sizeof(T));

  static constexpr std::size_t required_alignment = alignof(T*);

  basic_atomic_ref() = delete;
  basic_atomic_ref& operator=(basic_atomic_ref const&) = delete;

  explicit basic_atomic_ref(T*& arg) : _ptr(std::addressof(arg)) {}

  basic_atomic_ref(basic_atomic_ref const&) = default;

  T* operator=(T* desired) const noexcept {
    this->store(desired);
    return desired;
  }

  operator T*() const noexcept { return this->load(); }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION void store(T* desired,
                            _MemoryOrder order = _MemoryOrder()) const noexcept {
    atomic_store(_ptr, desired, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T* load(_MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION T* exchange(T* desired,
                             _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_load(_ptr, desired, order, MemoryScope());
  }

  DESUL_FUNCTION bool is_lock_free() const noexcept {
    return atomic_is_lock_free<sizeof(T*), required_alignment>();
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(T*& expected,
                                            T* desired,
                                            SuccessMemoryOrder success,
                                            FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_weak(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_weak(
      T*& expected, T* desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_weak(expected,
                                 desired,
                                 order,
                                 cmpexch_failure_memory_order<_MemoryOrder>(),
                                 MemoryScope());
  }

  template <typename SuccessMemoryOrder, typename FailureMemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T*& expected,
      T* desired,
      SuccessMemoryOrder success,
      FailureMemoryOrder failure) const noexcept {
    return atomic_compare_exchange_strong(
        _ptr, expected, desired, success, failure, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION bool compare_exchange_strong(
      T*& expected, T* desired, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return compare_exchange_strong(expected,
                                   desired,
                                   order,
                                   cmpexch_failure_memory_order<_MemoryOrder>(),
                                   MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_add(difference_type d, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_add(_ptr, _type_size(d), order, MemoryScope());
  }

  template <typename _MemoryOrder = MemoryOrder>
  DESUL_FUNCTION value_type
  fetch_sub(difference_type d, _MemoryOrder order = _MemoryOrder()) const noexcept {
    return atomic_fetch_sub(_ptr, _type_size(d), order, MemoryScope());
  }

  DESUL_FUNCTION value_type operator++() const noexcept {
    return atomic_add_fetch(_ptr, _type_size(1), MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator++(int) const noexcept { return fetch_add(1); }

  DESUL_FUNCTION value_type operator--() const noexcept {
    return atomic_sub_fetch(_ptr, _type_size(1), MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator--(int) const noexcept { return fetch_sub(1); }

  DESUL_FUNCTION value_type operator+=(difference_type d) const noexcept {
    atomic_add_fetch(_ptr, _type_size(d), MemoryOrder(), MemoryScope());
  }

  DESUL_FUNCTION value_type operator-=(difference_type d) const noexcept {
    atomic_sub_fetch(_ptr, _type_size(d), MemoryOrder(), MemoryScope());
  }

 private:
  static constexpr std::ptrdiff_t _type_size(std::ptrdiff_t d) noexcept {
    static_assert(std::is_object<T>{}, "");
    return d * sizeof(T);
  }
};

}  // namespace Impl

template <typename T, typename MemoryOrder, typename MemoryScope>
struct scoped_atomic_ref : Impl::basic_atomic_ref<T, MemoryOrder, MemoryScope> {
  explicit scoped_atomic_ref(T& obj) noexcept
      : Impl::basic_atomic_ref<T, MemoryOrder, MemoryScope>(obj) {}

  scoped_atomic_ref& operator=(scoped_atomic_ref const&) = delete;

  scoped_atomic_ref(scoped_atomic_ref const&) = default;

  using Impl::basic_atomic_ref<T, MemoryOrder, MemoryScope>::operator=;
};

}  // namespace desul

#endif
