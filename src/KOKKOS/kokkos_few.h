// clang-format off
#ifndef KOKKOS_FEW_H
#define KOKKOS_FEW_H

#include <Kokkos_Core.hpp>

template <typename T, std::size_t n>
class Few {
  alignas(T) char array_[n * sizeof(T)];

 public:
  enum { size = n };
  Few(std::initializer_list<T> l) {
    std::size_t i = 0;
    for (auto it = l.begin(); it != l.end(); ++it) {
      new (data() + (i++)) T(*it);
    }
  }
  KOKKOS_INLINE_FUNCTION Few(T const a[]) {
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T(a[i]);
  }
  KOKKOS_INLINE_FUNCTION Few() {
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T();
  }
  KOKKOS_INLINE_FUNCTION ~Few() {
    for (std::size_t i = 0; i < n; ++i) (data()[i]).~T();
  }
  KOKKOS_INLINE_FUNCTION Few(Few<T, n> const& rhs) {
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
  KOKKOS_INLINE_FUNCTION Few(Few<T, n> const volatile& rhs) {
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
  KOKKOS_INLINE_FUNCTION void operator=(Few<T, n> const& rhs) volatile {
    for (std::size_t i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  KOKKOS_INLINE_FUNCTION void operator=(Few<T, n> const& rhs) {
    for (std::size_t i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  KOKKOS_INLINE_FUNCTION void operator=(Few<T, n> const volatile& rhs) {
    for (std::size_t i = 0; i < n; ++i) data()[i] = rhs[i];
  }
  KOKKOS_INLINE_FUNCTION T* data() {
    return reinterpret_cast<T*>(array_);
  }
  KOKKOS_INLINE_FUNCTION T const* data() const {
    return reinterpret_cast<T const*>(array_);
  }
  KOKKOS_INLINE_FUNCTION T volatile* data() volatile {
    return reinterpret_cast<T volatile*>(array_);
  }
  KOKKOS_INLINE_FUNCTION T const volatile* data() const volatile {
    return reinterpret_cast<T const volatile*>(array_);
  }
  KOKKOS_INLINE_FUNCTION T& operator[](std::size_t i) {
    return data()[i];
  }
  KOKKOS_INLINE_FUNCTION T const& operator[](std::size_t i) const {
    return data()[i];
  }
  KOKKOS_INLINE_FUNCTION T volatile& operator[](std::size_t i) volatile {
    return data()[i];
  }
  KOKKOS_INLINE_FUNCTION T const volatile& operator[](std::size_t i) const volatile {
    return data()[i];
  }
};

#endif
