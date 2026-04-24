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
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION Few(T const a[]) {
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T(a[i]);
  }
// NOLINTNEXTLINE
  template <typename T2>
  KOKKOS_INLINE_FUNCTION Few(T2 const a[]) {
    for (std::size_t i = 0; i < n; ++i) {
      new (data() + i) T(static_cast<T>(a[i]));
    }
  }
// NOLINTNEXTLINE
  template<typename... Args,
           typename = std::enable_if_t<sizeof...(Args) == n>>
  KOKKOS_INLINE_FUNCTION Few(Args... args) {
    T tmp[n] = { T(args)... };
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T(tmp[i]);
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION Few() {
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T();
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION ~Few() {
    for (std::size_t i = 0; i < n; ++i) (data()[i]).~T();
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION Few(Few<T, n> const& rhs) {
    for (std::size_t i = 0; i < n; ++i) new (data() + i) T(rhs[i]);
  }
// NOLINTNEXTLINE
  template <typename T2>
  KOKKOS_INLINE_FUNCTION Few(Few<T2, n> const& rhs) {
    for (std::size_t i = 0; i < n; ++i) {
      new (data() + i) T(static_cast<T>(rhs[i]));
    }
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION void operator=(Few<T, n> const& rhs) {
    for (std::size_t i = 0; i < n; ++i) data()[i] = rhs[i];
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION T* data() {
    return reinterpret_cast<T*>(array_);
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION T const* data() const {
    return reinterpret_cast<T const*>(array_);
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION T& operator[](std::size_t i) {
    return data()[i];
  }
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION T const& operator[](std::size_t i) const {
    return data()[i];
  }

/* ----------------------------------------------------------------------
  convenience operators
------------------------------------------------------------------------- */

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION Few<T, n>& operator+=(Few<T, n> const& rhs) {
    for (std::size_t i = 0; i < n; ++i) data()[i] += rhs[i];
    return *this;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION Few<T, n>& operator-=(Few<T, n> const& rhs) {
    for (std::size_t i = 0; i < n; ++i) data()[i] -= rhs[i];
    return *this;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION Few<T, n> operator+(Few<T, n> const& rhs) const {
    Few<T, n> res(*this); // Calls the copy constructor
    res += rhs;
    return res;
  }

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION Few<T, n> operator-(Few<T, n> const& rhs) const {
    Few<T, n> res(*this); // Calls the copy constructor
    res -= rhs;
    return res;
  }


};

#endif
