// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_BIT_MANIPULATION_HPP
#define KOKKOS_BIT_MANIPULATION_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_NumericTraits.hpp>
#include <climits>  // CHAR_BIT
#include <concepts>
#include <cstring>  //memcpy
#include <type_traits>

namespace Kokkos::Impl {

template <template <bool /*constant_evaluated*/, bool /*device*/> class Op,
          class T>
KOKKOS_FUNCTION constexpr auto dispatch_helper(T x) noexcept {
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  // __builtin_is_device_code() is non-constexpr
  return Op<true, true>::do_compute(x);
#else
#if defined(__cpp_if_consteval)  // since C++23
  if consteval {
    KOKKOS_IF_ON_HOST((return Op<true, false>::do_compute(x);))
    KOKKOS_IF_ON_DEVICE((return Op<true, true>::do_compute(x);))
  } else {
    KOKKOS_IF_ON_HOST((return Op<false, false>::do_compute(x);))
    KOKKOS_IF_ON_DEVICE((return Op<false, true>::do_compute(x);))
  }
#else
  KOKKOS_IF_ON_HOST((return Op<true, false>::do_compute(x);))
  KOKKOS_IF_ON_DEVICE((return Op<true, true>::do_compute(x);))
#endif
#endif
}

template <template <bool /*constant_evaluated*/, bool /*device*/> class Op,
          class T>
KOKKOS_FUNCTION constexpr auto dispatch_helper_builtin(T x) noexcept {
  KOKKOS_IF_ON_HOST((return Op<false, false>::do_compute(x);))
  KOKKOS_IF_ON_DEVICE((return Op<false, true>::do_compute(x);))

  // FIXME_NVHPC: erroneous warning about return from non-void function
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC)
  return T();
#endif
}

#if defined(KOKKOS_COMPILER_CLANG) || defined(KOKKOS_COMPILER_INTEL_LLVM) || \
    defined(KOKKOS_COMPILER_GNU) || defined(KOKKOS_COMPILER_APPLECC)
#define KOKKOS_IMPL_USE_GCC_BUILT_IN_FUNCTIONS
#endif

//<editor-fold desc="byteswap implementation details">
template <bool constant_evaluated, bool device>
struct ByteSwap {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    if constexpr (sizeof(T) > 1) {
      using U = std::make_unsigned_t<T>;

      size_t shift = CHAR_BIT * (sizeof(T) - 1);

      U lo_mask = static_cast<unsigned char>(~0);
      U hi_mask = lo_mask << shift;

      U val = x;

      for (size_t i = 0; i < sizeof(T) / 2; ++i) {
        U lo_val = val & lo_mask;
        U hi_val = val & hi_mask;

        val = (val & ~lo_mask) | (hi_val >> shift);
        val = (val & ~hi_mask) | (lo_val << shift);

        lo_mask <<= CHAR_BIT;
        hi_mask >>= CHAR_BIT;

        shift -= static_cast<size_t>(2) * CHAR_BIT;
      }
      return val;
    }
    // sizeof(T) == 1
    return x;
  }
};

#ifdef KOKKOS_IMPL_USE_GCC_BUILT_IN_FUNCTIONS
template <bool constant_evaluated>
struct ByteSwap<constant_evaluated, /*device=*/false> {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    if constexpr (sizeof(T) == 1) {
      return x;
    } else if constexpr (sizeof(T) == 2) {
      return __builtin_bswap16(x);
    } else if constexpr (sizeof(T) == 4) {
      return __builtin_bswap32(x);
    } else if constexpr (sizeof(T) == 8) {
      return __builtin_bswap64(x);
    } else if constexpr (sizeof(T) == 16) {
#if defined(__has_builtin)
#if __has_builtin(__builtin_bswap128)
      return __builtin_bswap128(x);
#endif
#endif
      return (__builtin_bswap64(x >> 64) |
              (static_cast<T>(__builtin_bswap64(x)) << 64));
    }
    return ByteSwap<true, false>::do_compute(x);
  }
};
#endif
//</editor-fold>

//<editor-fold desc="countl_zero implementation details">
template <bool constant_evaluated, bool device>
struct CountlZero {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    // From Hacker's Delight (2nd edition) section 5-3
    unsigned int y = 0;
    using ::Kokkos::Experimental::digits_v;
    int n = digits_v<T>;
    int c = digits_v<T> / 2;
    do {
      y = x >> c;
      if (y != 0) {
        n -= c;
        x = y;
      }
      c >>= 1;
    } while (c != 0);
    return n - static_cast<int>(x);
  }
};

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
template <>
struct CountlZero</*constant_evaluated=*/false, /*device=*/true> {
  template <class T>
  static KOKKOS_IMPL_DEVICE_FUNCTION T do_compute(T x) noexcept {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    if constexpr (sizeof(T) == sizeof(long long int))
      return __clzll(reinterpret_cast<long long int&>(x));
    if constexpr (sizeof(T) == sizeof(int))
      return __clz(reinterpret_cast<int&>(x));
    using ::Kokkos::Experimental::digits_v;
    constexpr int shift = digits_v<unsigned int> - digits_v<T>;
    return __clz(x) - shift;
#elif defined(KOKKOS_ENABLE_SYCL)
    return sycl::clz(x);
#else
    static_assert(false, "implementation bug");
#endif
  }
};
#endif

#ifdef KOKKOS_IMPL_USE_GCC_BUILT_IN_FUNCTIONS
template <bool constant_evaluated>
struct CountlZero<constant_evaluated, /*device=*/false> {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    using ::Kokkos::Experimental::digits_v;
    if (x == 0) return digits_v<T>;
    if constexpr (std::is_same_v<T, unsigned long long>) {
      return __builtin_clzll(x);
    } else if constexpr (std::is_same_v<T, unsigned long>) {
      return __builtin_clzl(x);
    } else if constexpr (std::is_same_v<T, unsigned int>) {
      return __builtin_clz(x);
    } else {
      constexpr int shift = digits_v<unsigned int> - digits_v<T>;
      return __builtin_clz(x) - shift;
    }
  }
};
#endif
//</editor-fold>

//<editor-fold desc="countr_zero implementation details">
template <bool constant_evaluated, bool device>
struct CountrZero {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    using ::Kokkos::Experimental::digits_v;
    return digits_v<T> -
           CountlZero<constant_evaluated, device>::do_compute(
               static_cast<T>(static_cast<T>(~x) & static_cast<T>(x - 1)));
  }
};

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
template <>
struct CountrZero</*constant_evaluated=*/false, /*device=*/true> {
  template <class T>
  static KOKKOS_IMPL_DEVICE_FUNCTION T do_compute(T x) noexcept {
    using ::Kokkos::Experimental::digits_v;
    if (x == 0) return digits_v<T>;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    if constexpr (sizeof(T) == sizeof(long long int))
      return __ffsll(reinterpret_cast<long long int&>(x)) - 1;
    return __ffs(reinterpret_cast<int&>(x)) - 1;
#elif defined(KOKKOS_ENABLE_SYCL)
    return sycl::ctz(x);
#else
    static_assert(false, "implementation bug");
#endif
  }
};
#endif

#ifdef KOKKOS_IMPL_USE_GCC_BUILT_IN_FUNCTIONS
template <bool constant_evaluated>
struct CountrZero<constant_evaluated, /*device=*/false> {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    using ::Kokkos::Experimental::digits_v;
    if (x == 0) return digits_v<T>;
    if constexpr (std::is_same_v<T, unsigned long long>) {
      return __builtin_ctzll(x);
    } else if constexpr (std::is_same_v<T, unsigned long>) {
      return __builtin_ctzl(x);
    } else {
      return __builtin_ctz(x);
    }
  }
};
#endif
//</editor-fold>

//<editor-fold desc="popcount implementation details">
template <bool constant_evaluated, bool device>
struct PopCount {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    int c = 0;
    for (; x != 0; x &= x - 1) {
      ++c;
    }
    return c;
  }
};

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
template <>
struct PopCount</*constant_evaluated=*/false, /*device=*/true> {
  template <class T>
  static KOKKOS_IMPL_DEVICE_FUNCTION T do_compute(T x) noexcept {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
    if constexpr (sizeof(T) == sizeof(long long int)) return __popcll(x);
    return __popc(x);
#elif defined(KOKKOS_ENABLE_SYCL)
    return sycl::popcount(x);
#else
    static_assert(false, "implementation bug");
#endif
  }
};
#endif

#ifdef KOKKOS_IMPL_USE_GCC_BUILT_IN_FUNCTIONS
template <bool constant_evaluated>
struct PopCount<constant_evaluated, /*device=*/false> {
  template <class T>
  static KOKKOS_FUNCTION constexpr T do_compute(T x) noexcept {
    if constexpr (std::is_same_v<T, unsigned long long>) {
      return __builtin_popcountll(x);
    } else if constexpr (std::is_same_v<T, unsigned long>) {
      return __builtin_popcountl(x);
    } else {
      return __builtin_popcount(x);
    }
  }
};
#endif
//</editor-fold>

#undef KOKKOS_IMPL_USE_GCC_BUILT_IN_FUNCTIONS

template <class T>
concept StandardUnsignedInteger =
    std::same_as<T, unsigned char> || std::same_as<T, unsigned short> ||
    std::same_as<T, unsigned int> || std::same_as<T, unsigned long> ||
    std::same_as<T, unsigned long long>;

}  // namespace Kokkos::Impl

namespace Kokkos {

//<editor-fold desc="[bit.cast], bit_cast">
template <class To, class From>
  requires(sizeof(To) == sizeof(From) && std::is_trivially_copyable_v<To> &&
           std::is_trivially_copyable_v<From>)
KOKKOS_FUNCTION To bit_cast(From const& from) noexcept {
#if defined(KOKKOS_ENABLE_SYCL)
  return sycl::bit_cast<To>(from);
#else
  To to;
  memcpy(static_cast<void*>(&to), static_cast<const void*>(&from), sizeof(To));
  return to;
#endif
}
//</editor-fold>

//<editor-fold desc="[bit.byteswap], byteswap">
template <std::integral T>
KOKKOS_FUNCTION constexpr T byteswap(T value) noexcept {
  return Impl::dispatch_helper<Impl::ByteSwap>(value);
}
//</editor-fold>

//<editor-fold desc="[bit.count], counting">
template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr int countl_zero(T x) noexcept {
  using ::Kokkos::Experimental::digits_v;
  if (x == 0) return digits_v<T>;
  return Impl::dispatch_helper<Impl::CountlZero>(x);
}

template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr int countl_one(T x) noexcept {
  using ::Kokkos::Experimental::digits_v;
  using ::Kokkos::Experimental::finite_max_v;
  if (x == finite_max_v<T>) return digits_v<T>;
  return countl_zero(static_cast<T>(~x));
}

template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr int countr_zero(T x) noexcept {
  using ::Kokkos::Experimental::digits_v;
  if (x == 0) return digits_v<T>;
  return Impl::dispatch_helper<Impl::CountrZero>(x);
}

template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr int countr_one(T x) noexcept {
  using ::Kokkos::Experimental::digits_v;
  using ::Kokkos::Experimental::finite_max_v;
  if (x == finite_max_v<T>) return digits_v<T>;
  return countr_zero(static_cast<T>(~x));
}

template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr int popcount(T x) noexcept {
  if (x == 0) return 0;
  return Impl::dispatch_helper<Impl::PopCount>(x);
}
//</editor-fold>

//<editor-fold desc="[bit.pow.two], integral powers of 2">
template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr bool has_single_bit(T x) noexcept {
  return x != 0 && (((x & (x - 1)) == 0));
}

template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr T bit_ceil(T x) noexcept {
  if (x <= 1) return 1;
  using ::Kokkos::Experimental::digits_v;
  return T{1} << (digits_v<T> - countl_zero(static_cast<T>(x - 1)));
}

template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr T bit_floor(T x) noexcept {
  if (x == 0) return 0;
  using ::Kokkos::Experimental::digits_v;
  return T{1} << (digits_v<T> - 1 - countl_zero(x));
}

template <Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION constexpr int bit_width(T x) noexcept {
  if (x == 0) return 0;
  using ::Kokkos::Experimental::digits_v;
  return digits_v<T> - countl_zero(x);
}
//</editor-fold>

//<editor-fold desc="[bit.rotate], rotating">
template <Impl::StandardUnsignedInteger T>
[[nodiscard]] KOKKOS_FUNCTION constexpr T rotl(T x, int s) noexcept {
  using Experimental::digits_v;
  constexpr auto dig = digits_v<T>;
  int const rem      = s % dig;
  if (rem == 0) return x;
  if (rem > 0) return (x << rem) | (x >> ((dig - rem) % dig));
  return (x >> -rem) | (x << ((dig + rem) % dig));  // rotr(x, -rem)
}

template <Impl::StandardUnsignedInteger T>
[[nodiscard]] KOKKOS_FUNCTION constexpr T rotr(T x, int s) noexcept {
  using Experimental::digits_v;
  constexpr auto dig = digits_v<T>;
  int const rem      = s % dig;
  if (rem == 0) return x;
  if (rem > 0) return (x >> rem) | (x << ((dig - rem) % dig));
  return (x << -rem) | (x >> ((dig + rem) % dig));  // rotl(x, -rem)
}
//</editor-fold>

}  // namespace Kokkos

namespace Kokkos::Experimental {

template <class To, class From>
  requires(sizeof(To) == sizeof(From) && std::is_trivially_copyable_v<To> &&
           std::is_trivially_copyable_v<From>)
KOKKOS_FUNCTION To bit_cast_builtin(From const& from) noexcept {
  // qualify the call to avoid ADL
  return Kokkos::bit_cast<To>(from);  // no benefit to call the _builtin variant
}

template <std::integral T>
KOKKOS_FUNCTION T byteswap_builtin(T x) noexcept {
  return Kokkos::Impl::dispatch_helper_builtin<Kokkos::Impl::ByteSwap>(x);
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION int countl_zero_builtin(T x) noexcept {
  if (x == 0) return digits_v<T>;
  return Kokkos::Impl::dispatch_helper_builtin<Kokkos::Impl::CountlZero>(x);
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION int countl_one_builtin(T x) noexcept {
  if (x == finite_max_v<T>) return digits_v<T>;
  return countl_zero_builtin(static_cast<T>(~x));
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION int countr_zero_builtin(T x) noexcept {
  if (x == 0) return digits_v<T>;
  return Kokkos::Impl::dispatch_helper_builtin<Kokkos::Impl::CountrZero>(x);
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION int countr_one_builtin(T x) noexcept {
  if (x == finite_max_v<T>) return digits_v<T>;
  return countr_zero_builtin(static_cast<T>(~x));
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION int popcount_builtin(T x) noexcept {
  return Kokkos::Impl::dispatch_helper_builtin<Kokkos::Impl::PopCount>(x);
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION bool has_single_bit_builtin(T x) noexcept {
  return has_single_bit(x);  // no benefit to call the _builtin variant
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION T bit_ceil_builtin(T x) noexcept {
  if (x <= 1) return 1;
  return T{1} << (digits_v<T> - countl_zero_builtin(static_cast<T>(x - 1)));
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION T bit_floor_builtin(T x) noexcept {
  if (x == 0) return 0;
  return T{1} << (digits_v<T> - 1 - countl_zero_builtin(x));
}

template <Kokkos::Impl::StandardUnsignedInteger T>
KOKKOS_FUNCTION int bit_width_builtin(T x) noexcept {
  if (x == 0) return 0;
  return digits_v<T> - countl_zero_builtin(x);
}

template <Kokkos::Impl::StandardUnsignedInteger T>
[[nodiscard]] KOKKOS_FUNCTION T rotl_builtin(T x, int s) noexcept {
  return rotl(x, s);  // no benefit to call the _builtin variant
}

template <Kokkos::Impl::StandardUnsignedInteger T>
[[nodiscard]] KOKKOS_FUNCTION T rotr_builtin(T x, int s) noexcept {
  return rotr(x, s);  // no benefit to call the _builtin variant
}

}  // namespace Kokkos::Experimental

#endif
