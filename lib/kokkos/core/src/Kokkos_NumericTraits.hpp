// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NUMERIC_TRAITS_HPP
#define KOKKOS_NUMERIC_TRAITS_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERIC_TRAITS
#endif

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#include <Kokkos_ReductionIdentity.hpp>
#endif
#include <type_traits>
#include <limits>

namespace Kokkos::Experimental {

#define KOKKOS_IMPL_DEFINE_TRAIT(TRAIT, NUMERIC_LIMITS_MEMBER, CONSTRAINT)  \
  namespace Impl {                                                          \
  template <class T, class Enable = void>                                   \
  struct TRAIT##_helper {};                                                 \
  template <class T>                                                        \
  struct TRAIT##_helper<T, std::enable_if_t<std::is_##CONSTRAINT##_v<T>>> { \
    static constexpr auto value =                                           \
        std::numeric_limits<T>::NUMERIC_LIMITS_MEMBER;                      \
  };                                                                        \
  }                                                                         \
  template <class T>                                                        \
  struct TRAIT : Impl::TRAIT##_helper<T> {};                                \
  template <class T>                                                        \
  inline constexpr auto TRAIT##_v = TRAIT<T>::value;

// clang-format off
// Numeric distinguished value traits
KOKKOS_IMPL_DEFINE_TRAIT(infinity,       infinity(),      floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(finite_min,     lowest(),        arithmetic    )
KOKKOS_IMPL_DEFINE_TRAIT(finite_max,     max(),           arithmetic    )
KOKKOS_IMPL_DEFINE_TRAIT(epsilon,        epsilon(),       floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(round_error,    round_error(),   floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(norm_min,       min(),           floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(denorm_min,     denorm_min(),    floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(quiet_NaN,      quiet_NaN(),     floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(signaling_NaN,  signaling_NaN(), floating_point)

// Numeric characteristics traits
KOKKOS_IMPL_DEFINE_TRAIT(digits,         digits,          arithmetic    )
KOKKOS_IMPL_DEFINE_TRAIT(digits10,       digits10,        arithmetic    )
KOKKOS_IMPL_DEFINE_TRAIT(max_digits10,   max_digits10,    floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(radix,          radix,           arithmetic    )
KOKKOS_IMPL_DEFINE_TRAIT(min_exponent,   min_exponent,    floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(min_exponent10, min_exponent10,  floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(max_exponent,   max_exponent,    floating_point)
KOKKOS_IMPL_DEFINE_TRAIT(max_exponent10, max_exponent10,  floating_point)
// clang-format on

#undef KOKKOS_IMPL_DEFINE_TRAIT

}  // namespace Kokkos::Experimental

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERIC_TRAITS
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_NUMERIC_TRAITS
#endif
#endif
