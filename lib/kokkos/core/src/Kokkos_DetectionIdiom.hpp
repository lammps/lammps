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
#ifndef KOKKOS_DETECTION_IDIOM_HPP
#define KOKKOS_DETECTION_IDIOM_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DETECTIONIDIOM
#endif

#include <Kokkos_Macros.hpp>  // FIXME doesn't actually need it if it wasn't
                              // for the header self-containment test

#include <type_traits>

// NOTE This header implements the detection idiom from Version 2 of the C++
// Extensions for Library Fundamentals, ISO/IEC TS 19568:2017

// I deliberately omitted detected_or which does not fit well with the rest
// of the specification. In my opinion, it should be removed from the TS.

namespace Kokkos {

namespace Impl {
// base class for nonesuch to inherit from so it is not an aggregate
struct nonesuch_base {};

// primary template handles all types not supporting the archetypal Op
template <class Default, class /*AlwaysVoid*/, template <class...> class Op,
          class... /*Args*/>
struct detector {
  using value_t = std::false_type;
  using type    = Default;
};

// specialization recognizes and handles only types supporting Op
template <class Default, template <class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
  using value_t = std::true_type;
  using type    = Op<Args...>;
};
}  // namespace Impl

struct nonesuch : private Impl::nonesuch_base {
  ~nonesuch()               = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

template <template <class...> class Op, class... Args>
using is_detected =
    typename Impl::detector<nonesuch, void, Op, Args...>::value_t;

template <template <class...> class Op, class... Args>
using detected_t = typename Impl::detector<nonesuch, void, Op, Args...>::type;

template <class Default, template <class...> class Op, class... Args>
using detected_or_t = typename Impl::detector<Default, void, Op, Args...>::type;

template <class Expected, template <class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

template <class To, template <class...> class Op, class... Args>
using is_detected_convertible =
    std::is_convertible<detected_t<Op, Args...>, To>;

template <template <class...> class Op, class... Args>
inline constexpr bool is_detected_v = is_detected<Op, Args...>::value;

template <class Expected, template <class...> class Op, class... Args>
inline constexpr bool is_detected_exact_v =
    is_detected_exact<Expected, Op, Args...>::value;

template <class Expected, template <class...> class Op, class... Args>
inline constexpr bool is_detected_convertible_v =
    is_detected_convertible<Expected, Op, Args...>::value;

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DETECTIONIDIOM
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DETECTIONIDIOM
#endif
#endif
