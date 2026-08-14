// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_VIEWCOMMONTYPE_HPP
#define KOKKOS_VIEWCOMMONTYPE_HPP

#include <type_traits>
#include <Kokkos_Macros.hpp>

namespace Kokkos {
namespace Impl {

template <class Specialize, typename A, typename B>
struct CommonViewValueType;

template <typename A, typename B>
struct CommonViewValueType<void, A, B> {
  using value_type = std::common_type_t<A, B>;
};

template <class Specialize, class ValueType>
struct CommonViewAllocProp;

template <class ValueType>
struct CommonViewAllocProp<void, ValueType> {
  using value_type = ValueType;
  using data_type  = ValueType;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_5
  using scalar_array_type KOKKOS_DEPRECATED_WITH_COMMENT(
      "Use data_type instead.") = data_type;
#endif

  template <class... Views>
  KOKKOS_INLINE_FUNCTION CommonViewAllocProp(const Views&...) {}
};

template <class... Views>
struct DeduceCommonViewAllocProp;

// Base case must provide types for:
// 1. specialize  2. value_type  3. is_view  4. prop_type
template <class FirstView>
struct DeduceCommonViewAllocProp<FirstView> {
  using specialize = typename FirstView::traits::specialize;

  using value_type = typename FirstView::traits::value_type;

  enum : bool { is_view = is_view<FirstView>::value };

  using prop_type = CommonViewAllocProp<specialize, value_type>;
};

template <class FirstView, class... NextViews>
struct DeduceCommonViewAllocProp<FirstView, NextViews...> {
  using NextTraits = DeduceCommonViewAllocProp<NextViews...>;

  using first_specialize = typename FirstView::traits::specialize;
  using first_value_type = typename FirstView::traits::value_type;

  enum : bool { first_is_view = is_view<FirstView>::value };

  using next_specialize = typename NextTraits::specialize;
  using next_value_type = typename NextTraits::value_type;

  enum : bool { next_is_view = NextTraits::is_view };

  // common types

  // determine specialize type
  // if first and next specialize differ, but are not the same specialize, error
  // out
  static_assert(!(!std::is_same_v<first_specialize, next_specialize> &&
                  !std::is_void_v<first_specialize> &&
                  !std::is_void_v<next_specialize>),
                "Kokkos DeduceCommonViewAllocProp ERROR: Only one non-void "
                "specialize trait allowed");

  // otherwise choose non-void specialize if either/both are non-void
  using specialize =
      std::conditional_t<std::is_same_v<first_specialize, next_specialize>,
                         first_specialize,
                         std::conditional_t<(std::is_void_v<first_specialize> &&
                                             !std::is_void_v<next_specialize>),
                                            next_specialize, first_specialize>>;

  using value_type = typename CommonViewValueType<specialize, first_value_type,
                                                  next_value_type>::value_type;

  enum : bool { is_view = (first_is_view && next_is_view) };

  using prop_type = CommonViewAllocProp<specialize, value_type>;
};

}  // end namespace Impl

template <class... Views>
using DeducedCommonPropsType =
    typename Impl::DeduceCommonViewAllocProp<Views...>::prop_type;

// This function is required in certain scenarios where users customize
// Kokkos View internals. One example are dynamic length embedded ensemble
// types. The function is used to propagate necessary information
// (like the ensemble size) when creating new views.
// However, most of the time it is called with a single view.
// Furthermore, the propagated information is not just for view allocations.
// From what I can tell, the type of functionality provided by
// common_view_alloc_prop is the equivalent of propagating accessors in mdspan,
// a mechanism we will eventually use to replace this clunky approach here, when
// we are finally mdspan based.
// TODO: get rid of this when we have mdspan
template <class... Views>
KOKKOS_INLINE_FUNCTION DeducedCommonPropsType<Views...> common_view_alloc_prop(
    Views const&... views) {
  return DeducedCommonPropsType<Views...>(views...);
}

}  // namespace Kokkos

#endif  // KOKKOS_VIEWCOMMONTYPE_HPP
