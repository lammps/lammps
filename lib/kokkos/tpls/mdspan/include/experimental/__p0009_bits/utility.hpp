#pragma once

#include <cstddef>
#include <type_traits>

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {
namespace detail {

// type alias used for rank-based tag dispatch
//
// this is used to enable alternatives to constexpr if when building for C++14
//
template <std::size_t N>
using with_rank = std::integral_constant<std::size_t, N>;

template <class I1, class I2>
MDSPAN_INLINE_FUNCTION
constexpr bool common_integral_compare(I1 x, I2 y)
{
  static_assert(std::is_integral<I1>::value &&
                std::is_integral<I2>::value, "");

  using I = std::common_type_t<I1, I2>;
  return static_cast<I>(x) == static_cast<I>(y);
}

template <class T1, class T2, class F>
MDSPAN_INLINE_FUNCTION
constexpr bool rankwise_equal(with_rank<0>, const T1&, const T2&, F)
{
  return true;
}

template <std::size_t N, class T1, class T2, class F>
MDSPAN_INLINE_FUNCTION
constexpr bool rankwise_equal(with_rank<N>, const T1& x, const T2& y, F func)
{
  bool match = true;

  for (std::size_t r = 0; r < N; r++) {
    match = match && common_integral_compare(func(x, r), func(y, r));
  }

  return match;
}

constexpr struct
{
  template <class T, class I>
  MDSPAN_INLINE_FUNCTION
  constexpr auto operator()(const T& x, I i) const
  {
    return x.extent(i);
  }
} extent;

constexpr struct
{
  template <class T, class I>
  MDSPAN_INLINE_FUNCTION
  constexpr auto operator()(const T& x, I i) const
  {
    return x.stride(i);
  }
} stride;

} // namespace detail

constexpr struct mdspan_non_standard_tag {
} mdspan_non_standard;

} // namespace MDSPAN_IMPL_STANDARD_NAMESPACE
