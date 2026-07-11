// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_BEGIN_END_HPP
#define KOKKOS_BEGIN_END_HPP

#include "impl/Kokkos_RandomAccessIterator.hpp"
#include "impl/Kokkos_Constraints.hpp"

/// \file Kokkos_BeginEnd.hpp
/// \brief Kokkos begin, end, cbegin, cend

namespace Kokkos {
namespace Experimental {

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto begin(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using it_t =
      Impl::RandomAccessIterator<Kokkos::View<DataType, Properties...>>;
  return it_t(v);
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto end(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using it_t =
      Impl::RandomAccessIterator<Kokkos::View<DataType, Properties...>>;
  return it_t(v, v.extent(0));
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto cbegin(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using ViewConstType =
      typename Kokkos::View<DataType, Properties...>::const_type;
  const ViewConstType cv = v;
  using it_t             = Impl::RandomAccessIterator<ViewConstType>;
  return it_t(cv);
}

template <class DataType, class... Properties>
KOKKOS_INLINE_FUNCTION auto cend(
    const Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  using ViewConstType =
      typename Kokkos::View<DataType, Properties...>::const_type;
  const ViewConstType cv = v;
  using it_t             = Impl::RandomAccessIterator<ViewConstType>;
  return it_t(cv, cv.extent(0));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
