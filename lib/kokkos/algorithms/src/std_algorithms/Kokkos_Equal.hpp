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

#ifndef KOKKOS_STD_ALGORITHMS_EQUAL_HPP
#define KOKKOS_STD_ALGORITHMS_EQUAL_HPP

#include "impl/Kokkos_Equal.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2) {
  return Impl::equal_impl(label, ex, first1, last1, first2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2, BinaryPredicateType predicate) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2, std::move(predicate));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2,
      BinaryPredicateType predicate) {
  return Impl::equal_impl(label, ex, first1, last1, first2,
                          std::move(predicate));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
bool equal(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl("Kokkos::equal_view_api_default", ex,
                          KE::cbegin(view1), KE::cend(view1),
                          KE::cbegin(view2));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
bool equal(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl(label, ex, KE::cbegin(view1), KE::cend(view1),
                          KE::cbegin(view2));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
bool equal(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2,
           BinaryPredicateType predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl("Kokkos::equal_view_api_default", ex,
                          KE::cbegin(view1), KE::cend(view1), KE::cbegin(view2),
                          std::move(predicate));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryPredicateType>
bool equal(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2,
           BinaryPredicateType predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_impl(label, ex, KE::cbegin(view1), KE::cend(view1),
                          KE::cbegin(view2), std::move(predicate));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2, IteratorType2 last2) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2, last2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2, IteratorType2 last2) {
  return Impl::equal_impl(label, ex, first1, last1, first2, last2);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
      IteratorType2 first2, IteratorType2 last2,
      BinaryPredicateType predicate) {
  return Impl::equal_impl("Kokkos::equal_iterator_api_default", ex, first1,
                          last1, first2, last2, std::move(predicate));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      IteratorType1, IteratorType2>::value,
                  bool>
equal(const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
      IteratorType1 last1, IteratorType2 first2, IteratorType2 last2,
      BinaryPredicateType predicate) {
  return Impl::equal_impl(label, ex, first1, last1, first2, last2,
                          std::move(predicate));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
