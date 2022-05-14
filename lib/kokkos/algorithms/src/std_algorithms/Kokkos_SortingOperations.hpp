/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_STD_SORTING_OPERATIONS_HPP
#define KOKKOS_STD_SORTING_OPERATIONS_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_BeginEnd.hpp"
#include "Kokkos_Constraints.hpp"
#include "Kokkos_NonModifyingSequenceOperations.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

// ------------------
//
// functors
//
// ------------------

template <class IteratorType, class IndicatorViewType, class ComparatorType>
struct StdIsSortedUntilFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  IndicatorViewType m_indicator;
  ComparatorType m_comparator;

  KOKKOS_FUNCTION
  void operator()(const index_type i, int& update, const bool final) const {
    const auto& val_i   = m_first[i];
    const auto& val_ip1 = m_first[i + 1];

    if (m_comparator(val_ip1, val_i)) {
      ++update;
    }

    if (final) {
      m_indicator(i) = update;
    }
  }

  KOKKOS_FUNCTION
  StdIsSortedUntilFunctor(IteratorType _first1, IndicatorViewType indicator,
                          ComparatorType comparator)
      : m_first(std::move(_first1)),
        m_indicator(std::move(indicator)),
        m_comparator(std::move(comparator)) {}
};

template <class IteratorType, class ComparatorType>
struct StdIsSortedFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  ComparatorType m_comparator;

  KOKKOS_FUNCTION
  void operator()(const index_type i, std::size_t& update) const {
    const auto& val_i   = m_first[i];
    const auto& val_ip1 = m_first[i + 1];

    if (m_comparator(val_ip1, val_i)) {
      ++update;
    }
  }

  KOKKOS_FUNCTION
  StdIsSortedFunctor(IteratorType _first1, ComparatorType comparator)
      : m_first(std::move(_first1)), m_comparator(std::move(comparator)) {}
};

// ------------------------------------------
// is_sorted_until_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class ComparatorType>
IteratorType is_sorted_until_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  IteratorType last, ComparatorType comp) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);

  // trivial case
  if (num_elements <= 1) {
    return last;
  }

  /*
    use scan and a helper "indicator" view
    such that we scan the data and fill the indicator with
    partial sum that is always 0 unless we find a pair that
    breaks the sorting, so in that case the indicator will
    have a 1 starting at the location where the sorting breaks.
    So finding that 1 means finding the location we want.
   */

  // aliases
  using indicator_value_type = std::size_t;
  using indicator_view_type =
      ::Kokkos::View<indicator_value_type*, ExecutionSpace>;
  using functor_type =
      StdIsSortedUntilFunctor<IteratorType, indicator_view_type,
                              ComparatorType>;

  // do scan
  // use num_elements-1 because each index handles i and i+1
  const auto num_elements_minus_one = num_elements - 1;
  indicator_view_type indicator("is_sorted_until_indicator_helper",
                                num_elements_minus_one);
  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_minus_one),
      functor_type(first, indicator, std::move(comp)));

  // try to find the first sentinel value, which indicates
  // where the sorting condition breaks
  namespace KE                                  = ::Kokkos::Experimental;
  constexpr indicator_value_type sentinel_value = 1;
  auto r =
      KE::find(ex, KE::cbegin(indicator), KE::cend(indicator), sentinel_value);
  const auto shift = r - ::Kokkos::Experimental::cbegin(indicator);

  return first + (shift + 1);
}

template <class ExecutionSpace, class IteratorType>
IteratorType is_sorted_until_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  IteratorType last) {
  using value_type = typename IteratorType::value_type;
  using pred_t     = Impl::StdAlgoLessThanBinaryPredicate<value_type>;
  return is_sorted_until_impl(label, ex, first, last, pred_t());
}

// ------------------------------------------
// is_sorted_impl
// ------------------------------------------
template <class ExecutionSpace, class IteratorType, class ComparatorType>
bool is_sorted_impl(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last,
                    ComparatorType comp) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);
  if (num_elements <= 1) {
    return true;
  }

  // use num_elements-1 because each index handles i and i+1
  const auto num_elements_minus_one = num_elements - 1;
  using functor_type = StdIsSortedFunctor<IteratorType, ComparatorType>;

  // result is incremented by one if sorting breaks at index i
  std::size_t result = 0;
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_minus_one),
      functor_type(first, std::move(comp)), result);

  return result == 0;
}

template <class ExecutionSpace, class IteratorType>
bool is_sorted_impl(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last) {
  using value_type = typename IteratorType::value_type;
  using pred_t     = Impl::StdAlgoLessThanBinaryPredicate<value_type>;
  return is_sorted_impl(label, ex, first, last, pred_t());
}

}  // namespace Impl

// ----------------------------------
// is_sorted_until public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType>
IteratorType is_sorted_until(const ExecutionSpace& ex, IteratorType first,
                             IteratorType last) {
  return Impl::is_sorted_until_impl(
      "Kokkos::is_sorted_until_iterator_api_default", ex, first, last);
}

template <class ExecutionSpace, class IteratorType>
IteratorType is_sorted_until(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, IteratorType last) {
  return Impl::is_sorted_until_impl(label, ex, first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto is_sorted_until(const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_until_impl("Kokkos::is_sorted_until_view_api_default",
                                    ex, KE::begin(view), KE::end(view));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto is_sorted_until(const std::string& label, const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_until_impl(label, ex, KE::begin(view), KE::end(view));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
IteratorType is_sorted_until(const ExecutionSpace& ex, IteratorType first,
                             IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);
  return Impl::is_sorted_until_impl(
      "Kokkos::is_sorted_until_iterator_api_default", ex, first, last,
      std::move(comp));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
IteratorType is_sorted_until(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, IteratorType last,
                             ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::is_sorted_until_impl(label, ex, first, last, std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ComparatorType>
auto is_sorted_until(const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& view,
                     ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_not_openmptarget(ex);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_until_impl("Kokkos::is_sorted_until_view_api_default",
                                    ex, KE::begin(view), KE::end(view),
                                    std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ComparatorType>
auto is_sorted_until(const std::string& label, const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& view,
                     ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_not_openmptarget(ex);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_until_impl(label, ex, KE::begin(view), KE::end(view),
                                    std::move(comp));
}

// ----------------------------------
// is_sorted public API
// ----------------------------------
template <class ExecutionSpace, class IteratorType>
bool is_sorted(const ExecutionSpace& ex, IteratorType first,
               IteratorType last) {
  return Impl::is_sorted_impl("Kokkos::is_sorted_iterator_api_default", ex,
                              first, last);
}

template <class ExecutionSpace, class IteratorType>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               IteratorType first, IteratorType last) {
  return Impl::is_sorted_impl(label, ex, first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
bool is_sorted(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_impl("Kokkos::is_sorted_view_api_default", ex,
                              KE::cbegin(view), KE::cend(view));
}

template <class ExecutionSpace, class DataType, class... Properties>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_impl(label, ex, KE::cbegin(view), KE::cend(view));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
bool is_sorted(const ExecutionSpace& ex, IteratorType first, IteratorType last,
               ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);
  return Impl::is_sorted_impl("Kokkos::is_sorted_iterator_api_default", ex,
                              first, last, std::move(comp));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               IteratorType first, IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);
  return Impl::is_sorted_impl(label, ex, first, last, std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ComparatorType>
bool is_sorted(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_not_openmptarget(ex);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_impl("Kokkos::is_sorted_view_api_default", ex,
                              KE::cbegin(view), KE::cend(view),
                              std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ComparatorType>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_not_openmptarget(ex);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_impl(label, ex, KE::cbegin(view), KE::cend(view),
                              std::move(comp));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
