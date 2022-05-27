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

#ifndef KOKKOS_STD_MIN_MAX_ELEMENT_OPERATIONS_HPP
#define KOKKOS_STD_MIN_MAX_ELEMENT_OPERATIONS_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_BeginEnd.hpp"
#include "Kokkos_Constraints.hpp"
#include "Kokkos_Distance.hpp"
#include "Kokkos_ModifyingOperations.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ReducerType>
struct StdMinOrMaxElemFunctor {
  using index_type     = typename IteratorType::difference_type;
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  ReducerType m_reducer;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& red_value) const {
    m_reducer.join(red_value, red_value_type{m_first[i], i});
  }

  KOKKOS_FUNCTION
  StdMinOrMaxElemFunctor(IteratorType first, ReducerType reducer)
      : m_first(std::move(first)), m_reducer(std::move(reducer)) {}
};

template <class IteratorType, class ReducerType>
struct StdMinMaxElemFunctor {
  using index_type     = typename IteratorType::difference_type;
  using red_value_type = typename ReducerType::value_type;
  IteratorType m_first;
  ReducerType m_reducer;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& red_value) const {
    const auto& my_value = m_first[i];
    m_reducer.join(red_value, red_value_type{my_value, my_value, i, i});
  }

  KOKKOS_FUNCTION
  StdMinMaxElemFunctor(IteratorType first, ReducerType reducer)
      : m_first(std::move(first)), m_reducer(std::move(reducer)) {}
};

// ------------------------------------------
// min_or_max_element_impl
// ------------------------------------------
template <template <class... Args> class ReducerType, class ExecutionSpace,
          class IteratorType, class... Args>
IteratorType min_or_max_element_impl(const std::string& label,
                                     const ExecutionSpace& ex,
                                     IteratorType first, IteratorType last,
                                     Args&&... args) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return last;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using value_type           = typename IteratorType::value_type;
  using reducer_type         = ReducerType<value_type, index_type, Args...>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdMinOrMaxElemFunctor<IteratorType, reducer_type>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result, std::forward<Args>(args)...);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer), reducer);

  // fence not needed because reducing into scalar

  // return
  return first + red_result.loc;
}

// ------------------------------------------
// minmax_element_impl
// ------------------------------------------
template <template <class... Args> class ReducerType, class ExecutionSpace,
          class IteratorType, class... Args>
::Kokkos::pair<IteratorType, IteratorType> minmax_element_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType first,
    IteratorType last, Args&&... args) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return {first, first};
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using value_type           = typename IteratorType::value_type;
  using reducer_type         = ReducerType<value_type, index_type, Args...>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t               = StdMinMaxElemFunctor<IteratorType, reducer_type>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result, std::forward<Args>(args)...);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer), reducer);

  // fence not needed because reducing into scalar

  // return
  return {first + red_result.min_loc, first + red_result.max_loc};
}

}  // end namespace Impl

// ----------------------
// min_element public API
// ----------------------
template <class ExecutionSpace, class IteratorType>
auto min_element(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last) {
  return Impl::min_or_max_element_impl<MinFirstLoc>(
      "Kokkos::min_element_iterator_api_default", ex, first, last);
}

template <class ExecutionSpace, class IteratorType>
auto min_element(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last) {
  return Impl::min_or_max_element_impl<MinFirstLoc>(label, ex, first, last);
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto min_element(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MinFirstLocCustomComparator>(
      "Kokkos::min_element_iterator_api_default", ex, first, last,
      std::move(comp));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto min_element(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MinFirstLocCustomComparator>(
      label, ex, first, last, std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto min_element(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::min_or_max_element_impl<MinFirstLoc>(
      "Kokkos::min_element_view_api_default", ex, begin(v), end(v));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto min_element(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MinFirstLocCustomComparator>(
      "Kokkos::min_element_view_api_default", ex, begin(v), end(v),
      std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto min_element(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::min_or_max_element_impl<MinFirstLoc>(label, ex, begin(v),
                                                    end(v));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto min_element(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MinFirstLocCustomComparator>(
      label, ex, begin(v), end(v), std::move(comp));
}

// ----------------------
// max_element public API
// ----------------------
template <class ExecutionSpace, class IteratorType>
auto max_element(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last) {
  return Impl::min_or_max_element_impl<MaxFirstLoc>(
      "Kokkos::max_element_iterator_api_default", ex, first, last);
}

template <class ExecutionSpace, class IteratorType>
auto max_element(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last) {
  return Impl::min_or_max_element_impl<MaxFirstLoc>(label, ex, first, last);
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto max_element(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MaxFirstLocCustomComparator>(
      "Kokkos::max_element_iterator_api_default", ex, first, last,
      std::move(comp));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto max_element(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MaxFirstLocCustomComparator>(
      label, ex, first, last, std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto max_element(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::min_or_max_element_impl<MaxFirstLoc>(
      "Kokkos::max_element_view_api_default", ex, begin(v), end(v));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto max_element(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::min_or_max_element_impl<MaxFirstLoc>(label, ex, begin(v),
                                                    end(v));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto max_element(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MaxFirstLocCustomComparator>(
      "Kokkos::max_element_view_api_default", ex, begin(v), end(v),
      std::move(comp));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto max_element(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::min_or_max_element_impl<MaxFirstLocCustomComparator>(
      label, ex, begin(v), end(v), std::move(comp));
}

// -------------------------
// minmax_element public API
// -------------------------
template <class ExecutionSpace, class IteratorType>
auto minmax_element(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last) {
  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(
      "Kokkos::minmax_element_iterator_api_default", ex, first, last);
}

template <class ExecutionSpace, class IteratorType>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last) {
  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(label, ex, first, last);
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto minmax_element(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      "Kokkos::minmax_element_iterator_api_default", ex, first, last,
      std::move(comp));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last,
                    ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      label, ex, first, last, std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto minmax_element(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(
      "Kokkos::minmax_element_view_api_default", ex, begin(v), end(v));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(label, ex, begin(v),
                                                       end(v));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto minmax_element(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      "Kokkos::minmax_element_view_api_default", ex, begin(v), end(v),
      std::move(comp));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      label, ex, begin(v), end(v), std::move(comp));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
