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

#ifndef KOKKOS_STD_NUMERICS_EXCLUSIVE_SCAN_HPP
#define KOKKOS_STD_NUMERICS_EXCLUSIVE_SCAN_HPP

#include <Kokkos_Core.hpp>
#include "../Kokkos_BeginEnd.hpp"
#include "../Kokkos_Constraints.hpp"
#include "../Kokkos_Distance.hpp"
#include "../Kokkos_ModifyingOperations.hpp"
#include "../Kokkos_ValueWrapperForNoNeutralElement.hpp"
#include "Kokkos_IdentityReferenceUnaryFunctor.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest>
struct ExclusiveScanDefaultFunctor {
  using execution_space = ExeSpace;
  using value_type =
      ::Kokkos::Experimental::Impl::ValueWrapperForNoNeutralElement<ValueType>;

  ValueType m_init_value;
  FirstFrom m_first_from;
  FirstDest m_first_dest;

  KOKKOS_FUNCTION
  ExclusiveScanDefaultFunctor(ValueType init, FirstFrom first_from,
                              FirstDest first_dest)
      : m_init_value(std::move(init)),
        m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, value_type& update,
                  const bool final_pass) const {
    if (final_pass) {
      if (i == 0) {
        m_first_dest[i] = m_init_value;
      } else {
        m_first_dest[i] = update.val + m_init_value;
      }
    }

    const auto tmp = value_type{m_first_from[i], false};
    this->join(update, tmp);
  }

  KOKKOS_FUNCTION
  void init(value_type& update) const {
    update.val        = {};
    update.is_initial = true;
  }

  KOKKOS_FUNCTION
  void join(volatile value_type& update,
            volatile const value_type& input) const {
    if (update.is_initial) {
      update.val        = input.val;
      update.is_initial = false;
    } else {
      update.val = update.val + input.val;
    }
  }
};

template <class ExeSpace, class IndexType, class ValueType, class FirstFrom,
          class FirstDest, class BinaryOpType, class UnaryOpType>
struct TransformExclusiveScanFunctor {
  using execution_space = ExeSpace;
  using value_type =
      ::Kokkos::Experimental::Impl::ValueWrapperForNoNeutralElement<ValueType>;

  ValueType m_init_value;
  FirstFrom m_first_from;
  FirstDest m_first_dest;
  BinaryOpType m_binary_op;
  UnaryOpType m_unary_op;

  KOKKOS_FUNCTION
  TransformExclusiveScanFunctor(ValueType init, FirstFrom first_from,
                                FirstDest first_dest, BinaryOpType bop,
                                UnaryOpType uop)
      : m_init_value(std::move(init)),
        m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_binary_op(std::move(bop)),
        m_unary_op(std::move(uop)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, value_type& update,
                  const bool final_pass) const {
    if (final_pass) {
      if (i == 0) {
        // for both ExclusiveScan and TransformExclusiveScan,
        // init is unmodified
        m_first_dest[i] = m_init_value;
      } else {
        m_first_dest[i] = m_binary_op(update.val, m_init_value);
      }
    }

    const auto tmp = value_type{m_unary_op(m_first_from[i]), false};
    this->join(update, tmp);
  }

  KOKKOS_FUNCTION
  void init(value_type& update) const {
    update.val        = {};
    update.is_initial = true;
  }

  KOKKOS_FUNCTION
  void join(volatile value_type& update,
            volatile const value_type& input) const {
    if (update.is_initial) {
      update.val = input.val;
    } else {
      update.val = m_binary_op(update.val, input.val);
    }
    update.is_initial = false;
  }
};

// --------------------------------------------------
// exclusive_scan_custom_op_impl
// --------------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
OutputIteratorType exclusive_scan_custom_op_impl(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType first_from, InputIteratorType last_from,
    OutputIteratorType first_dest, ValueType init_value, BinaryOpType bop) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type    = typename InputIteratorType::difference_type;
  using unary_op_type = StdNumericScanIdentityReferenceUnaryFunctor<ValueType>;
  using func_type =
      TransformExclusiveScanFunctor<ExecutionSpace, index_type, ValueType,
                                    InputIteratorType, OutputIteratorType,
                                    BinaryOpType, unary_op_type>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_type(init_value, first_from, first_dest, bop, unary_op_type()));
  ex.fence("Kokkos::exclusive_scan_custom_op: fence after operation");

  // return
  return first_dest + num_elements;
}

// --------------------------------------------------
// transform_exclusive_scan_impl
// --------------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType,
          class UnaryOpType>
OutputIteratorType transform_exclusive_scan_impl(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType first_from, InputIteratorType last_from,
    OutputIteratorType first_dest, ValueType init_value, BinaryOpType bop,
    UnaryOpType uop) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using func_type =
      TransformExclusiveScanFunctor<ExecutionSpace, index_type, ValueType,
                                    InputIteratorType, OutputIteratorType,
                                    BinaryOpType, UnaryOpType>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_type(init_value, first_from, first_dest, bop, uop));
  ex.fence("Kokkos::transform_exclusive_scan: fence after operation");

  // return
  return first_dest + num_elements;
}

// --------------------------------------------------
// exclusive_scan_default_op_impl
// --------------------------------------------------
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
OutputIteratorType exclusive_scan_default_op_impl(const std::string& label,
                                                  const ExecutionSpace& ex,
                                                  InputIteratorType first_from,
                                                  InputIteratorType last_from,
                                                  OutputIteratorType first_dest,
                                                  ValueType init_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // does it make sense to do this static_assert too?
  // using input_iterator_value_type = typename InputIteratorType::value_type;
  // static_assert
  //   (std::is_convertible<std::remove_cv_t<input_iterator_value_type>,
  //   ValueType>::value,
  //    "exclusive_scan: InputIteratorType::value_type not convertible to
  //    ValueType");

  // we are unnecessarily duplicating code, but this is on purpose
  // so that we can use the default_op for OpenMPTarget.
  // Originally, I had this implemented as:
  // '''
  // using bop_type   = StdExclusiveScanDefaultJoinFunctor<ValueType>;
  // call exclusive_scan_custom_op_impl(..., bop_type());
  // '''
  // which avoids duplicating the functors, but for OpenMPTarget
  // I cannot use a custom binary op.
  // This is the same problem that occurs for reductions.

  // aliases
  using index_type = typename InputIteratorType::difference_type;
  using func_type =
      ExclusiveScanDefaultFunctor<ExecutionSpace, index_type, ValueType,
                                  InputIteratorType, OutputIteratorType>;

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_scan(label,
                          RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                          func_type(init_value, first_from, first_dest));
  ex.fence("Kokkos::exclusive_scan_default_op: fence after operation");

  return first_dest + num_elements;
}

}  // end namespace Impl

///////////////////////////////
//
// exclusive scan API
//
///////////////////////////////

// overload set 1
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIteratorType, OutputIteratorType>::value,
                  OutputIteratorType>
exclusive_scan(const ExecutionSpace& ex, InputIteratorType first,
               InputIteratorType last, OutputIteratorType first_dest,
               ValueType init_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_default_op_impl(
      "Kokkos::exclusive_scan_default_functors_iterator_api", ex, first, last,
      first_dest, init_value);
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIteratorType, OutputIteratorType>::value,
                  OutputIteratorType>
exclusive_scan(const std::string& label, const ExecutionSpace& ex,
               InputIteratorType first, InputIteratorType last,
               OutputIteratorType first_dest, ValueType init_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_default_op_impl(label, ex, first, last,
                                              first_dest, init_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto exclusive_scan(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_default_op_impl(
      "Kokkos::exclusive_scan_default_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      init_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto exclusive_scan(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_default_op_impl(label, ex, KE::cbegin(view_from),
                                              KE::cend(view_from),
                                              KE::begin(view_dest), init_value);
}

// overload set 2
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIteratorType, OutputIteratorType>::value,
                  OutputIteratorType>
exclusive_scan(const ExecutionSpace& ex, InputIteratorType first,
               InputIteratorType last, OutputIteratorType first_dest,
               ValueType init_value, BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_custom_op_impl(
      "Kokkos::exclusive_scan_custom_functors_iterator_api", ex, first, last,
      first_dest, init_value, bop);
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIteratorType, OutputIteratorType>::value,
                  OutputIteratorType>
exclusive_scan(const std::string& label, const ExecutionSpace& ex,
               InputIteratorType first, InputIteratorType last,
               OutputIteratorType first_dest, ValueType init_value,
               BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_custom_op_impl(label, ex, first, last, first_dest,
                                             init_value, bop);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryOpType>
auto exclusive_scan(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value, BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_custom_op_impl(
      "Kokkos::exclusive_scan_custom_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      init_value, bop);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryOpType>
auto exclusive_scan(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value, BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_custom_op_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), init_value, bop);
}

//////////////////////////////////////
//
// transform_exclusive_scan public API
//
//////////////////////////////////////

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType,
          class UnaryOpType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIteratorType, OutputIteratorType>::value,
                  OutputIteratorType>
transform_exclusive_scan(const ExecutionSpace& ex, InputIteratorType first,
                         InputIteratorType last, OutputIteratorType first_dest,
                         ValueType init_value, BinaryOpType binary_op,
                         UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::transform_exclusive_scan_impl(
      "Kokkos::transform_exclusive_scan_custom_functors_iterator_api", ex,
      first, last, first_dest, init_value, binary_op, unary_op);
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType,
          class UnaryOpType>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIteratorType, OutputIteratorType>::value,
                  OutputIteratorType>
transform_exclusive_scan(const std::string& label, const ExecutionSpace& ex,
                         InputIteratorType first, InputIteratorType last,
                         OutputIteratorType first_dest, ValueType init_value,
                         BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::transform_exclusive_scan_impl(label, ex, first, last, first_dest,
                                             init_value, binary_op, unary_op);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryOpType, class UnaryOpType>
auto transform_exclusive_scan(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    ValueType init_value, BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_exclusive_scan_impl(
      "Kokkos::transform_exclusive_scan_custom_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      init_value, binary_op, unary_op);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryOpType, class UnaryOpType>
auto transform_exclusive_scan(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    ValueType init_value, BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_exclusive_scan_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), init_value, binary_op, unary_op);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
