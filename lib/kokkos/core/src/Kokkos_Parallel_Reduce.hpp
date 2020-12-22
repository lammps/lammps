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

#ifndef KOKKOS_PARALLEL_REDUCE_HPP
#define KOKKOS_PARALLEL_REDUCE_HPP

#include <Kokkos_NumericTraits.hpp>
#include <Kokkos_View.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>
#include <type_traits>

namespace Kokkos {

template <class T, class Enable = void>
struct is_reducer_type {
  enum { value = 0 };
};

template <class T>
struct is_reducer_type<
    T, typename std::enable_if<std::is_same<
           typename std::remove_cv<T>::type,
           typename std::remove_cv<typename T::reducer>::type>::value>::type> {
  enum { value = 1 };
};

template <class Scalar, class Space>
struct Sum {
 public:
  // Required
  using reducer    = Sum<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  Sum(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  Sum(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const { dest += src; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest += src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::sum();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Space>
struct Prod {
 public:
  // Required
  using reducer    = Prod<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  Prod(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  Prod(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const { dest *= src; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest *= src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::prod();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Space>
struct Min {
 public:
  // Required
  using reducer    = Min<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  Min(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  Min(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src < dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if (src < dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Space>
struct Max {
 public:
  // Required
  using reducer    = Max<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  Max(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  Max(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src > dest) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if (src > dest) dest = src;
  }

  // Required
  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::max();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Space>
struct LAnd {
 public:
  // Required
  using reducer    = LAnd<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  LAnd(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  LAnd(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::land();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Space>
struct LOr {
 public:
  // Required
  using reducer    = LOr<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  LOr(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  LOr(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::lor();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Space>
struct BAnd {
 public:
  // Required
  using reducer    = BAnd<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  BAnd(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  BAnd(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::band();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Space>
struct BOr {
 public:
  // Required
  using reducer    = BOr<Scalar, Space>;
  using value_type = typename std::remove_cv<Scalar>::type;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  BOr(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  BOr(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val = reduction_identity<value_type>::bor();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Index>
struct ValLocScalar {
  Scalar val;
  Index loc;

  KOKKOS_INLINE_FUNCTION
  void operator=(const ValLocScalar& rhs) {
    val = rhs.val;
    loc = rhs.loc;
  }

  KOKKOS_INLINE_FUNCTION
  void operator=(const volatile ValLocScalar& rhs) volatile {
    val = rhs.val;
    loc = rhs.loc;
  }
};

template <class Scalar, class Index, class Space>
struct MinLoc {
 private:
  using scalar_type = typename std::remove_cv<Scalar>::type;
  using index_type  = typename std::remove_cv<Index>::type;

 public:
  // Required
  using reducer    = MinLoc<Scalar, Index, Space>;
  using value_type = ValLocScalar<scalar_type, index_type>;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  MinLoc(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  MinLoc(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.val < dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if (src.val < dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::min();
    val.loc = reduction_identity<index_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Index, class Space>
struct MaxLoc {
 private:
  using scalar_type = typename std::remove_cv<Scalar>::type;
  using index_type  = typename std::remove_cv<Index>::type;

 public:
  // Required
  using reducer    = MaxLoc<Scalar, Index, Space>;
  using value_type = ValLocScalar<scalar_type, index_type>;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  MaxLoc(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  MaxLoc(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.val > dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if (src.val > dest.val) dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.val = reduction_identity<scalar_type>::max();
    val.loc = reduction_identity<index_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar>
struct MinMaxScalar {
  Scalar min_val, max_val;

  KOKKOS_INLINE_FUNCTION
  void operator=(const MinMaxScalar& rhs) {
    min_val = rhs.min_val;
    max_val = rhs.max_val;
  }

  KOKKOS_INLINE_FUNCTION
  void operator=(const volatile MinMaxScalar& rhs) volatile {
    min_val = rhs.min_val;
    max_val = rhs.max_val;
  }
};

template <class Scalar, class Space>
struct MinMax {
 private:
  using scalar_type = typename std::remove_cv<Scalar>::type;

 public:
  // Required
  using reducer    = MinMax<Scalar, Space>;
  using value_type = MinMaxScalar<scalar_type>;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  MinMax(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  MinMax(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_val = reduction_identity<scalar_type>::max();
    val.min_val = reduction_identity<scalar_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};

template <class Scalar, class Index>
struct MinMaxLocScalar {
  Scalar min_val, max_val;
  Index min_loc, max_loc;

  KOKKOS_INLINE_FUNCTION
  void operator=(const MinMaxLocScalar& rhs) {
    min_val = rhs.min_val;
    min_loc = rhs.min_loc;
    max_val = rhs.max_val;
    max_loc = rhs.max_loc;
  }

  KOKKOS_INLINE_FUNCTION
  void operator=(const volatile MinMaxLocScalar& rhs) volatile {
    min_val = rhs.min_val;
    min_loc = rhs.min_loc;
    max_val = rhs.max_val;
    max_loc = rhs.max_loc;
  }
};

template <class Scalar, class Index, class Space>
struct MinMaxLoc {
 private:
  using scalar_type = typename std::remove_cv<Scalar>::type;
  using index_type  = typename std::remove_cv<Index>::type;

 public:
  // Required
  using reducer    = MinMaxLoc<Scalar, Index, Space>;
  using value_type = MinMaxLocScalar<scalar_type, index_type>;

  using result_view_type = Kokkos::View<value_type, Space>;

 private:
  result_view_type value;
  bool references_scalar_v;

 public:
  KOKKOS_INLINE_FUNCTION
  MinMaxLoc(value_type& value_) : value(&value_), references_scalar_v(true) {}

  KOKKOS_INLINE_FUNCTION
  MinMaxLoc(const result_view_type& value_)
      : value(value_), references_scalar_v(false) {}

  // Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if (src.min_val < dest.min_val) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if (src.max_val > dest.max_val) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const {
    val.max_val = reduction_identity<scalar_type>::max();
    val.min_val = reduction_identity<scalar_type>::min();
    val.max_loc = reduction_identity<index_type>::min();
    val.min_loc = reduction_identity<index_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return *value.data(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return value; }

  KOKKOS_INLINE_FUNCTION
  bool references_scalar() const { return references_scalar_v; }
};
}  // namespace Kokkos
namespace Kokkos {
namespace Impl {

template <class T, class ReturnType, class ValueTraits>
struct ParallelReduceReturnValue;

template <class ReturnType, class FunctorType>
struct ParallelReduceReturnValue<
    typename std::enable_if<Kokkos::is_view<ReturnType>::value>::type,
    ReturnType, FunctorType> {
  using return_type  = ReturnType;
  using reducer_type = InvalidType;

  using value_type_scalar = typename return_type::value_type;
  using value_type_array  = typename return_type::value_type* const;

  using value_type = typename if_c<return_type::rank == 0, value_type_scalar,
                                   value_type_array>::type;

  static return_type& return_value(ReturnType& return_val, const FunctorType&) {
    return return_val;
  }
};

template <class ReturnType, class FunctorType>
struct ParallelReduceReturnValue<
    typename std::enable_if<!Kokkos::is_view<ReturnType>::value &&
                            (!std::is_array<ReturnType>::value &&
                             !std::is_pointer<ReturnType>::value) &&
                            !Kokkos::is_reducer_type<ReturnType>::value>::type,
    ReturnType, FunctorType> {
  using return_type =
      Kokkos::View<ReturnType, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

  using reducer_type = InvalidType;

  using value_type = typename return_type::value_type;

  static return_type return_value(ReturnType& return_val, const FunctorType&) {
    return return_type(&return_val);
  }
};

template <class ReturnType, class FunctorType>
struct ParallelReduceReturnValue<
    typename std::enable_if<(std::is_array<ReturnType>::value ||
                             std::is_pointer<ReturnType>::value)>::type,
    ReturnType, FunctorType> {
  using return_type = Kokkos::View<typename std::remove_const<ReturnType>::type,
                                   Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

  using reducer_type = InvalidType;

  using value_type = typename return_type::value_type[];

  static return_type return_value(ReturnType& return_val,
                                  const FunctorType& functor) {
    if (std::is_array<ReturnType>::value)
      return return_type(return_val);
    else
      return return_type(return_val, functor.value_count);
  }
};

template <class ReturnType, class FunctorType>
struct ParallelReduceReturnValue<
    typename std::enable_if<Kokkos::is_reducer_type<ReturnType>::value>::type,
    ReturnType, FunctorType> {
  using return_type  = ReturnType;
  using reducer_type = ReturnType;
  using value_type   = typename return_type::value_type;

  static return_type return_value(ReturnType& return_val, const FunctorType&) {
    return return_val;
  }
};

template <class T, class ReturnType, class FunctorType>
struct ParallelReducePolicyType;

template <class PolicyType, class FunctorType>
struct ParallelReducePolicyType<
    typename std::enable_if<
        Kokkos::Impl::is_execution_policy<PolicyType>::value>::type,
    PolicyType, FunctorType> {
  using policy_type = PolicyType;
  static PolicyType policy(const PolicyType& policy_) { return policy_; }
};

template <class PolicyType, class FunctorType>
struct ParallelReducePolicyType<
    typename std::enable_if<std::is_integral<PolicyType>::value>::type,
    PolicyType, FunctorType> {
  using execution_space =
      typename Impl::FunctorPolicyExecutionSpace<FunctorType,
                                                 void>::execution_space;

  using policy_type = Kokkos::RangePolicy<execution_space>;

  static policy_type policy(const PolicyType& policy_) {
    return policy_type(0, policy_);
  }
};

template <class FunctorType, class ExecPolicy, class ValueType,
          class ExecutionSpace>
struct ParallelReduceFunctorType {
  using functor_type = FunctorType;
  static const functor_type& functor(const functor_type& functor) {
    return functor;
  }
};

template <class PolicyType, class FunctorType, class ReturnType>
struct ParallelReduceAdaptor {
  using return_value_adapter =
      Impl::ParallelReduceReturnValue<void, ReturnType, FunctorType>;
#ifdef KOKKOS_IMPL_NEED_FUNCTOR_WRAPPER
  using functor_adaptor =
      Impl::ParallelReduceFunctorType<FunctorType, PolicyType,
                                      typename return_value_adapter::value_type,
                                      typename PolicyType::execution_space>;
#endif
  static inline void execute(const std::string& label, const PolicyType& policy,
                             const FunctorType& functor,
                             ReturnType& return_value) {
    uint64_t kpID = 0;

    PolicyType inner_policy = policy;
    Kokkos::Tools::Impl::begin_parallel_reduce<
        typename return_value_adapter::reducer_type>(inner_policy, functor,
                                                     label, kpID);

    Kokkos::Impl::shared_allocation_tracking_disable();
#ifdef KOKKOS_IMPL_NEED_FUNCTOR_WRAPPER
    Impl::ParallelReduce<typename functor_adaptor::functor_type, PolicyType,
                         typename return_value_adapter::reducer_type>
        closure(functor_adaptor::functor(functor), inner_policy,
                return_value_adapter::return_value(return_value, functor));
#else
    Impl::ParallelReduce<FunctorType, PolicyType,
                         typename return_value_adapter::reducer_type>
        closure(functor, inner_policy,
                return_value_adapter::return_value(return_value, functor));
#endif
    Kokkos::Impl::shared_allocation_tracking_enable();
    closure.execute();

    Kokkos::Tools::Impl::end_parallel_reduce<
        typename return_value_adapter::reducer_type>(inner_policy, functor,
                                                     label, kpID);
  }
};
}  // namespace Impl

//----------------------------------------------------------------------------

/*! \fn void parallel_reduce(label,policy,functor,return_argument)
    \brief Perform a parallel reduction.
    \param label An optional Label giving the call name. Must be able to
   construct a std::string from the argument. \param policy A Kokkos Execution
   Policy, such as an integer, a RangePolicy or a TeamPolicy. \param functor A
   functor with a reduction operator, and optional init, join and final
   functions. \param return_argument A return argument which can be a scalar, a
   View, or a ReducerStruct. This argument can be left out if the functor has a
   final function.
*/

// Parallel Reduce Blocking behavior

namespace Impl {
template <typename T>
struct ReducerHasTestReferenceFunction {
  template <typename E>
  static std::true_type test_func(decltype(&E::references_scalar));
  template <typename E>
  static std::false_type test_func(...);

  enum {
    value = std::is_same<std::true_type, decltype(test_func<T>(nullptr))>::value
  };
};

template <class ExecutionSpace, class Arg>
constexpr std::enable_if_t<
    // constraints only necessary because SFINAE lacks subsumption
    !ReducerHasTestReferenceFunction<Arg>::value &&
        !Kokkos::is_view<Arg>::value,
    // return type:
    bool>
parallel_reduce_needs_fence(ExecutionSpace const&, Arg const&) {
  return true;
}

template <class ExecutionSpace, class Reducer>
constexpr std::enable_if_t<
    // equivalent to:
    // (requires (Reducer const& r) {
    //   { reducer.references_scalar() } -> std::convertible_to<bool>;
    // })
    ReducerHasTestReferenceFunction<Reducer>::value,
    // return type:
    bool>
parallel_reduce_needs_fence(ExecutionSpace const&, Reducer const& reducer) {
  return reducer.references_scalar();
}

template <class ExecutionSpace, class ViewLike>
constexpr std::enable_if_t<
    // requires Kokkos::ViewLike<ViewLike>
    Kokkos::is_view<ViewLike>::value,
    // return type:
    bool>
parallel_reduce_needs_fence(ExecutionSpace const&, ViewLike const&) {
  return false;
}

template <class ExecutionSpace, class... Args>
struct ParallelReduceFence {
  template <class... ArgsDeduced>
  static void fence(const ExecutionSpace& ex, ArgsDeduced&&... args) {
    if (Impl::parallel_reduce_needs_fence(ex, (ArgsDeduced &&) args...)) {
      ex.fence();
    }
  }
};

}  // namespace Impl

/** \brief  Parallel reduction
 *
 * parallel_reduce performs parallel reductions with arbitrary functions - i.e.
 * it is not solely data based. The call expects up to 4 arguments:
 *
 *
 * Example of a parallel_reduce functor for a POD (plain old data) value type:
 * \code
 *  class FunctorType { // For POD value type
 *  public:
 *    using execution_space = ...;
 *    using value_type = <podType>;
 *    void operator()( <intType> iwork , <podType> & update ) const ;
 *    void init( <podType> & update ) const ;
 *    void join( volatile       <podType> & update ,
 *               volatile const <podType> & input ) const ;
 *
 *    using has_final = true_type;
 *    void final( <podType> & update ) const ;
 *  };
 * \endcode
 *
 * Example of a parallel_reduce functor for an array of POD (plain old data)
 * values:
 * \code
 *  class FunctorType { // For array of POD value
 *  public:
 *    using execution_space = ...;
 *    using value_type = <podType>[];
 *    void operator()( <intType> , <podType> update[] ) const ;
 *    void init( <podType> update[] ) const ;
 *    void join( volatile       <podType> update[] ,
 *               volatile const <podType> input[] ) const ;
 *
 *    using has_final = true_type;
 *    void final( <podType> update[] ) const ;
 *  };
 * \endcode
 */

// ReturnValue is scalar or array: take by reference

template <class PolicyType, class FunctorType, class ReturnType>
inline typename std::enable_if<
    Kokkos::Impl::is_execution_policy<PolicyType>::value>::type
parallel_reduce(const std::string& label, const PolicyType& policy,
                const FunctorType& functor, ReturnType& return_value) {
  Impl::ParallelReduceAdaptor<PolicyType, FunctorType, ReturnType>::execute(
      label, policy, functor, return_value);
  Impl::ParallelReduceFence<typename PolicyType::execution_space,
                            ReturnType>::fence(policy.space(), return_value);
}

template <class PolicyType, class FunctorType, class ReturnType>
inline typename std::enable_if<
    Kokkos::Impl::is_execution_policy<PolicyType>::value>::type
parallel_reduce(const PolicyType& policy, const FunctorType& functor,
                ReturnType& return_value) {
  Impl::ParallelReduceAdaptor<PolicyType, FunctorType, ReturnType>::execute(
      "", policy, functor, return_value);
  Impl::ParallelReduceFence<typename PolicyType::execution_space,
                            ReturnType>::fence(policy.space(), return_value);
}

template <class FunctorType, class ReturnType>
inline void parallel_reduce(const size_t& policy, const FunctorType& functor,
                            ReturnType& return_value) {
  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  Impl::ParallelReduceAdaptor<policy_type, FunctorType, ReturnType>::execute(
      "", policy_type(0, policy), functor, return_value);
  Impl::ParallelReduceFence<typename policy_type::execution_space, ReturnType>::
      fence(typename policy_type::execution_space(), return_value);
}

template <class FunctorType, class ReturnType>
inline void parallel_reduce(const std::string& label, const size_t& policy,
                            const FunctorType& functor,
                            ReturnType& return_value) {
  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  Impl::ParallelReduceAdaptor<policy_type, FunctorType, ReturnType>::execute(
      label, policy_type(0, policy), functor, return_value);
  Impl::ParallelReduceFence<typename policy_type::execution_space, ReturnType>::
      fence(typename policy_type::execution_space(), return_value);
}

// ReturnValue as View or Reducer: take by copy to allow for inline construction

template <class PolicyType, class FunctorType, class ReturnType>
inline typename std::enable_if<
    Kokkos::Impl::is_execution_policy<PolicyType>::value>::type
parallel_reduce(const std::string& label, const PolicyType& policy,
                const FunctorType& functor, const ReturnType& return_value) {
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<PolicyType, FunctorType, ReturnType>::execute(
      label, policy, functor, return_value_impl);
  Impl::ParallelReduceFence<typename PolicyType::execution_space,
                            ReturnType>::fence(policy.space(), return_value);
}

template <class PolicyType, class FunctorType, class ReturnType>
inline typename std::enable_if<
    Kokkos::Impl::is_execution_policy<PolicyType>::value>::type
parallel_reduce(const PolicyType& policy, const FunctorType& functor,
                const ReturnType& return_value) {
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<PolicyType, FunctorType, ReturnType>::execute(
      "", policy, functor, return_value_impl);
  Impl::ParallelReduceFence<typename PolicyType::execution_space,
                            ReturnType>::fence(policy.space(), return_value);
}

template <class FunctorType, class ReturnType>
inline void parallel_reduce(const size_t& policy, const FunctorType& functor,
                            const ReturnType& return_value) {
  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<policy_type, FunctorType, ReturnType>::execute(
      "", policy_type(0, policy), functor, return_value_impl);
  Impl::ParallelReduceFence<typename policy_type::execution_space, ReturnType>::
      fence(typename policy_type::execution_space(), return_value);
}

template <class FunctorType, class ReturnType>
inline void parallel_reduce(const std::string& label, const size_t& policy,
                            const FunctorType& functor,
                            const ReturnType& return_value) {
  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<policy_type, FunctorType, ReturnType>::execute(
      label, policy_type(0, policy), functor, return_value_impl);
  Impl::ParallelReduceFence<typename policy_type::execution_space, ReturnType>::
      fence(typename policy_type::execution_space(), return_value);
}

// No Return Argument

template <class PolicyType, class FunctorType>
inline void parallel_reduce(
    const std::string& label, const PolicyType& policy,
    const FunctorType& functor,
    typename std::enable_if<
        Kokkos::Impl::is_execution_policy<PolicyType>::value>::type* =
        nullptr) {
  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, void>;
  using value_type =
      typename Kokkos::Impl::if_c<(ValueTraits::StaticValueSize != 0),
                                  typename ValueTraits::value_type,
                                  typename ValueTraits::pointer_type>::type;

  static_assert(
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE, PolicyType,
                            FunctorType>::has_final_member_function,
      "Calling parallel_reduce without either return value or final function.");

  using result_view_type =
      Kokkos::View<value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  result_view_type result_view;

  Impl::ParallelReduceAdaptor<PolicyType, FunctorType,
                              result_view_type>::execute(label, policy, functor,
                                                         result_view);
}

template <class PolicyType, class FunctorType>
inline void parallel_reduce(
    const PolicyType& policy, const FunctorType& functor,
    typename std::enable_if<
        Kokkos::Impl::is_execution_policy<PolicyType>::value>::type* =
        nullptr) {
  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, void>;
  using value_type =
      typename Kokkos::Impl::if_c<(ValueTraits::StaticValueSize != 0),
                                  typename ValueTraits::value_type,
                                  typename ValueTraits::pointer_type>::type;

  static_assert(
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE, PolicyType,
                            FunctorType>::has_final_member_function,
      "Calling parallel_reduce without either return value or final function.");

  using result_view_type =
      Kokkos::View<value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  result_view_type result_view;

  Impl::ParallelReduceAdaptor<PolicyType, FunctorType,
                              result_view_type>::execute("", policy, functor,
                                                         result_view);
}

template <class FunctorType>
inline void parallel_reduce(const size_t& policy, const FunctorType& functor) {
  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, void>;
  using value_type =
      typename Kokkos::Impl::if_c<(ValueTraits::StaticValueSize != 0),
                                  typename ValueTraits::value_type,
                                  typename ValueTraits::pointer_type>::type;

  static_assert(
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                            RangePolicy<>,
                            FunctorType>::has_final_member_function,
      "Calling parallel_reduce without either return value or final function.");

  using result_view_type =
      Kokkos::View<value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  result_view_type result_view;

  Impl::ParallelReduceAdaptor<policy_type, FunctorType,
                              result_view_type>::execute("",
                                                         policy_type(0, policy),
                                                         functor, result_view);
}

template <class FunctorType>
inline void parallel_reduce(const std::string& label, const size_t& policy,
                            const FunctorType& functor) {
  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, void>;
  using value_type =
      typename Kokkos::Impl::if_c<(ValueTraits::StaticValueSize != 0),
                                  typename ValueTraits::value_type,
                                  typename ValueTraits::pointer_type>::type;

  static_assert(
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                            RangePolicy<>,
                            FunctorType>::has_final_member_function,
      "Calling parallel_reduce without either return value or final function.");

  using result_view_type =
      Kokkos::View<value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  result_view_type result_view;

  Impl::ParallelReduceAdaptor<policy_type, FunctorType,
                              result_view_type>::execute(label,
                                                         policy_type(0, policy),
                                                         functor, result_view);
}

}  // namespace Kokkos

#endif  // KOKKOS_PARALLEL_REDUCE_HPP
