/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_PARALLEL_REDUCE_HPP
#define KOKKOS_PARALLEL_REDUCE_HPP

#include <Kokkos_NumericTraits.hpp>

namespace Kokkos {

template<class T, class Enable = void>
struct is_reducer_type {
  enum { value = 0 };
};


template<class T>
struct is_reducer_type<T,typename std::enable_if<
                       std::is_same<typename std::remove_cv<T>::type,
                                    typename std::remove_cv<typename T::reducer>::type>::value
                      >::type> {
  enum { value = 1 };
};

namespace Experimental {


template<class Scalar, class Space>
struct Sum {
public:
  //Required
  typedef Sum reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  Sum(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  Sum(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest += src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest += src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::sum();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Space>
struct Prod {
public:
  //Required
  typedef Prod reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  Prod(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  Prod(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest *= src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest *= src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::prod();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Space>
struct Min {
public:
  //Required
  typedef Min reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  Min(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  Min(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src < dest )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src < dest )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Space>
struct Max {
public:
  //Required
  typedef Max reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  Max(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  Max(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src > dest )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src > dest )
      dest = src;
  }

  //Required
  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::max();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Space>
struct LAnd {
public:
  //Required
  typedef LAnd reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  LAnd(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  LAnd(const result_view_type& value_): value(value_.data()) {}

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest && src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::land();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Space>
struct LOr {
public:
  //Required
  typedef LOr reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  LOr(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  LOr(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest || src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::lor();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Space>
struct BAnd {
public:
  //Required
  typedef BAnd reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  BAnd(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  BAnd(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
      dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest & src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::band();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Space>
struct BOr {
public:
  //Required
  typedef BOr reducer;
  typedef typename std::remove_cv<Scalar>::type value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  BOr(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  BOr(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
      dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest = dest | src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val = reduction_identity<value_type>::bor();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Index>
struct ValLocScalar {
  Scalar val;
  Index loc;

  KOKKOS_INLINE_FUNCTION
  void operator = (const ValLocScalar& rhs) {
    val = rhs.val;
    loc = rhs.loc;
  }

  KOKKOS_INLINE_FUNCTION
  void operator = (const volatile ValLocScalar& rhs) volatile {
    val = rhs.val;
    loc = rhs.loc;
  }
};

template<class Scalar, class Index, class Space>
struct MinLoc {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;
  typedef typename std::remove_cv<Index>::type index_type;

public:
  //Required
  typedef MinLoc reducer;
  typedef ValLocScalar<scalar_type,index_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  MinLoc(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  MinLoc(const result_view_type& value_): value(value_.data()) {}


  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.val < dest.val )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.val < dest.val )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.val = reduction_identity<scalar_type>::min();
    val.loc = reduction_identity<index_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Index, class Space>
struct MaxLoc {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;
  typedef typename std::remove_cv<Index>::type index_type;

public:
  //Required
  typedef MaxLoc reducer;
  typedef ValLocScalar<scalar_type,index_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  MaxLoc(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  MaxLoc(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.val > dest.val )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.val > dest.val )
      dest = src;
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.val = reduction_identity<scalar_type>::max();;
    val.loc = reduction_identity<index_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar>
struct MinMaxScalar {
  Scalar min_val,max_val;

  KOKKOS_INLINE_FUNCTION
  void operator = (const MinMaxScalar& rhs) {
    min_val = rhs.min_val;
    max_val = rhs.max_val;
  }

  KOKKOS_INLINE_FUNCTION
  void operator = (const volatile MinMaxScalar& rhs) volatile {
    min_val = rhs.min_val;
    max_val = rhs.max_val;
  }
};

template<class Scalar, class Space>
struct MinMax {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;

public:
  //Required
  typedef MinMax reducer;
  typedef MinMaxScalar<scalar_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  MinMax(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  MinMax(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.max_val = reduction_identity<scalar_type>::max();;
    val.min_val = reduction_identity<scalar_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};

template<class Scalar, class Index>
struct MinMaxLocScalar {
  Scalar min_val,max_val;
  Index min_loc,max_loc;

  KOKKOS_INLINE_FUNCTION
  void operator = (const MinMaxLocScalar& rhs) {
    min_val = rhs.min_val;
    min_loc = rhs.min_loc;
    max_val = rhs.max_val;
    max_loc = rhs.max_loc;
  }

  KOKKOS_INLINE_FUNCTION
  void operator = (const volatile MinMaxLocScalar& rhs) volatile {
    min_val = rhs.min_val;
    min_loc = rhs.min_loc;
    max_val = rhs.max_val;
    max_loc = rhs.max_loc;
  }
};

template<class Scalar, class Index, class Space>
struct MinMaxLoc {
private:
  typedef typename std::remove_cv<Scalar>::type scalar_type;
  typedef typename std::remove_cv<Index>::type index_type;

public:
  //Required
  typedef MinMaxLoc reducer;
  typedef MinMaxLocScalar<scalar_type,index_type> value_type;

  typedef Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  value_type* value;

public:

  KOKKOS_INLINE_FUNCTION
  MinMaxLoc(value_type& value_): value(&value_) {}

  KOKKOS_INLINE_FUNCTION
  MinMaxLoc(const result_view_type& value_): value(value_.data()) {}

  //Required
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src)  const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    if ( src.min_val < dest.min_val ) {
      dest.min_val = src.min_val;
      dest.min_loc = src.min_loc;
    }
    if ( src.max_val > dest.max_val ) {
      dest.max_val = src.max_val;
      dest.max_loc = src.max_loc;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type& val)  const {
    val.max_val = reduction_identity<scalar_type>::max();;
    val.min_val = reduction_identity<scalar_type>::min();
    val.max_loc = reduction_identity<index_type>::min();
    val.min_loc = reduction_identity<index_type>::min();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() {
    return *value;
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result_view_type(value);
  }
};
}
}


namespace Kokkos {
namespace Impl {

template< class T, class ReturnType , class ValueTraits>
struct ParallelReduceReturnValue;

template< class ReturnType , class FunctorType >
struct ParallelReduceReturnValue<typename std::enable_if<Kokkos::is_view<ReturnType>::value>::type, ReturnType, FunctorType> {
  typedef ReturnType return_type;
  typedef InvalidType reducer_type;

  typedef typename return_type::value_type value_type_scalar;
  typedef typename return_type::value_type* const value_type_array;

  typedef typename if_c<return_type::rank==0,value_type_scalar,value_type_array>::type value_type;

  static return_type& return_value(ReturnType& return_val, const FunctorType&) {
    return return_val;
  }
};

template< class ReturnType , class FunctorType>
struct ParallelReduceReturnValue<typename std::enable_if<
                                   !Kokkos::is_view<ReturnType>::value &&
                                  (!std::is_array<ReturnType>::value && !std::is_pointer<ReturnType>::value) &&
                                   !Kokkos::is_reducer_type<ReturnType>::value
                                 >::type, ReturnType, FunctorType> {
  typedef Kokkos::View<  ReturnType
                       , Kokkos::HostSpace
                       , Kokkos::MemoryUnmanaged
      > return_type;

  typedef InvalidType reducer_type;

  typedef typename return_type::value_type value_type;

  static return_type return_value(ReturnType& return_val, const FunctorType&) {
    return return_type(&return_val);
  }
};

template< class ReturnType , class FunctorType>
struct ParallelReduceReturnValue<typename std::enable_if<
                                  (is_array<ReturnType>::value || std::is_pointer<ReturnType>::value)
                                >::type, ReturnType, FunctorType> {
  typedef Kokkos::View<  typename std::remove_const<ReturnType>::type
                       , Kokkos::HostSpace
                       , Kokkos::MemoryUnmanaged
      > return_type;

  typedef InvalidType reducer_type;

  typedef typename return_type::value_type value_type[];

  static return_type return_value(ReturnType& return_val,
                                  const FunctorType& functor) {
    return return_type(return_val,functor.value_count);
  }
};

template< class ReturnType , class FunctorType>
struct ParallelReduceReturnValue<typename std::enable_if<
                                   Kokkos::is_reducer_type<ReturnType>::value
                                >::type, ReturnType, FunctorType> {
  typedef ReturnType return_type;
  typedef ReturnType reducer_type;
  typedef typename return_type::value_type value_type;

  static return_type return_value(ReturnType& return_val,
                                  const FunctorType& functor) {
    return return_val;
  }
};
}

namespace Impl {
template< class T, class ReturnType , class FunctorType>
struct ParallelReducePolicyType;

template< class PolicyType , class FunctorType >
struct ParallelReducePolicyType<typename std::enable_if<Kokkos::Impl::is_execution_policy<PolicyType>::value>::type, PolicyType,FunctorType> {

  typedef PolicyType policy_type;
  static PolicyType policy(const PolicyType& policy_) {
    return policy_;
  }
};

template< class PolicyType , class FunctorType >
struct ParallelReducePolicyType<typename std::enable_if<std::is_integral<PolicyType>::value>::type, PolicyType,FunctorType> {
  typedef typename
    Impl::FunctorPolicyExecutionSpace< FunctorType , void >::execution_space
      execution_space ;

  typedef Kokkos::RangePolicy<execution_space> policy_type;

  static policy_type policy(const PolicyType& policy_) {
    return policy_type(0,policy_);
  }
};

}

namespace Impl {
  template< class FunctorType, class ExecPolicy, class ValueType, class ExecutionSpace>
  struct ParallelReduceFunctorType {
    typedef FunctorType functor_type;
    static const functor_type& functor(const functor_type& functor) {
      return functor;
    }
  };
}

namespace Impl {

  template< class PolicyType, class FunctorType, class ReturnType >
  struct ParallelReduceAdaptor {
    typedef Impl::ParallelReduceReturnValue<void,ReturnType,FunctorType> return_value_adapter;
    #ifdef KOKKOS_IMPL_NEED_FUNCTOR_WRAPPER
    typedef Impl::ParallelReduceFunctorType<FunctorType,PolicyType,
                                            typename return_value_adapter::value_type,
                                            typename PolicyType::execution_space> functor_adaptor;
    #endif
    static inline
    void execute(const std::string& label,
        const PolicyType& policy,
        const FunctorType& functor,
        ReturnType& return_value) {
          #if defined(KOKKOS_ENABLE_PROFILING)
          uint64_t kpID = 0;
          if(Kokkos::Profiling::profileLibraryLoaded()) {
            Kokkos::Impl::ParallelConstructName<FunctorType, typename PolicyType::work_tag> name(label);
            Kokkos::Profiling::beginParallelReduce(name.get(), 0, &kpID);
          }
          #endif

          Kokkos::Impl::shared_allocation_tracking_disable();
          #ifdef KOKKOS_IMPL_NEED_FUNCTOR_WRAPPER
          Impl::ParallelReduce<typename functor_adaptor::functor_type, PolicyType, typename return_value_adapter::reducer_type >
             closure(functor_adaptor::functor(functor),
                     policy,
                     return_value_adapter::return_value(return_value,functor));
          #else
          Impl::ParallelReduce<FunctorType, PolicyType, typename return_value_adapter::reducer_type >
             closure(functor,
                     policy,
                     return_value_adapter::return_value(return_value,functor));
          #endif
          Kokkos::Impl::shared_allocation_tracking_enable();
          closure.execute();

          #if defined(KOKKOS_ENABLE_PROFILING)
          if(Kokkos::Profiling::profileLibraryLoaded()) {
            Kokkos::Profiling::endParallelReduce(kpID);
          }
          #endif
        }

  };
}

//----------------------------------------------------------------------------

#if 0

//----------------------------------------------------------------------------

namespace Impl {

template< class OutType , class InType >
inline
typename std::enable_if
  < std::is_same< OutType , InType >::value , InType const & >::type
forward_execution_policy( InType const & p ) { return p ; }

template< class OutType , class InType >
inline
typename std::enable_if
  < ! std::is_same< OutType , InType >::value , OutType >::type
forward_execution_policy( InType const & p ) { return OutType(p); }


template< class OutType , class InType >
inline
typename std::enable_if
  < std::is_same< OutType , InType >::value , InType const & >::type
forward_reducer( InType const & r ) { return r ; }

template< class OutType , class InType >
inline
typename std::enable_if< Kokkos::is_view< InType >::value , OutType >::type
forward_reducer( InType const & v )
{ return OutType( v.data() ); }

template< class OutType >
inline
OutType
forward_reducer( typename OutType::reference ref )
{ return OutType( ref ); }

} /* namespace Impl */

//----------------------------------------------------------------------------
// parallel_reduce with 4 args: label, policy, closure, and reducer

/**\brief  Parallel reduce with an explicit Reducer */
template< class PolicyType , class ClosureType , class ReduceType >
inline
typename std::enable_if< Kokkos::is_reducer< ReduceType >::value >::type
parallel_reduce( std::string const  & arg_label
               , PolicyType        && arg_policy
               , ClosureType       && arg_closure
               , ReduceType        && arg_reduce
               )
{
  //------------------------------

  using input_policy_type =
    typename std::remove_const<
      typename std::remove_reference< PolicyType >::type >::type ;

  using input_reduce_type =
    typename std::remove_const<
      typename std::remove_reference< ReduceType >::type >::type ;

  using Analysis = Kokkos::Impl::FunctorAnalysis
    < Kokkos::Impl::FunctorPatternInterface::REDUCE
    , input_policy_type
    , ClosureType
    > ;

  //------------------------------
  // Policy is either given or an integer value
  // If an integer value then is a RangePolicy with queried execution space

  enum { is_policy = Kokkos::is_execution_policy< input_policy_type >::value };
  enum { is_intval = std::is_integral< input_policy_type >::value };

  static_assert( is_policy || is_intval ,
    "Kokkos::parallel_reduce 2nd argument must be execution policy or integral value" );

  using policy_type = typename std::conditional
    < is_policy , input_policy_type
    , Kokkos::RangePolicy< typename Analysis::execution_space >
    >::type ;

  //------------------------------
  // ReduceType is either a reducer, view, or value reference

  enum { is_reducer = Kokkos::is_reducer< input_reduce_type >::value };
  enum { is_view    = Kokkos::is_view< input_reduce_type >::value };
  enum { is_ref     = std::is_same< ReduceType
                                  , typename Analysis::reference_type
                                  >::value };

  static_assert( is_reducer || is_view || is_ref ,
    "Kokkos::parallel_reduce 4th argument must be reducer, output View, or output variable" );

  // If input_reducer_type is_view or is_ref then need its memory_space.
  // A View has a memory_space, a reference is in the HostSpace.

  using has_space = typename std::conditional
    < is_view , input_reduce_type , Kokkos::HostSpace >::type ;

  using memory_space = typename has_space::memory_space ;

  using reduce_type = typename std::conditional
    < is_reducer , input_reduce_type
    , typename Analysis::Reducer< memory_space >
    >::type ;

  //------------------------------

  #if defined(KOKKOS_ENABLE_PROFILING)
  uint64_t kpID = 0;
  if(Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::beginParallelReduce(arg_label, 0, &kpID);
  }
  #endif

  //------------------------------
  // Disable tracking while creating the closure:

  Kokkos::Impl::shared_allocation_tracking_claim_and_disable();

  Kokkos::Impl::ParallelReduce< ClosureType , policy_type, reduce_type
                              , typename Analysis::execution_space >
    closure( arg_closure
           , forward_execution_policy< policy_type >( arg_policy )
           , forward_reducer< reduce_type >( arg_reduce ) );

  Kokkos::Impl::shared_allocation_tracking_release_and_enable();

  // Enable tracking after creating the closure

  closure.execute();

  //------------------------------

  #if defined(KOKKOS_ENABLE_PROFILING)
  if(Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::endParallelReduce(kpID);
  }
  #endif
}

//----------------------------------------------------------------------------
// parallel_reduce with 3 args: policy, closure, and reducer

template< class PolicyType , class ClosureType , class ReduceType >
inline
typename std::enable_if
  < Kokkos::is_execution_policy<
      typename std::remove_const<
      typename std::remove_reference< PolicyType >::type >::type
    >::value
    ||
    std::is_integral<
      typename std::remove_const<
      typename std::remove_reference< PolicyType >::type >::type
    >::value
  >::type ;
parallel_reduce( PolicyType   && arg_policy
               , ClosureType  && arg_closure
               , ReduceType   && arg_reduce
               )
{
  parallel_reduce( typeid(ClosureType).name()
                 , std::forward< PolicyType  >( arg_policy )
                 , std::forward< ClosureType >( arg_closure )
                 , std::forward< ReduceType  >( arg_reduce ) );
}

// parallel_reduce with 3 args: label, policy, and closure

template< class PolicyType , class ClosureType >
inline
void
parallel_reduce( std::string const & arg_label
               , PolicyType   && arg_policy
               , ClosureType  && arg_closure
               )
{
  // Deduce a Reducer from the Closure

  using input_policy_type =
    typename std::remove_const<
      typename std::remove_reference< PolicyType >::type >::type ;

  using Analysis = Kokkos::Impl::FunctorAnalysis
    < Kokkos::Impl::FunctorPatternInterface::REDUCE
    , input_policy_type
    , ClosureType
    > ;

  static_assert( Analysis::has_final_member_function ,
    "Kokkos::parallel_reduce functor does not have a final member function" );

  parallel_reduce( arg_label
                 , std::forward< PolicyType  >( arg_policy )
                 , std::forward< ClosureType >( arg_closure )
                 , typename Analysis::Reducer<>() );
}

//----------------------------------------------------------------------------
// parallel_reduce with 2 arguments: policy and closure:

/**\brief  Parallel reduce processed by ClosureType::final */
template< class PolicyType , class ClosureType >
inline
parallel_reduce( PolicyType  && arg_policy
               , ClosureType && arg_closure )
{
  // Deduce a Reducer from the Closure

  using input_policy_type =
    typename std::remove_const<
      typename std::remove_reference< PolicyType >::type >::type ;

  using Analysis = Kokkos::Impl::FunctorAnalysis
    < Kokkos::Impl::FunctorPatternInterface::REDUCE
    , input_policy_type
    , ClosureType
    > ;

  static_assert( Analysis::has_final_member_function ,
    "Kokkos::parallel_reduce functor does not have a final member function" );

  parallel_reduce( typeid(ClosureType).name()
                 , std::forward< PolicyType  >( arg_policy )
                 , std::forward< ClosureType >( arg_closure )
                 , typename Analysis::Reducer<>() );
}

#endif

//----------------------------------------------------------------------------

/*! \fn void parallel_reduce(label,policy,functor,return_argument)
    \brief Perform a parallel reduction.
    \param label An optional Label giving the call name. Must be able to construct a std::string from the argument.
    \param policy A Kokkos Execution Policy, such as an integer, a RangePolicy or a TeamPolicy.
    \param functor A functor with a reduction operator, and optional init, join and final functions.
    \param return_argument A return argument which can be a scalar, a View, or a ReducerStruct. This argument can be left out if the functor has a final function.
*/

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
 *    typedef    ...     execution_space ;
 *    typedef <podType>  value_type ;
 *    void operator()( <intType> iwork , <podType> & update ) const ;
 *    void init( <podType> & update ) const ;
 *    void join( volatile       <podType> & update ,
 *               volatile const <podType> & input ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> & update ) const ;
 *  };
 * \endcode
 *
 * Example of a parallel_reduce functor for an array of POD (plain old data) values:
 * \code
 *  class FunctorType { // For array of POD value
 *  public:
 *    typedef    ...     execution_space ;
 *    typedef <podType>  value_type[] ;
 *    void operator()( <intType> , <podType> update[] ) const ;
 *    void init( <podType> update[] ) const ;
 *    void join( volatile       <podType> update[] ,
 *               volatile const <podType> input[] ) const ;
 *
 *    typedef true_type has_final ;
 *    void final( <podType> update[] ) const ;
 *  };
 * \endcode
 */

// ReturnValue is scalar or array: take by reference

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const PolicyType& policy,
                     const FunctorType& functor,
                     ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,ReturnType>::execute(label,policy,functor,return_value);
}

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const PolicyType& policy,
                     const FunctorType& functor,
                     ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,ReturnType>::execute("",policy,functor,return_value);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const size_t& policy,
                     const FunctorType& functor,
                     ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute("",policy_type(0,policy),functor,return_value);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const size_t& policy,
                     const FunctorType& functor,
                     ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute(label,policy_type(0,policy),functor,return_value);
}

// ReturnValue as View or Reducer: take by copy to allow for inline construction

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const PolicyType& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,const ReturnType>::execute(label,policy,functor,return_value);
}

template< class PolicyType, class FunctorType, class ReturnType >
inline
void parallel_reduce(const PolicyType& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,ReturnType>::execute("",policy,functor,return_value_impl);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const size_t& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute("",policy_type(0,policy),functor,return_value_impl);
}

template< class FunctorType, class ReturnType >
inline
void parallel_reduce(const std::string& label,
                     const size_t& policy,
                     const FunctorType& functor,
                     const ReturnType& return_value) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<policy_type,FunctorType,ReturnType>::execute(label,policy_type(0,policy),functor,return_value_impl);
}

// No Return Argument

template< class PolicyType, class FunctorType>
inline
void parallel_reduce(const std::string& label,
                     const PolicyType& policy,
                     const FunctorType& functor,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,result_view_type>::execute(label,policy,functor,result_view);
}

template< class PolicyType, class FunctorType >
inline
void parallel_reduce(const PolicyType& policy,
                     const FunctorType& functor,
                     typename Impl::enable_if<
                       Kokkos::Impl::is_execution_policy<PolicyType>::value
                     >::type * = 0) {
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<PolicyType,FunctorType,result_view_type>::execute("",policy,functor,result_view);
}

template< class FunctorType >
inline
void parallel_reduce(const size_t& policy,
                     const FunctorType& functor) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<policy_type,FunctorType,result_view_type>::execute("",policy_type(0,policy),functor,result_view);
}

template< class FunctorType>
inline
void parallel_reduce(const std::string& label,
                     const size_t& policy,
                     const FunctorType& functor) {
  typedef typename Impl::ParallelReducePolicyType<void,size_t,FunctorType>::policy_type policy_type;
  typedef Kokkos::Impl::FunctorValueTraits< FunctorType , void >  ValueTraits ;
  typedef typename Kokkos::Impl::if_c< (ValueTraits::StaticValueSize != 0)
                                     , typename ValueTraits::value_type
                                     , typename ValueTraits::pointer_type
                                     >::type value_type ;

  typedef Kokkos::View< value_type
              , Kokkos::HostSpace
              , Kokkos::MemoryUnmanaged
              > result_view_type;
  result_view_type result_view ;

  Impl::ParallelReduceAdaptor<policy_type,FunctorType,result_view_type>::execute(label,policy_type(0,policy),functor,result_view);
}

} //namespace Kokkos

#endif // KOKKOS_PARALLEL_REDUCE_HPP

