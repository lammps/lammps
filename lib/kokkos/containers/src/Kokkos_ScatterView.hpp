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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/


/// \file Kokkos_ScatterView.hpp
/// \brief Declaration and definition of Kokkos::ScatterView.
///
/// This header file declares and defines Kokkos::ScatterView and its
/// related nonmember functions.

#ifndef KOKKOS_SCATTER_VIEW_HPP
#define KOKKOS_SCATTER_VIEW_HPP

#include <Kokkos_Core.hpp>
#include <utility>

namespace Kokkos {
namespace Experimental {

/*
 * Reduction Type list
 *  - These corresponds to subset of the reducers in parallel_reduce
 *  - See Implementations of ScatterValue for details.
 */
enum : int {
  ScatterSum,
  ScatterProd,
  ScatterMax,
  ScatterMin,
};

enum : int {
  ScatterNonDuplicated = 0,
  ScatterDuplicated    = 1
};

enum : int {
  ScatterNonAtomic = 0,
  ScatterAtomic    = 1
};

}} // Kokkos::Experimental

namespace Kokkos {
namespace Impl {
namespace Experimental {

template <typename ExecSpace>
struct DefaultDuplication;

template <typename ExecSpace, int duplication>
struct DefaultContribution;

#ifdef KOKKOS_ENABLE_SERIAL
template <>
struct DefaultDuplication<Kokkos::Serial> {
  enum : int { value = Kokkos::Experimental::ScatterNonDuplicated };
};
template <>
struct DefaultContribution<Kokkos::Serial, Kokkos::Experimental::ScatterNonDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterNonAtomic };
};
template <>
struct DefaultContribution<Kokkos::Serial, Kokkos::Experimental::ScatterDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterNonAtomic };
};
#endif

#ifdef KOKKOS_ENABLE_OPENMP
template <>
struct DefaultDuplication<Kokkos::OpenMP> {
  enum : int { value = Kokkos::Experimental::ScatterDuplicated };
};
template <>
struct DefaultContribution<Kokkos::OpenMP, Kokkos::Experimental::ScatterNonDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterAtomic };
};
template <>
struct DefaultContribution<Kokkos::OpenMP, Kokkos::Experimental::ScatterDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterNonAtomic };
};
#endif

#ifdef KOKKOS_ENABLE_HPX
template <>
struct DefaultDuplication<Kokkos::Experimental::HPX> {
  enum : int { value = Kokkos::Experimental::ScatterDuplicated };
};
template <>
struct DefaultContribution<Kokkos::Experimental::HPX, Kokkos::Experimental::ScatterNonDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterAtomic };
};
template <>
struct DefaultContribution<Kokkos::Experimental::HPX, Kokkos::Experimental::ScatterDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterNonAtomic };
};
#endif

#ifdef KOKKOS_ENABLE_THREADS
template <>
struct DefaultDuplication<Kokkos::Threads> {
  enum : int { value = Kokkos::Experimental::ScatterDuplicated };
};
template <>
struct DefaultContribution<Kokkos::Threads, Kokkos::Experimental::ScatterNonDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterAtomic };
};
template <>
struct DefaultContribution<Kokkos::Threads, Kokkos::Experimental::ScatterDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterNonAtomic };
};
#endif

#ifdef KOKKOS_ENABLE_CUDA
template <>
struct DefaultDuplication<Kokkos::Cuda> {
  enum : int { value = Kokkos::Experimental::ScatterNonDuplicated };
};
template <>
struct DefaultContribution<Kokkos::Cuda, Kokkos::Experimental::ScatterNonDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterAtomic };
};
template <>
struct DefaultContribution<Kokkos::Cuda, Kokkos::Experimental::ScatterDuplicated> {
  enum : int { value = Kokkos::Experimental::ScatterAtomic };
};
#endif

/* ScatterValue <Op=ScatterSum, contribution=ScatterNonAtomic> is the object returned by the access operator() of ScatterAccess,
   This class inherits from the Sum<> reducer and it wraps join(dest, src) with convenient operator+=, etc. 
   Note the addition of update(ValueType const& rhs) and reset()  so that all reducers can have common functions
   See ReduceDuplicates and ResetDuplicates ) */
template <typename ValueType, int Op, int contribution>
struct ScatterValue;

template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonAtomic> :
  Sum<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Sum<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ScatterValue&& other) : 
       Sum<ValueType,Kokkos::DefaultExecutionSpace>(other.reference())
    {}
    KOKKOS_FORCEINLINE_FUNCTION void operator+=(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void operator-=(ValueType const& rhs) {
      this->join( this->reference(), -rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }
};

/* ScatterValue <Op=ScatterSum, contribution=ScatterAtomic> is the object returned by the access operator() 
 * of ScatterAccess, similar to that returned by an Atomic View, it wraps Kokkos::atomic_add with convenient
   operator+=, etc. This version also has the update(rhs) and reset() functions. */
template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterAtomic> :
  Sum<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Sum<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}

    KOKKOS_FORCEINLINE_FUNCTION void operator+=(ValueType const& rhs) {
     this->join(this->reference(), rhs);
    }
    KOKKOS_FORCEINLINE_FUNCTION void operator-=(ValueType const& rhs) {
      this->join(this->reference(), -rhs);
    }
    
    KOKKOS_INLINE_FUNCTION
    void join(ValueType& dest, const ValueType& src)  const {
      Kokkos::atomic_add(&dest, src);
    }

    KOKKOS_INLINE_FUNCTION
    void join(volatile ValueType& dest, const volatile ValueType& src) const {
      Kokkos::atomic_add(&dest, src);
    } 

    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }

    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }
};

/* ScatterValue <Op=ScatterProd, contribution=ScatterNonAtomic> is the object returned by the access operator() of ScatterAccess,
   This class inherits from the Prod<> reducer and it wraps join(dest, src) with convenient operator*=, etc. 
   Note the addition of update(ValueType const& rhs) and reset()  so that all reducers can have common functions
   See ReduceDuplicates and ResetDuplicates ) */
template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterProd, Kokkos::Experimental::ScatterNonAtomic> :
  Prod<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Prod<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ScatterValue&& other) : 
       Prod<ValueType,Kokkos::DefaultExecutionSpace>(other.reference())
    {}
    KOKKOS_FORCEINLINE_FUNCTION void operator*=(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void operator/=(ValueType const& rhs) {
      this->join( this->reference(), static_cast<ValueType>(1)/rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }
};

/* ScatterValue <Op=ScatterProd, contribution=ScatterAtomic> is the object returned by the access operator() 
 * of ScatterAccess, similar to that returned by an Atomic View, it wraps and atomic_prod with convenient
   operator*=, etc. atomic_prod uses the atomic_compare_exchange. This version also has the update(rhs) and reset() functions. */
template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterProd, Kokkos::Experimental::ScatterAtomic> :
  Prod<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Prod<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}

    KOKKOS_FORCEINLINE_FUNCTION void operator*=(ValueType const& rhs) {
     this->join(this->reference(), rhs);
    }
    KOKKOS_FORCEINLINE_FUNCTION void operator/=(ValueType const& rhs) {
      this->join(this->reference(), static_cast<ValueType>(1)/rhs);
    }

    KOKKOS_FORCEINLINE_FUNCTION 
    void atomic_prod(ValueType & dest, const ValueType& src) const {

        bool success = false;
        while(!success) {
            ValueType dest_old = dest;
            ValueType dest_new = dest_old * src;
            dest_new = Kokkos::atomic_compare_exchange<ValueType>(&dest,dest_old,dest_new);
            success = ( (dest_new - dest_old)/dest_old <= 1e-15 );
        }
    }
    
    KOKKOS_INLINE_FUNCTION
    void join(ValueType& dest, const ValueType& src)  const {
      atomic_prod(dest, src);
    }

    KOKKOS_INLINE_FUNCTION
    void join(volatile ValueType& dest, const volatile ValueType& src) const {
      atomic_prod(dest, src);
    } 

    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }

};

/* ScatterValue <Op=ScatterMin, contribution=ScatterNonAtomic> is the object returned by the access operator() of ScatterAccess,
   This class inherits from the Min<> reducer and it wraps join(dest, src) with convenient update(rhs). 
   Note the addition of update(ValueType const& rhs) and reset() are so that all reducers can have a common update function
   See ReduceDuplicates and ResetDuplicates ) */
template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterMin, Kokkos::Experimental::ScatterNonAtomic> :
  Min<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Min<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ScatterValue&& other) : 
       Min<ValueType,Kokkos::DefaultExecutionSpace>(other.reference())
    {}
    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }
};

/* ScatterValue <Op=ScatterMin, contribution=ScatterAtomic> is the object returned by the access operator() 
 * of ScatterAccess, similar to that returned by an Atomic View, it wraps and atomic_min with the update(rhs)
   function. atomic_min uses the atomic_compare_exchange. This version also has the reset() function */
template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterMin, Kokkos::Experimental::ScatterAtomic> :
  Min<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Min<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}

    KOKKOS_FORCEINLINE_FUNCTION 
    void atomic_min(ValueType & dest, const ValueType& src) const {

        bool success = false;
        while(!success) {
            ValueType dest_old = dest;
            ValueType dest_new = ( dest_old > src ) ? src : dest_old;
            dest_new = Kokkos::atomic_compare_exchange<ValueType>(&dest,dest_old,dest_new);
            success = ( (dest_new - dest_old)/dest_old <= 1e-15 );
        }
    }
    
    KOKKOS_INLINE_FUNCTION
    void join(ValueType& dest, const ValueType& src)  const {
      atomic_min(dest, src);
    }

    KOKKOS_INLINE_FUNCTION
    void join(volatile ValueType& dest, const volatile ValueType& src) const {
      atomic_min(dest, src);
    } 

    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }

};

/* ScatterValue <Op=ScatterMax, contribution=ScatterNonAtomic> is the object returned by the access operataor() of ScatterAccess,
   This class inherits from the Max<> reducer and it wraps join(dest, src) with convenient update(rhs). 
   Note the addition of update(ValueType const& rhs) and reset() are so that all reducers can have a common update function
   See ReduceDuplicates and ResetDuplicates ) */
template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterMax, Kokkos::Experimental::ScatterNonAtomic> :
  Max<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Max<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ScatterValue&& other) : 
       Max<ValueType,Kokkos::DefaultExecutionSpace>(other.reference())
    {}
    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }
};

/* ScatterValue <Op=ScatterMax, contribution=ScatterAtomic> is the object returned by the access operator() 
 * of ScatterAccess, similar to that returned by an Atomic View, it wraps and atomic_max with the update(rhs)
   function. atomic_max uses the atomic_compare_exchange. This version also has the reset() function  */
template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterMax, Kokkos::Experimental::ScatterAtomic> :
  Max<ValueType,Kokkos::DefaultExecutionSpace> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : 
       Max<ValueType,Kokkos::DefaultExecutionSpace>(value_in)
    {}

    KOKKOS_FORCEINLINE_FUNCTION 
    void atomic_max(ValueType & dest, const ValueType& src) const {

        bool success = false;
        while(!success) {
            ValueType dest_old = dest;
            ValueType dest_new = ( dest_old < src ) ? src : dest_old;
            dest_new = Kokkos::atomic_compare_exchange<ValueType>(&dest,dest_old,dest_new);
            success = ( (dest_new - dest_old)/dest_old <= 1e-15 );
        }
    }
    
    KOKKOS_INLINE_FUNCTION
    void join(ValueType& dest, const ValueType& src)  const {
      atomic_max(dest, src);
    }

    KOKKOS_INLINE_FUNCTION
    void join(volatile ValueType& dest, const volatile ValueType& src) const {
      atomic_max(dest, src);
    } 

    KOKKOS_FORCEINLINE_FUNCTION void update(ValueType const& rhs) {
      this->join( this->reference(), rhs );
    }
    KOKKOS_FORCEINLINE_FUNCTION void reset() {
      this->init( this->reference() );
    }

};

/* DuplicatedDataType, given a View DataType, will create a new DataType
   that has a new runtime dimension which becomes the largest-stride dimension.
   In the case of LayoutLeft, due to the limitation induced by the design of DataType
   itself, it must convert any existing compile-time dimensions into runtime dimensions. */
template <typename T, typename Layout>
struct DuplicatedDataType;

template <typename T>
struct DuplicatedDataType<T, Kokkos::LayoutRight> {
  typedef T* value_type; // For LayoutRight, add a star all the way on the left
};

template <typename T, size_t N>
struct DuplicatedDataType<T[N], Kokkos::LayoutRight> {
  typedef typename DuplicatedDataType<T, Kokkos::LayoutRight>::value_type value_type[N];
};

template <typename T>
struct DuplicatedDataType<T[], Kokkos::LayoutRight> {
  typedef typename DuplicatedDataType<T, Kokkos::LayoutRight>::value_type value_type[];
};

template <typename T>
struct DuplicatedDataType<T*, Kokkos::LayoutRight> {
  typedef typename DuplicatedDataType<T, Kokkos::LayoutRight>::value_type* value_type;
};

template <typename T>
struct DuplicatedDataType<T, Kokkos::LayoutLeft> {
  typedef T* value_type;
};

template <typename T, size_t N>
struct DuplicatedDataType<T[N], Kokkos::LayoutLeft> {
  typedef typename DuplicatedDataType<T, Kokkos::LayoutLeft>::value_type* value_type;
};

template <typename T>
struct DuplicatedDataType<T[], Kokkos::LayoutLeft> {
  typedef typename DuplicatedDataType<T, Kokkos::LayoutLeft>::value_type* value_type;
};

template <typename T>
struct DuplicatedDataType<T*, Kokkos::LayoutLeft> {
  typedef typename DuplicatedDataType<T, Kokkos::LayoutLeft>::value_type* value_type;
};

/* Insert integer argument pack into array */

template<class T>
void args_to_array(size_t* array, int pos, T dim0) {
  array[pos] = dim0;
}
template<class T, class ... Dims>
void args_to_array(size_t* array, int pos, T dim0, Dims ... dims) {
  array[pos] = dim0;
  args_to_array(array,pos+1,dims...);
}

/* Slice is just responsible for stuffing the correct number of Kokkos::ALL
   arguments on the correct side of the index in a call to subview() to get a
   subview where the index specified is the largest-stride one. */
template <typename Layout, int rank, typename V, typename ... Args>
struct Slice {
  typedef Slice<Layout, rank - 1, V, Kokkos::Impl::ALL_t, Args...> next;
  typedef typename next::value_type value_type;

  static
  value_type get(V const& src, const size_t i, Args ... args) {
    return next::get(src, i, Kokkos::ALL, args...);
  }
};

template <typename V, typename ... Args>
struct Slice<Kokkos::LayoutRight, 1, V, Args...> {
  typedef typename Kokkos::Impl::ViewMapping
                          < void
                          , V
                          , const size_t
                          , Args ...
                          >::type value_type;
  static
  value_type get(V const& src, const size_t i, Args ... args) {
    return Kokkos::subview(src, i, args...);
  }
};

template <typename V, typename ... Args>
struct Slice<Kokkos::LayoutLeft, 1, V, Args...> {
  typedef typename Kokkos::Impl::ViewMapping
                          < void
                          , V
                          , Args ...
                          , const size_t
                          >::type value_type;
  static
  value_type get(V const& src, const size_t i, Args ... args) {
    return Kokkos::subview(src, args..., i);
  }
};

template <typename ExecSpace, typename ValueType, int Op>
struct ReduceDuplicates;

template <typename ExecSpace, typename ValueType, int Op>
struct ReduceDuplicatesBase {
  typedef ReduceDuplicates<ExecSpace, ValueType, Op> Derived;
  ValueType const* src;
  ValueType* dst;
  size_t stride;
  size_t start;
  size_t n;
  ReduceDuplicatesBase(ValueType const* src_in, ValueType* dest_in, size_t stride_in, size_t start_in, size_t n_in, std::string const& name)
    : src(src_in)
    , dst(dest_in)
    , stride(stride_in)
    , start(start_in)
    , n(n_in)
  {
#if defined(KOKKOS_ENABLE_PROFILING)
    uint64_t kpID = 0;
    if(Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::beginParallelFor(std::string("reduce_") + name, 0, &kpID);
    }
#endif
    typedef RangePolicy<ExecSpace, size_t> policy_type;
    typedef Kokkos::Impl::ParallelFor<Derived, policy_type> closure_type;
    const closure_type closure(*(static_cast<Derived*>(this)), policy_type(0, stride));
    closure.execute();
#if defined(KOKKOS_ENABLE_PROFILING)
    if(Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endParallelFor(kpID);
    }
#endif
  }
};

/* ReduceDuplicates -- Perform reduction on destination array using strided source 
 *    Use ScatterValue<> specific to operation to wrap destination array so that
 *    the reduction operation can be accessed via the update(rhs) function */
template <typename ExecSpace, typename ValueType, int Op>
struct ReduceDuplicates :
  public ReduceDuplicatesBase<ExecSpace, ValueType, Op>
{
  typedef ReduceDuplicatesBase<ExecSpace, ValueType, Op> Base;
  ReduceDuplicates(ValueType const* src_in, ValueType* dst_in, size_t stride_in, size_t start_in, size_t n_in, std::string const& name):
    Base(src_in, dst_in, stride_in, start_in, n_in, name)
  {}
  KOKKOS_FORCEINLINE_FUNCTION void operator()(size_t i) const {
    for (size_t j = Base::start; j < Base::n; ++j) {
      ScatterValue<ValueType, Op, Kokkos::Experimental::ScatterNonAtomic> sv(Base::dst[i]);
      sv.update( Base::src[i + Base::stride * j] );
    }
  }
};


template <typename ExecSpace, typename ValueType, int Op>
struct ResetDuplicates;

template <typename ExecSpace, typename ValueType, int Op>
struct ResetDuplicatesBase {
  typedef ResetDuplicates<ExecSpace, ValueType, Op> Derived;
  ValueType* data;
  ResetDuplicatesBase(ValueType* data_in, size_t size_in, std::string const& name)
    : data(data_in)
  {
#if defined(KOKKOS_ENABLE_PROFILING)
    uint64_t kpID = 0;
    if(Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::beginParallelFor(std::string("reduce_") + name, 0, &kpID);
    }
#endif
    typedef RangePolicy<ExecSpace, size_t> policy_type;
    typedef Kokkos::Impl::ParallelFor<Derived, policy_type> closure_type;
    const closure_type closure(*(static_cast<Derived*>(this)), policy_type(0, size_in));
    closure.execute();
#if defined(KOKKOS_ENABLE_PROFILING)
    if(Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::endParallelFor(kpID);
    }
#endif
  }
};

/* ResetDuplicates -- Perform reset on destination array
 *    Use ScatterValue<> specific to operation to wrap destination array so that
 *    the reset operation can be accessed via the reset() function */
template <typename ExecSpace, typename ValueType, int Op>
struct ResetDuplicates :
  public ResetDuplicatesBase<ExecSpace, ValueType, Op>
{
  typedef ResetDuplicatesBase<ExecSpace, ValueType, Op> Base;
  ResetDuplicates(ValueType* data_in, size_t size_in, std::string const& name):
    Base(data_in, size_in, name)
  {}
  KOKKOS_FORCEINLINE_FUNCTION void operator()(size_t i) const {
    ScatterValue<ValueType, Op, Kokkos::Experimental::ScatterNonAtomic> sv(Base::data[i]);
    sv.reset();
  }
};


}}} // Kokkos::Impl::Experimental

namespace Kokkos {
namespace Experimental {

template <typename DataType
         ,typename Layout = Kokkos::DefaultExecutionSpace::array_layout
         ,typename ExecSpace = Kokkos::DefaultExecutionSpace
         ,int Op = ScatterSum
         ,int duplication = Kokkos::Impl::Experimental::DefaultDuplication<ExecSpace>::value
         ,int contribution = Kokkos::Impl::Experimental::DefaultContribution<ExecSpace, duplication>::value
         >
class ScatterView;

template <typename DataType
         ,int Op
         ,typename ExecSpace
         ,typename Layout
         ,int duplication
         ,int contribution
         ,int override_contribution
         >
class ScatterAccess;

// non-duplicated implementation
template <typename DataType
         ,int Op
         ,typename ExecSpace
         ,typename Layout
         ,int contribution
         >
class ScatterView<DataType
                   ,Layout
                   ,ExecSpace
                   ,Op
                   ,ScatterNonDuplicated
                   ,contribution>
{
public:
  typedef Kokkos::View<DataType, Layout, ExecSpace> original_view_type;
  typedef typename original_view_type::value_type original_value_type;
  typedef typename original_view_type::reference_type original_reference_type;
  friend class ScatterAccess<DataType, Op, ExecSpace, Layout, ScatterNonDuplicated, contribution, ScatterNonAtomic>;
  friend class ScatterAccess<DataType, Op, ExecSpace, Layout, ScatterNonDuplicated, contribution, ScatterAtomic>;

  ScatterView()
  {
  }

  template <typename RT, typename ... RP>
  ScatterView(View<RT, RP...> const& original_view)
  : internal_view(original_view)
  {
  }

  template <typename ... Dims>
  ScatterView(std::string const& name, Dims ... dims)
  : internal_view(name, dims ...)
  {
  }

  template <int override_contrib = contribution>
  KOKKOS_FORCEINLINE_FUNCTION
  ScatterAccess<DataType, Op, ExecSpace, Layout, ScatterNonDuplicated, contribution, override_contrib>
  access() const {
    return ScatterAccess<DataType, Op, ExecSpace, Layout, ScatterNonDuplicated, contribution, override_contrib>{*this};
  }

  original_view_type subview() const {
    return internal_view;
  }

  template <typename DT, typename ... RP>
  void contribute_into(View<DT, RP...> const& dest) const
  {
    typedef View<DT, RP...> dest_type;
    static_assert(std::is_same<
        typename dest_type::array_layout,
        Layout>::value,
        "ScatterView contribute destination has different layout");
    static_assert(Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
        typename ExecSpace::memory_space,
        typename dest_type::memory_space>::value,
        "ScatterView contribute destination memory space not accessible");
    if (dest.data() == internal_view.data()) return;
    Kokkos::Impl::Experimental::ReduceDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data(),
        dest.data(),
        0,
        0,
        1,
        internal_view.label());
  }

  void reset() {
    Kokkos::Impl::Experimental::ResetDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data(),
        internal_view.size(),
        internal_view.label());
  }
  template <typename DT, typename ... RP>
  void reset_except(View<DT, RP...> const& view) {
    if (view.data() != internal_view.data()) reset();
  }

  void resize(const size_t n0 = 0,
           const size_t n1 = 0,
           const size_t n2 = 0,
           const size_t n3 = 0,
           const size_t n4 = 0,
           const size_t n5 = 0,
           const size_t n6 = 0,
           const size_t n7 = 0) {
    ::Kokkos::resize(internal_view,n0,n1,n2,n3,n4,n5,n6,n7);
  }

  void realloc(const size_t n0 = 0,
           const size_t n1 = 0,
           const size_t n2 = 0,
           const size_t n3 = 0,
           const size_t n4 = 0,
           const size_t n5 = 0,
           const size_t n6 = 0,
           const size_t n7 = 0) {
    ::Kokkos::realloc(internal_view,n0,n1,n2,n3,n4,n5,n6,n7);
  }

protected:
  template <typename ... Args>
  KOKKOS_FORCEINLINE_FUNCTION
  original_reference_type at(Args ... args) const {
    return internal_view(args...);
  }
private:
  typedef original_view_type internal_view_type;
  internal_view_type internal_view;
};

template <typename DataType
         ,int Op
         ,typename ExecSpace
         ,typename Layout
         ,int contribution
         ,int override_contribution
         >
class ScatterAccess<DataType
                   ,Op
                   ,ExecSpace
                   ,Layout
                   ,ScatterNonDuplicated
                   ,contribution
                   ,override_contribution>
{
public:
  typedef ScatterView<DataType, Layout, ExecSpace, Op, ScatterNonDuplicated, contribution> view_type;
  typedef typename view_type::original_value_type original_value_type;
  typedef Kokkos::Impl::Experimental::ScatterValue<
      original_value_type, Op, override_contribution> value_type;

  KOKKOS_INLINE_FUNCTION
  ScatterAccess() :
    view(view_type())  {
  }

  KOKKOS_INLINE_FUNCTION
  ScatterAccess(view_type const& view_in)
    : view(view_in)
  {
  }

  KOKKOS_INLINE_FUNCTION
  ~ScatterAccess()
  {
  }

  template <typename ... Args>
  KOKKOS_FORCEINLINE_FUNCTION
  value_type operator()(Args ... args) const {
    return view.at(args...);
  }

  template <typename Arg>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<view_type::original_view_type::rank == 1 &&
  std::is_integral<Arg>::value, value_type>::type
  operator[](Arg arg) const {
    return view.at(arg);
  }

private:
  view_type const& view;
};

// duplicated implementation
// LayoutLeft and LayoutRight are different enough that we'll just specialize each

template <typename DataType
         ,int Op
         ,typename ExecSpace
         ,int contribution
         >
class ScatterView<DataType
                   ,Kokkos::LayoutRight
                   ,ExecSpace
                   ,Op
                   ,ScatterDuplicated
                   ,contribution>
{
public:
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, ExecSpace> original_view_type;
  typedef typename original_view_type::value_type original_value_type;
  typedef typename original_view_type::reference_type original_reference_type;
  friend class ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutRight, ScatterDuplicated, contribution, ScatterNonAtomic>;
  friend class ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutRight, ScatterDuplicated, contribution, ScatterAtomic>;
  typedef typename Kokkos::Impl::Experimental::DuplicatedDataType<DataType, Kokkos::LayoutRight> data_type_info;
  typedef typename data_type_info::value_type internal_data_type;
  typedef Kokkos::View<internal_data_type, Kokkos::LayoutRight, ExecSpace> internal_view_type;

  ScatterView()
  {
  }

  template <typename RT, typename ... RP >
  ScatterView(View<RT, RP...> const& original_view)
  : unique_token()
  , internal_view(Kokkos::ViewAllocateWithoutInitializing(
                    std::string("duplicated_") + original_view.label()),
                  unique_token.size(),
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
                  original_view.extent(0),
                  original_view.extent(1),
                  original_view.extent(2),
                  original_view.extent(3),
                  original_view.extent(4),
                  original_view.extent(5),
                  original_view.extent(6) )
#else
                  original_view.rank_dynamic > 0 ? original_view.extent(0): KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                  original_view.rank_dynamic > 1 ? original_view.extent(1): KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                  original_view.rank_dynamic > 2 ? original_view.extent(2): KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                  original_view.rank_dynamic > 3 ? original_view.extent(3): KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                  original_view.rank_dynamic > 4 ? original_view.extent(4): KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                  original_view.rank_dynamic > 5 ? original_view.extent(5): KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                  original_view.rank_dynamic > 6 ? original_view.extent(6): KOKKOS_IMPL_CTOR_DEFAULT_ARG)

#endif
  {
    reset();
  }

  template <typename ... Dims>
  ScatterView(std::string const& name, Dims ... dims)
  : internal_view(Kokkos::ViewAllocateWithoutInitializing(name), unique_token.size(), dims ...)
  {
    reset();
  }

  template <int override_contribution = contribution>
  KOKKOS_FORCEINLINE_FUNCTION
  ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutRight, ScatterDuplicated, contribution, override_contribution>
  access() const {
    return ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutRight, ScatterDuplicated, contribution, override_contribution>{*this};
  }

  typename Kokkos::Impl::Experimental::Slice<
    Kokkos::LayoutRight, internal_view_type::rank, internal_view_type>::value_type
  subview() const
  {
    return Kokkos::Impl::Experimental::Slice<
      Kokkos::LayoutRight, internal_view_type::Rank, internal_view_type>::get(internal_view, 0);
  }

  template <typename DT, typename ... RP>
  void contribute_into(View<DT, RP...> const& dest) const
  {
    typedef View<DT, RP...> dest_type;
    static_assert(std::is_same<
        typename dest_type::array_layout,
        Kokkos::LayoutRight>::value,
        "ScatterView deep_copy destination has different layout");
    static_assert(Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
        typename ExecSpace::memory_space,
        typename dest_type::memory_space>::value,
        "ScatterView deep_copy destination memory space not accessible");
    bool is_equal = (dest.data() == internal_view.data());
    size_t start = is_equal ? 1 : 0;
    Kokkos::Impl::Experimental::ReduceDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data(),
        dest.data(),
        internal_view.stride(0),
        start,
        internal_view.extent(0),
        internal_view.label());
  }

  void reset() {
    Kokkos::Impl::Experimental::ResetDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data(),
        internal_view.size(),
        internal_view.label());
  }
  template <typename DT, typename ... RP>
  void reset_except(View<DT, RP...> const& view) {
    if (view.data() != internal_view.data()) {
      reset();
      return;
    }
    Kokkos::Impl::Experimental::ResetDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data() + view.size(),
        internal_view.size() - view.size(),
        internal_view.label());
  }

  void resize(const size_t n0 = 0,
           const size_t n1 = 0,
           const size_t n2 = 0,
           const size_t n3 = 0,
           const size_t n4 = 0,
           const size_t n5 = 0,
           const size_t n6 = 0) {
    ::Kokkos::resize(internal_view,unique_token.size(),n0,n1,n2,n3,n4,n5,n6);
  }

  void realloc(const size_t n0 = 0,
           const size_t n1 = 0,
           const size_t n2 = 0,
           const size_t n3 = 0,
           const size_t n4 = 0,
           const size_t n5 = 0,
           const size_t n6 = 0) {
    ::Kokkos::realloc(internal_view,unique_token.size(),n0,n1,n2,n3,n4,n5,n6);
  }

protected:
  template <typename ... Args>
  KOKKOS_FORCEINLINE_FUNCTION
  original_reference_type at(int rank, Args ... args) const {
    return internal_view(rank, args...);
  }

protected:
  typedef Kokkos::Experimental::UniqueToken<
      ExecSpace, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;

  unique_token_type unique_token;
  internal_view_type internal_view;
};

template <typename DataType
         ,int Op
         ,typename ExecSpace
         ,int contribution
         >
class ScatterView<DataType
                   ,Kokkos::LayoutLeft
                   ,ExecSpace
                   ,Op
                   ,ScatterDuplicated
                   ,contribution>
{
public:
  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, ExecSpace> original_view_type;
  typedef typename original_view_type::value_type original_value_type;
  typedef typename original_view_type::reference_type original_reference_type;
  friend class ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutLeft, ScatterDuplicated, contribution, ScatterNonAtomic>;
  friend class ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutLeft, ScatterDuplicated, contribution, ScatterAtomic>;
  typedef typename Kokkos::Impl::Experimental::DuplicatedDataType<DataType, Kokkos::LayoutLeft> data_type_info;
  typedef typename data_type_info::value_type internal_data_type;
  typedef Kokkos::View<internal_data_type, Kokkos::LayoutLeft, ExecSpace> internal_view_type;

  ScatterView()
  {
  }

  template <typename RT, typename ... RP >
  ScatterView(View<RT, RP...> const& original_view)
  : unique_token()
  {
    size_t arg_N[8] = {
      original_view.rank>0?original_view.extent(0):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>1?original_view.extent(1):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>2?original_view.extent(2):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>3?original_view.extent(3):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>4?original_view.extent(4):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>5?original_view.extent(5):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>6?original_view.extent(6):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      KOKKOS_IMPL_CTOR_DEFAULT_ARG
    };
    arg_N[internal_view_type::rank - 1] = unique_token.size();
    internal_view = internal_view_type(
        Kokkos::ViewAllocateWithoutInitializing(
          std::string("duplicated_") + original_view.label()),
        arg_N[0], arg_N[1], arg_N[2], arg_N[3],
        arg_N[4], arg_N[5], arg_N[6], arg_N[7]);
    reset();
  }

  template <typename ... Dims>
  ScatterView(std::string const& name, Dims ... dims) {
    original_view_type original_view;
    size_t arg_N[8] = {
      original_view.rank>0?original_view.static_extent(0):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>1?original_view.static_extent(1):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>2?original_view.static_extent(2):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>3?original_view.static_extent(3):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>4?original_view.static_extent(4):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>5?original_view.static_extent(5):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      original_view.rank>6?original_view.static_extent(6):KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      KOKKOS_IMPL_CTOR_DEFAULT_ARG
    };
    Kokkos::Impl::Experimental::args_to_array(arg_N,0,dims ...);
    arg_N[internal_view_type::rank - 1] = unique_token.size();
    internal_view = internal_view_type(Kokkos::ViewAllocateWithoutInitializing(name),
     arg_N[0], arg_N[1], arg_N[2], arg_N[3],
     arg_N[4], arg_N[5], arg_N[6], arg_N[7]);
    reset();
  }

  template <int override_contribution = contribution>
  KOKKOS_FORCEINLINE_FUNCTION
  ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutLeft, ScatterDuplicated, contribution, override_contribution>
  access() const {
    return ScatterAccess<DataType, Op, ExecSpace, Kokkos::LayoutLeft, ScatterDuplicated, contribution, override_contribution>{*this};
  }

  typename Kokkos::Impl::Experimental::Slice<
    Kokkos::LayoutLeft, internal_view_type::rank, internal_view_type>::value_type
  subview() const
  {
    return Kokkos::Impl::Experimental::Slice<
      Kokkos::LayoutLeft, internal_view_type::rank, internal_view_type>::get(internal_view, 0);
  }

  template <typename ... RP>
  void contribute_into(View<RP...> const& dest) const
  {
    typedef View<RP...> dest_type;
    static_assert(std::is_same<
        typename dest_type::value_type,
        typename original_view_type::non_const_value_type>::value,
        "ScatterView deep_copy destination has wrong value_type");
    static_assert(std::is_same<
        typename dest_type::array_layout,
        Kokkos::LayoutLeft>::value,
        "ScatterView deep_copy destination has different layout");
    static_assert(Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
        typename ExecSpace::memory_space,
        typename dest_type::memory_space>::value,
        "ScatterView deep_copy destination memory space not accessible");
    auto extent = internal_view.extent(
        internal_view_type::rank - 1);
    bool is_equal = (dest.data() == internal_view.data());
    size_t start = is_equal ? 1 : 0;
    Kokkos::Impl::Experimental::ReduceDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data(),
        dest.data(),
        internal_view.stride(internal_view_type::rank - 1),
        start,
        extent,
        internal_view.label());
  }

  void reset() {
    Kokkos::Impl::Experimental::ResetDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data(),
        internal_view.size(),
        internal_view.label());
  }
  template <typename DT, typename ... RP>
  void reset_except(View<DT, RP...> const& view) {
    if (view.data() != internal_view.data()) {
      reset();
      return;
    }
    Kokkos::Impl::Experimental::ResetDuplicates<ExecSpace, original_value_type, Op>(
        internal_view.data() + view.size(),
        internal_view.size() - view.size(),
        internal_view.label());
  }

  void resize(const size_t n0 = 0,
           const size_t n1 = 0,
           const size_t n2 = 0,
           const size_t n3 = 0,
           const size_t n4 = 0,
           const size_t n5 = 0,
           const size_t n6 = 0) {

    size_t arg_N[8] = {n0,n1,n2,n3,n4,n5,n6,0};
    const int i = internal_view.rank-1;
    arg_N[i] = unique_token.size();

    ::Kokkos::resize(internal_view,
        arg_N[0], arg_N[1], arg_N[2], arg_N[3],
        arg_N[4], arg_N[5], arg_N[6], arg_N[7]);
  }

  void realloc(const size_t n0 = 0,
           const size_t n1 = 0,
           const size_t n2 = 0,
           const size_t n3 = 0,
           const size_t n4 = 0,
           const size_t n5 = 0,
           const size_t n6 = 0) {

    size_t arg_N[8] = {n0,n1,n2,n3,n4,n5,n6,0};
    const int i = internal_view.rank-1;
    arg_N[i] = unique_token.size();

    ::Kokkos::realloc(internal_view,
        arg_N[0], arg_N[1], arg_N[2], arg_N[3],
        arg_N[4], arg_N[5], arg_N[6], arg_N[7]);
  }

protected:
  template <typename ... Args>
  inline original_reference_type at(int thread_id, Args ... args) const {
    return internal_view(args..., thread_id);
  }

protected:
  typedef Kokkos::Experimental::UniqueToken<
      ExecSpace, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;

  unique_token_type unique_token;
  internal_view_type internal_view;
};


/* This object has to be separate in order to store the thread ID, which cannot
   be obtained until one is inside a parallel construct, and may be relatively
   expensive to obtain at every contribution
   (calls a non-inlined function, looks up a thread-local variable).
   Due to the expense, it is sensible to query it at most once per parallel iterate
   (ideally once per thread, but parallel_for doesn't expose that)
   and then store it in a stack variable.
   ScatterAccess serves as a non-const object on the stack which can store the thread ID */

template <typename DataType
         ,int Op
         ,typename ExecSpace
         ,typename Layout
         ,int contribution
         ,int override_contribution
         >
class ScatterAccess<DataType
                   ,Op
                   ,ExecSpace
                   ,Layout
                   ,ScatterDuplicated
                   ,contribution
                   ,override_contribution>
{
public:
  typedef ScatterView<DataType, Layout, ExecSpace, Op, ScatterDuplicated, contribution> view_type;
  typedef typename view_type::original_value_type original_value_type;
  typedef Kokkos::Impl::Experimental::ScatterValue<
      original_value_type, Op, override_contribution> value_type;

  KOKKOS_FORCEINLINE_FUNCTION
  ScatterAccess(view_type const& view_in)
    : view(view_in)
    , thread_id(view_in.unique_token.acquire()) {
  }

  KOKKOS_FORCEINLINE_FUNCTION
  ~ScatterAccess() {
    if (thread_id != ~thread_id_type(0)) view.unique_token.release(thread_id);
  }

  template <typename ... Args>
  KOKKOS_FORCEINLINE_FUNCTION
  value_type operator()(Args ... args) const {
    return view.at(thread_id, args...);
  }

  template <typename Arg>
  KOKKOS_FORCEINLINE_FUNCTION
  typename std::enable_if<view_type::original_view_type::rank == 1 &&
  std::is_integral<Arg>::value, value_type>::type
  operator[](Arg arg) const {
    return view.at(thread_id, arg);
  }

private:

  view_type const& view;

  // simplify RAII by disallowing copies
  ScatterAccess(ScatterAccess const& other) = delete;
  ScatterAccess& operator=(ScatterAccess const& other) = delete;
  ScatterAccess& operator=(ScatterAccess&& other) = delete;

public:
  // do need to allow moves though, for the common
  // auto b = a.access();
  // that assignments turns into a move constructor call  
  KOKKOS_FORCEINLINE_FUNCTION
  ScatterAccess(ScatterAccess&& other)
    : view(other.view)
    , thread_id(other.thread_id)
  {
    other.thread_id = ~thread_id_type(0);
  }

private:

  typedef typename view_type::unique_token_type unique_token_type;
  typedef typename unique_token_type::size_type thread_id_type;
  thread_id_type thread_id;
};

template <int Op = Kokkos::Experimental::ScatterSum,
          int duplication = -1,
          int contribution = -1,
          typename RT, typename ... RP>
ScatterView
  < RT
  , typename ViewTraits<RT, RP...>::array_layout
  , typename ViewTraits<RT, RP...>::execution_space
  , Op
  /* just setting defaults if not specified... things got messy because the view type
     does not come before the duplication/contribution settings in the
     template parameter list */
  , duplication == -1 ? Kokkos::Impl::Experimental::DefaultDuplication<typename ViewTraits<RT, RP...>::execution_space>::value : duplication
  , contribution == -1 ?
      Kokkos::Impl::Experimental::DefaultContribution<
                        typename ViewTraits<RT, RP...>::execution_space,
                        (duplication == -1 ?
                           Kokkos::Impl::Experimental::DefaultDuplication<
                             typename ViewTraits<RT, RP...>::execution_space
                             >::value
                                           : duplication
                        )
                        >::value
                       : contribution
  >
create_scatter_view(View<RT, RP...> const& original_view) {
  return original_view; // implicit ScatterView constructor call
}

}} // namespace Kokkos::Experimental

namespace Kokkos {
namespace Experimental {

template <typename DT1, typename DT2, typename LY, typename ES,  int OP, int CT, int DP, typename ... VP>
void
contribute(View<DT1, VP...>& dest, Kokkos::Experimental::ScatterView<DT2, LY, ES, OP, CT, DP> const& src)
{
  src.contribute_into(dest);
}

}} // namespace Kokkos::Experimental

namespace Kokkos {

template <typename DT, typename LY, typename ES,  int OP, int CT, int DP, typename ... IS>
void
realloc(Kokkos::Experimental::ScatterView<DT, LY, ES, OP, CT, DP>& scatter_view, IS ... is)
{
  scatter_view.realloc(is ...);
}

template <typename DT, typename LY, typename ES,  int OP, int CT, int DP, typename ... IS>
void
resize(Kokkos::Experimental::ScatterView<DT, LY, ES, OP, CT, DP>& scatter_view, IS ... is)
{
  scatter_view.resize(is ...);
}

} // namespace Kokkos

#endif
