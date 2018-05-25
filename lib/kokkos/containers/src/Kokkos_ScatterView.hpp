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

//TODO: replace this enum with the Kokkos::Sum, etc reducers for parallel_reduce
enum : int {
  ScatterSum,
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

/* ScatterValue is the object returned by the access operator() of ScatterAccess,
   similar to that returned by an Atomic View, it wraps Kokkos::atomic_add with convenient
   operator+=, etc. */
template <typename ValueType, int Op, int contribution>
struct ScatterValue;

template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonAtomic> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : value( value_in ) {}
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ScatterValue&& other) : value( other.value ) {}
    KOKKOS_FORCEINLINE_FUNCTION void operator+=(ValueType const& rhs) {
      value += rhs;
    }
    KOKKOS_FORCEINLINE_FUNCTION void operator-=(ValueType const& rhs) {
      value -= rhs;
    }
  private:
    ValueType& value;
};

template <typename ValueType>
struct ScatterValue<ValueType, Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterAtomic> {
  public:
    KOKKOS_FORCEINLINE_FUNCTION ScatterValue(ValueType& value_in) : value( value_in ) {}
    KOKKOS_FORCEINLINE_FUNCTION void operator+=(ValueType const& rhs) {
      Kokkos::atomic_add(&value, rhs);
    }
    KOKKOS_FORCEINLINE_FUNCTION void operator-=(ValueType const& rhs) {
      Kokkos::atomic_add(&value, -rhs);
    }
  private:
    ValueType& value;
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

template <typename ExecSpace, typename ValueType>
struct ReduceDuplicates<ExecSpace, ValueType, Kokkos::Experimental::ScatterSum> :
  public ReduceDuplicatesBase<ExecSpace, ValueType, Kokkos::Experimental::ScatterSum>
{
  typedef ReduceDuplicatesBase<ExecSpace, ValueType, Kokkos::Experimental::ScatterSum> Base;
  ReduceDuplicates(ValueType const* src_in, ValueType* dst_in, size_t stride_in, size_t start_in, size_t n_in, std::string const& name):
    Base(src_in, dst_in, stride_in, start_in, n_in, name)
  {}
  KOKKOS_FORCEINLINE_FUNCTION void operator()(size_t i) const {
    for (size_t j = Base::start; j < Base::n; ++j) {
      Base::dst[i] += Base::src[i + Base::stride * j];
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

template <typename ExecSpace, typename ValueType>
struct ResetDuplicates<ExecSpace, ValueType, Kokkos::Experimental::ScatterSum> :
  public ResetDuplicatesBase<ExecSpace, ValueType, Kokkos::Experimental::ScatterSum>
{
  typedef ResetDuplicatesBase<ExecSpace, ValueType, Kokkos::Experimental::ScatterSum> Base;
  ResetDuplicates(ValueType* data_in, size_t size_in, std::string const& name):
    Base(data_in, size_in, name)
  {}
  KOKKOS_FORCEINLINE_FUNCTION void operator()(size_t i) const {
    Base::data[i] = Kokkos::reduction_identity<ValueType>::sum();
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
  ScatterAccess(view_type const& view_in)
    : view(view_in)
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
  inline
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
      original_view.extent(0),
      original_view.extent(1),
      original_view.extent(2),
      original_view.extent(3),
      original_view.extent(4),
      original_view.extent(5),
      original_view.extent(6),
      0
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
  ScatterView(std::string const& name, Dims ... dims)
  : internal_view(Kokkos::ViewAllocateWithoutInitializing(name), dims ..., unique_token.size())
  {
    reset();
  }

  template <int override_contribution = contribution>
  inline
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
  void contribute_into(View<DataType, RP...> const& dest) const
  {
    typedef View<DataType, RP...> dest_type;
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

  inline ScatterAccess(view_type const& view_in)
    : view(view_in)
    , thread_id(view_in.unique_token.acquire()) {
  }

  inline ~ScatterAccess() {
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
  inline ScatterAccess(ScatterAccess&& other)
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
