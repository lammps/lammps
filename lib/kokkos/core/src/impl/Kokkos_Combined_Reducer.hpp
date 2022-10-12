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

#ifndef KOKKOS_COMBINED_REDUCER_HPP
#define KOKKOS_COMBINED_REDUCER_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_Parallel_Reduce.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_AnonymousSpace.hpp>
#include <impl/Kokkos_Utilities.hpp>  // comma operator fold emulation

#include <utility>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="CombinedReducer reducer and value storage helpers"> {{{1

// Note: the index is only to avoid repeating the same base class multiple times
template <size_t /*Idx*/, class ValueType>
struct CombinedReducerValueItemImpl {
 public:
  using value_type = ValueType;

 private:
  value_type m_value;

 public:
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerValueItemImpl() = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerValueItemImpl(
      CombinedReducerValueItemImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerValueItemImpl(
      CombinedReducerValueItemImpl&&) = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerValueItemImpl& operator=(
      CombinedReducerValueItemImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerValueItemImpl& operator=(
      CombinedReducerValueItemImpl&&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  ~CombinedReducerValueItemImpl() = default;
  explicit KOKKOS_FUNCTION CombinedReducerValueItemImpl(value_type arg_value)
      : m_value(std::move(arg_value)) {}

  KOKKOS_FORCEINLINE_FUNCTION
  constexpr value_type& ref() & noexcept { return m_value; }
  KOKKOS_FORCEINLINE_FUNCTION
  constexpr value_type const& ref() const& noexcept { return m_value; }
};

//==============================================================================

template <class IdxSeq, class... ValueTypes>
struct CombinedReducerValueImpl;

template <size_t... Idxs, class... ValueTypes>
struct CombinedReducerValueImpl<std::integer_sequence<size_t, Idxs...>,
                                ValueTypes...>
    : CombinedReducerValueItemImpl<Idxs, ValueTypes>... {
 public:
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReducerValueImpl() = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReducerValueImpl(CombinedReducerValueImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReducerValueImpl(CombinedReducerValueImpl&&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReducerValueImpl& operator=(
      CombinedReducerValueImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReducerValueImpl& operator=(CombinedReducerValueImpl&&) =
      default;
  KOKKOS_DEFAULTED_FUNCTION
  ~CombinedReducerValueImpl() = default;

  KOKKOS_FUNCTION
  explicit CombinedReducerValueImpl(ValueTypes... arg_values)
      : CombinedReducerValueItemImpl<Idxs, ValueTypes>(
            std::move(arg_values))... {}

  template <size_t Idx, class ValueType>
      KOKKOS_INLINE_FUNCTION ValueType& get() & noexcept {
    return this->CombinedReducerValueItemImpl<Idx, ValueType>::ref();
  }
  template <size_t Idx, class ValueType>
  KOKKOS_INLINE_FUNCTION ValueType const& get() const& noexcept {
    return this->CombinedReducerValueItemImpl<Idx, ValueType>::ref();
  }
};

//==============================================================================

// TODO Empty base optmization?
template <size_t /*Idx*/, class Reducer>
// requires Kokkos::is_reducer<Reducer>
struct CombinedReducerStorageImpl {
 public:
  using value_type = typename Reducer::value_type;

 private:
  Reducer m_reducer;

 public:
  KOKKOS_INLINE_FUNCTION
  explicit constexpr CombinedReducerStorageImpl(Reducer arg_reducer)
      : m_reducer(std::move(arg_reducer)) {}

  // Leading underscores to make it clear that this class is not intended to
  // model Reducer

  KOKKOS_INLINE_FUNCTION
  constexpr _fold_comma_emulation_return _init(value_type& val) const {
    m_reducer.init(val);
    return _fold_comma_emulation_return{};
  }

  KOKKOS_INLINE_FUNCTION constexpr _fold_comma_emulation_return _join(
      value_type& dest, value_type const& src) const {
    m_reducer.join(dest, src);
    return _fold_comma_emulation_return{};
  }
};

// </editor-fold> end CombinedReducerStorage }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="CombinedReducer"> {{{1

struct _construct_combined_reducer_from_args_tag {};

template <class T>
KOKKOS_INLINE_FUNCTION auto _get_value_from_combined_reducer_ctor_arg(
    T&& arg) noexcept
    -> std::enable_if_t<!is_view<std::decay_t<T>>::value &&
                            !is_reducer<std::decay_t<T>>::value,
                        std::decay_t<T>> {
  return arg;
}

template <class T>
KOKKOS_INLINE_FUNCTION auto _get_value_from_combined_reducer_ctor_arg(
    T&&) noexcept ->
    typename std::enable_if_t<is_view<std::decay_t<T>>::value ||
                                  is_reducer<std::decay_t<T>>::value,
                              std::decay_t<T>>::value_type {
  return typename std::decay_t<T>::value_type{};
}

template <class IdxSeq, class Space, class...>
struct CombinedReducerImpl;

template <size_t... Idxs, class Space, class... Reducers>
struct CombinedReducerImpl<std::integer_sequence<size_t, Idxs...>, Space,
                           Reducers...>
    : private CombinedReducerStorageImpl<Idxs, Reducers>... {
 public:
  using reducer = CombinedReducerImpl<std::integer_sequence<size_t, Idxs...>,
                                      Space, Reducers...>;
  using value_type =
      CombinedReducerValueImpl<std::integer_sequence<size_t, Idxs...>,
                               typename Reducers::value_type...>;
  using result_view_type =
      Kokkos::View<value_type, Space, Kokkos::MemoryUnmanaged>;

 private:
  result_view_type m_value_view;

 public:
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerImpl() = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerImpl(
      CombinedReducerImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerImpl(
      CombinedReducerImpl&&)                                       = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerImpl& operator=(
      CombinedReducerImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION constexpr CombinedReducerImpl& operator=(
      CombinedReducerImpl&&) = default;

  KOKKOS_DEFAULTED_FUNCTION ~CombinedReducerImpl() = default;

  template <class... ReducersDeduced>
  KOKKOS_FUNCTION constexpr explicit CombinedReducerImpl(
      value_type& value, ReducersDeduced&&... reducers) noexcept
      : CombinedReducerStorageImpl<Idxs, Reducers>((ReducersDeduced &&)
                                                       reducers)...,
        m_value_view(&value) {}

  KOKKOS_FUNCTION constexpr void join(value_type& dest,
                                      value_type const& src) const noexcept {
    emulate_fold_comma_operator(
        this->CombinedReducerStorageImpl<Idxs, Reducers>::_join(
            dest.template get<Idxs, typename Reducers::value_type>(),
            src.template get<Idxs, typename Reducers::value_type>())...);
  }

  KOKKOS_FUNCTION constexpr void init(value_type& dest) const noexcept {
    emulate_fold_comma_operator(
        this->CombinedReducerStorageImpl<Idxs, Reducers>::_init(
            dest.template get<Idxs, typename Reducers::value_type>())...);
  }

  // TODO figure out if we also need to call through to final

  KOKKOS_FUNCTION
  constexpr bool references_scalar() const noexcept {
    // For now, always pretend that we reference a scalar since we need to
    // block to do the write-back because the references may not be contiguous
    // in memory and the backends currently assume this and just do a single
    // deep copy back to a chunk of memory associated with the output argument
    return true;
  }

  KOKKOS_FUNCTION
  constexpr result_view_type const& view() const noexcept {
    return m_value_view;
  }

  template <class ExecutionSpace, int Idx, class View>
  static void write_one_value_back(
      const ExecutionSpace& exec_space, View const& view,
      typename View::const_value_type& value) noexcept {
    if (Kokkos::SpaceAccessibility<typename View::memory_space,
                                   Space>::assignable)
      view() = value;
    else
      Kokkos::deep_copy(exec_space, view, value);
  }

  template <class ExecutionSpace>
  static void write_value_back_to_original_references(
      const ExecutionSpace& exec_space, value_type const& value,
      Reducers const&... reducers_that_reference_original_values) noexcept {
    emulate_fold_comma_operator(
        (write_one_value_back<ExecutionSpace, Idxs>(
             exec_space, reducers_that_reference_original_values.view(),
             value.template get<Idxs, typename Reducers::value_type>()),
         0)...);
  }
};

// Apparently this can't be an alias template because of a bug/unimplemented
// feature in GCC's name mangler.  But in this case, this amounts to the same
// thing.
template <class Space, class... Reducers>
struct CombinedReducer
    : CombinedReducerImpl<std::make_index_sequence<sizeof...(Reducers)>, Space,
                          Reducers...> {
  using base_t =
      CombinedReducerImpl<std::make_index_sequence<sizeof...(Reducers)>, Space,
                          Reducers...>;
  using base_t::base_t;
  using reducer = CombinedReducer<Space, Reducers...>;
};

// </editor-fold> end CombinedReducer }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="CombinedReductionFunctorWrapper"> {{{1

template <class IdxSeq, class Functor, class Space, class... Reducers>
struct CombinedReductionFunctorWrapperImpl;

template <size_t... Idxs, class Functor, class Space, class... Reducers>
struct CombinedReductionFunctorWrapperImpl<
    std::integer_sequence<size_t, Idxs...>, Functor, Space, Reducers...> {
 private:
  Functor m_functor;

 public:
  //------------------------------------------------------------------------------
  // <editor-fold desc="type aliases"> {{{2

  using reducer_type = CombinedReducer<Space, Reducers...>;

  // Prevent Kokkos from attempting to deduce value_type
  using value_type = typename reducer_type::value_type;

  // </editor-fold> end type aliases }}}2
  //------------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="Ctors, destructor, and assignment"> {{{2

  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReductionFunctorWrapperImpl() noexcept = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReductionFunctorWrapperImpl(
      CombinedReductionFunctorWrapperImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReductionFunctorWrapperImpl(
      CombinedReductionFunctorWrapperImpl&&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReductionFunctorWrapperImpl& operator=(
      CombinedReductionFunctorWrapperImpl const&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReductionFunctorWrapperImpl& operator=(
      CombinedReductionFunctorWrapperImpl&&) = default;
  KOKKOS_DEFAULTED_FUNCTION
  ~CombinedReductionFunctorWrapperImpl() = default;

  KOKKOS_INLINE_FUNCTION
  constexpr explicit CombinedReductionFunctorWrapperImpl(Functor arg_functor)
      : m_functor(std::move(arg_functor)) {}

  // </editor-fold> end Ctors, destructor, and assignment }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="call operator"> {{{2

  // Variadic version for MDRangePolicy
  // There are a number of ways to do this, but most of them that involve
  // not assuming an implementation of tuple is available are gross.
  // Unfortunately, that's what we have to do here
  template <class IndexOrMemberOrTagType1,
            class... IndexOrMemberTypesThenValueType>
  KOKKOS_FUNCTION void operator()(
      IndexOrMemberOrTagType1&& arg_first,
      IndexOrMemberTypesThenValueType&&... args) const {
    this->template _call_op_impl<IndexOrMemberOrTagType1&&>(
        (IndexOrMemberOrTagType1 &&) arg_first,
        (IndexOrMemberTypesThenValueType &&) args...);
  }

  // </editor-fold> end call operator }}}2
  //----------------------------------------------------------------------------

  // These are things that need to be done if we decide to ever support
  // functor-customized join/init/final hooks with combined reducers. For now,
  // they are explicitly not supported.
  // TODO: forward join() function to user functor hook, or just ignore it?
  // TODO: forward init() function to user functor hook, or just ignore it?
  // TODO: forward final() function to user functor hook, or just ignore it?

 private:
  // variadic forwarding for MDRangePolicy
  // see comment above for why this has to be so gross
  // recursive case
  template <class... IdxOrMemberTypes, class IdxOrMemberType1,
            class... IdxOrMemberTypesThenValueType>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      !std::is_same<remove_cvref_t<IdxOrMemberType1>, value_type>::value>
  _call_op_impl(IdxOrMemberTypes&&... idxs, IdxOrMemberType1&& idx,
                IdxOrMemberTypesThenValueType&&... args) const {
    this->template _call_op_impl<IdxOrMemberTypes&&..., IdxOrMemberType1&&>(
        (IdxOrMemberTypes &&) idxs..., (IdxOrMemberType1 &&) idx,
        (IdxOrMemberTypesThenValueType &&) args...);
  }

  // base case
  template <class... IdxOrMemberTypes>
  KOKKOS_FORCEINLINE_FUNCTION void _call_op_impl(IdxOrMemberTypes&&... idxs,
                                                 value_type& out) const {
    m_functor((IdxOrMemberTypes &&) idxs...,
              out.template get<Idxs, typename Reducers::value_type>()...);
  }
};

template <class Functor, class Space, class... Reducers>
struct CombinedReductionFunctorWrapper
    : CombinedReductionFunctorWrapperImpl<
          std::make_index_sequence<sizeof...(Reducers)>, Functor, Space,
          Reducers...> {
  using base_t = CombinedReductionFunctorWrapperImpl<
      std::make_index_sequence<sizeof...(Reducers)>, Functor, Space,
      Reducers...>;
  using base_t::base_t;
};

// </editor-fold> end CombinedReductionFunctorWrapper }}}1
//==============================================================================

//------------------------------------------------------------------------------
// <editor-fold desc="_make_reducer_from_arg"> {{{2

template <class Space, class Reducer>
KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<
    Kokkos::is_reducer<std::decay_t<Reducer>>::value, std::decay_t<Reducer>>
_make_reducer_from_arg(Reducer&& arg_reducer) noexcept {
  return arg_reducer;
}

// Two purposes: SFINAE-safety for the `View` case and laziness for the return
// value otherwise to prevent extra instantiations of the Kokkos::Sum template
template <class Space, class T, class Enable = void>
struct _wrap_with_kokkos_sum {
  using type = Kokkos::Sum<T, Space>;
};

template <class Space, class T>
struct _wrap_with_kokkos_sum<Space, T,
                             std::enable_if_t<Kokkos::is_view<T>::value>> {
  using type = Kokkos::Sum<typename T::value_type, typename T::memory_space>;
};

// TODO better error message for the case when a const& to a scalar is passed in
//      (this is needed in general, though)
template <class Space, class T>
KOKKOS_INLINE_FUNCTION constexpr typename std::enable_if_t<
    !Kokkos::is_reducer<std::decay_t<T>>::value,
    _wrap_with_kokkos_sum<Space, std::decay_t<T>>>::type
_make_reducer_from_arg(T&& arg_scalar) noexcept {
  return
      typename _wrap_with_kokkos_sum<Space, std::decay_t<T>>::type{arg_scalar};
}

// This can't be an alias template because GCC doesn't know how to mangle
// decltype expressions in return statements (and, even though every compiler
// is supposed to, GCC is the only one that does dependent alias template
// substitution correctly and tries to do the mangling, aparently).
template <class Space, class ReferenceOrViewOrReducer, class = void>
struct _reducer_from_arg {
  using type = decltype(Impl::_make_reducer_from_arg<Space>(
      std::declval<ReferenceOrViewOrReducer&&>()));
};
template <class Space, class ReferenceOrViewOrReducer>
using _reducer_from_arg_t =
    typename _reducer_from_arg<Space, ReferenceOrViewOrReducer>::type;

// </editor-fold> end _make_reducer_from_arg }}}2
//------------------------------------------------------------------------------

template <class Space, class... ReferencesOrViewsOrReducers>
KOKKOS_INLINE_FUNCTION constexpr auto make_combined_reducer_value(
    ReferencesOrViewsOrReducers&&... args) {
  //----------------------------------------
  // This is a bit round-about and we should make sure it doesn't have
  // any performance implications. Basically, we make a reducer out of anything
  // just to get the value back out here (for the sake of uniformity). Most
  // compilers should figure out what's going on, but we should double-check
  // that.
  return CombinedReducerValueImpl<
      std::make_index_sequence<sizeof...(ReferencesOrViewsOrReducers)>,
      typename _reducer_from_arg_t<Space,
                                   ReferencesOrViewsOrReducers>::value_type...>{
      // This helper function is now poorly named after refactoring.
      _get_value_from_combined_reducer_ctor_arg((ReferencesOrViewsOrReducers &&)
                                                    args)...};
  //----------------------------------------
}

template <class Space, class ValueType, class... ReferencesOrViewsOrReducers>
KOKKOS_INLINE_FUNCTION constexpr auto make_combined_reducer(
    ValueType& value, ReferencesOrViewsOrReducers&&... args) {
  //----------------------------------------
  // This is doing more or less the same thing of making every argument into
  // a reducer, just in a different place than in `make_combined_reducer_value`,
  // so we should probably eventually make this read a little more similarly
  using reducer_type = CombinedReducer<
      Space, _reducer_from_arg_t<Space, ReferencesOrViewsOrReducers>...>;
  return reducer_type(value,
                      _reducer_from_arg_t<Space, ReferencesOrViewsOrReducers>{
                          (ReferencesOrViewsOrReducers &&) args}...);
  //----------------------------------------
}

template <class Functor, class Space, class... ReferencesOrViewsOrReducers>
KOKKOS_INLINE_FUNCTION constexpr auto make_wrapped_combined_functor(
    Functor const& functor, Space, ReferencesOrViewsOrReducers&&...) {
  //----------------------------------------
  return CombinedReductionFunctorWrapper<
      Functor, Space,
      _reducer_from_arg_t<Space, ReferencesOrViewsOrReducers>...>(functor);
  //----------------------------------------
}

template <typename FunctorType>
using functor_has_value_t = typename FunctorType::value_type;
}  // end namespace Impl

//==============================================================================
// <editor-fold desc="Overloads of parallel_reduce for multiple outputs"> {{{1

// These need to be forwarding references so that we can deduce const-ness,
// but none of them should be forwarded (and, indeed, none of them should be
// rvalue references)
template <class PolicyType, class Functor, class ReturnType1, class ReturnType2,
          class... ReturnTypes>
auto parallel_reduce(std::string const& label, PolicyType const& policy,
                     Functor const& functor, ReturnType1&& returnType1,
                     ReturnType2&& returnType2,
                     ReturnTypes&&... returnTypes) noexcept
    -> std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value> {
  //----------------------------------------
  // Since we don't support asynchronous combined reducers yet for various
  // reasons, we actually just want to work with the pointers and references
  // directly
  using space_type = Kokkos::DefaultHostExecutionSpace::memory_space;

  auto value = Impl::make_combined_reducer_value<space_type>(
      returnType1, returnType2, returnTypes...);

  using combined_reducer_type = Impl::CombinedReducer<
      space_type, Impl::_reducer_from_arg_t<space_type, ReturnType1>,
      Impl::_reducer_from_arg_t<space_type, ReturnType2>,
      Impl::_reducer_from_arg_t<space_type, ReturnTypes>...>;
  auto combined_reducer = Impl::make_combined_reducer<space_type>(
      value, returnType1, returnType2, returnTypes...);

  auto combined_functor = Impl::make_wrapped_combined_functor(
      functor, space_type{}, returnType1, returnType2, returnTypes...);

  using combined_functor_type = decltype(combined_functor);
  static_assert(
      is_detected<Impl::functor_has_value_t, combined_functor_type>::value,
      "value_type not properly detected");
  using reduce_adaptor_t =
      Impl::ParallelReduceAdaptor<PolicyType, combined_functor_type,
                                  combined_reducer_type>;

  reduce_adaptor_t::execute(label, policy, combined_functor, combined_reducer);
  Impl::ParallelReduceFence<typename PolicyType::execution_space,
                            combined_reducer_type>::
      fence(
          policy.space(),
          "Kokkos::parallel_reduce: fence due to result being value, not view",
          combined_reducer);
  combined_reducer.write_value_back_to_original_references(
      policy.space(), value,
      Impl::_make_reducer_from_arg<space_type>(returnType1),
      Impl::_make_reducer_from_arg<space_type>(returnType2),
      Impl::_make_reducer_from_arg<space_type>(returnTypes)...);
  policy.space().fence(
      "Kokkos::parallel_reduce: fence after copying values back");
  //----------------------------------------
}

template <class PolicyType, class Functor, class ReturnType1, class ReturnType2,
          class... ReturnTypes>
auto parallel_reduce(PolicyType const& policy, Functor const& functor,
                     ReturnType1&& returnType1, ReturnType2&& returnType2,
                     ReturnTypes&&... returnTypes) noexcept
    -> std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value> {
  //----------------------------------------
  Kokkos::parallel_reduce("", policy, functor,
                          std::forward<ReturnType1>(returnType1),
                          std::forward<ReturnType2>(returnType2),
                          std::forward<ReturnTypes>(returnTypes)...);
  //----------------------------------------
}

template <class Functor, class ReturnType1, class ReturnType2,
          class... ReturnTypes>
void parallel_reduce(std::string const& label, size_t n, Functor const& functor,
                     ReturnType1&& returnType1, ReturnType2&& returnType2,
                     ReturnTypes&&... returnTypes) noexcept {
  Kokkos::parallel_reduce(label,
                          RangePolicy<Kokkos::DefaultExecutionSpace>(0, n),
                          functor, std::forward<ReturnType1>(returnType1),
                          std::forward<ReturnType2>(returnType2),
                          std::forward<ReturnTypes>(returnTypes)...);
}

template <class Functor, class ReturnType1, class ReturnType2,
          class... ReturnTypes>
void parallel_reduce(size_t n, Functor const& functor,
                     ReturnType1&& returnType1, ReturnType2&& returnType2,
                     ReturnTypes&&... returnTypes) noexcept {
  Kokkos::parallel_reduce("", n, functor,
                          std::forward<ReturnType1>(returnType1),
                          std::forward<ReturnType2>(returnType2),
                          std::forward<ReturnTypes>(returnTypes)...);
}

//------------------------------------------------------------------------------
// <editor-fold desc="Team overloads"> {{{2

// Copied three times because that's the best way we have right now to match
// Impl::TeamThreadRangeBoundariesStruct,
// Impl::ThreadVectorRangeBoundariesStruct, and
// Impl::TeamVectorRangeBoundariesStruct.
// TODO make these work after restructuring

// template <class iType, class MemberType, class Functor, class ReturnType1,
//          class ReturnType2, class... ReturnTypes>
// KOKKOS_INLINE_FUNCTION void parallel_reduce(
//    std::string const& label,
//    Impl::TeamThreadRangeBoundariesStruct<iType, MemberType> const&
//    boundaries, Functor const& functor, ReturnType1&& returnType1,
//    ReturnType2&& returnType2, ReturnTypes&&... returnTypes) noexcept {
//  const auto combined_reducer =
//      Impl::make_combined_reducer<Kokkos::AnonymousSpace>(
//          returnType1, returnType2, returnTypes...);
//
//  auto combined_functor = Impl::make_wrapped_combined_functor(
//      functor, Kokkos::AnonymousSpace{}, returnType1, returnType2,
//      returnTypes...);
//
//  parallel_reduce(label, boundaries, combined_functor, combined_reducer);
//}
//
// template <class iType, class MemberType, class Functor, class ReturnType1,
//          class ReturnType2, class... ReturnTypes>
// KOKKOS_INLINE_FUNCTION void parallel_reduce(
//    std::string const& label,
//    Impl::ThreadVectorRangeBoundariesStruct<iType, MemberType> const&
//        boundaries,
//    Functor const& functor, ReturnType1&& returnType1,
//    ReturnType2&& returnType2, ReturnTypes&&... returnTypes) noexcept {
//  const auto combined_reducer =
//      Impl::make_combined_reducer<Kokkos::AnonymousSpace>(
//          returnType1, returnType2, returnTypes...);
//
//  auto combined_functor = Impl::make_wrapped_combined_functor(
//      functor, Kokkos::AnonymousSpace{}, returnType1, returnType2,
//      returnTypes...);
//
//  parallel_reduce(label, boundaries, combined_functor, combined_reducer);
//}

// template <class iType, class MemberType, class Functor, class ReturnType1,
//          class ReturnType2, class... ReturnTypes>
// KOKKOS_INLINE_FUNCTION void parallel_reduce(
//    std::string const& label,
//    Impl::TeamVectorRangeBoundariesStruct<iType, MemberType> const&
//    boundaries, Functor const& functor, ReturnType1&& returnType1,
//    ReturnType2&& returnType2, ReturnTypes&&... returnTypes) noexcept {
//  const auto combined_reducer =
//      Impl::make_combined_reducer<Kokkos::AnonymousSpace>(
//          returnType1, returnType2, returnTypes...);
//
//  auto combined_functor = Impl::make_wrapped_combined_functor(
//      functor, Kokkos::AnonymousSpace{}, returnType1, returnType2,
//      returnTypes...);
//
//  parallel_reduce(label, boundaries, combined_functor, combined_reducer);
//}

// </editor-fold> end Team overloads }}}2
//------------------------------------------------------------------------------

// </editor-fold> end Overloads of parallel_reduce for multiple outputs }}}1
//==============================================================================

}  // namespace Kokkos

#endif  // KOKKOS_COMBINED_REDUCER_HPP
