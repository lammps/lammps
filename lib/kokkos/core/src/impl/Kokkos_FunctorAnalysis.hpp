// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_FUNCTORANALYSIS_HPP
#define KOKKOS_FUNCTORANALYSIS_HPP

#include <cstddef>
#include <new>
#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class F, class WTag, class ref_type>
concept has_init_tag_function = requires(const F f) {
  { f.init(WTag{}, std::declval<ref_type>()) } -> std::same_as<void>;
};

template <class F, class ref_type>
concept has_init_no_tag_function = requires(const F f) {
  { f.init(std::declval<ref_type>()) } -> std::same_as<void>;
};

template <class F, class WTag, class ref_type, class cref_type>
concept has_join_tag_function = requires(const F f) {
  {
    f.join(WTag{}, std::declval<ref_type>(), std::declval<cref_type>())
  } -> std::same_as<void>;
};

template <class F, class ref_type, class cref_type>
concept has_join_no_tag_function = requires(const F f) {
  {
    f.join(std::declval<ref_type>(), std::declval<cref_type>())
  } -> std::same_as<void>;
};

template <class F, class WTag, class ref_type>
concept has_final_tag_function = requires(const F f) {
  { f.final(WTag{}, std::declval<ref_type>()) } -> std::same_as<void>;
};

template <class F, class ref_type>
concept has_final_no_tag_function = requires(const F f) {
  { f.final(std::declval<ref_type>()) } -> std::same_as<void>;
};

struct FunctorPatternInterface {
  struct FOR {};
  struct REDUCE {};
  struct SCAN {};
};

template <typename T>
struct DeduceFunctorPatternInterface;

template <class FunctorType, class ExecPolicy, class ExecutionSpace>
struct DeduceFunctorPatternInterface<
    ParallelFor<FunctorType, ExecPolicy, ExecutionSpace>> {
  using type = FunctorPatternInterface::FOR;
};

template <class CombinedFunctorReducerType, class ExecPolicy,
          class ExecutionSpace>
struct DeduceFunctorPatternInterface<
    ParallelReduce<CombinedFunctorReducerType, ExecPolicy, ExecutionSpace>> {
  using type = FunctorPatternInterface::REDUCE;
};

template <class FunctorType, class ExecPolicy, class ExecutionSpace>
struct DeduceFunctorPatternInterface<
    ParallelScan<FunctorType, ExecPolicy, ExecutionSpace>> {
  using type = FunctorPatternInterface::SCAN;
};

template <class FunctorType, class ExecPolicy, class ReturnType,
          class ExecutionSpace>
struct DeduceFunctorPatternInterface<ParallelScanWithTotal<
    FunctorType, ExecPolicy, ReturnType, ExecutionSpace>> {
  using type = FunctorPatternInterface::SCAN;
};

/** \brief  Query Functor and execution policy argument tag for value type.
 *
 *  If 'value_type' is not explicitly declared in the functor and
 * OverrideValueType is void, then attempt to deduce the type from
 * FunctorType::operator() interface used by the pattern and policy.
 *
 *  For the REDUCE pattern generate a Reducer and finalization function
 *  derived from what is available within the functor.
 */
template <typename PatternInterface, class Policy, class Functor,
          typename OverrideValueType>
struct FunctorAnalysis {
 private:
  using FOR    = FunctorPatternInterface::FOR;
  using REDUCE = FunctorPatternInterface::REDUCE;
  using SCAN   = FunctorPatternInterface::SCAN;

  //----------------------------------------

  struct void_tag {};

  template <typename P = Policy>
  struct has_work_tag {
    using type = void;
    using wtag = void_tag;
  };

  template <typename P>
    requires requires { typename P::work_tag; }
  struct has_work_tag<P> {
    using type = typename P::work_tag;
    using wtag = typename P::work_tag;
  };

  using Tag  = typename has_work_tag<>::type;
  using WTag = typename has_work_tag<>::wtag;

  //----------------------------------------
  // Check for T::execution_space

  template <typename T>
  struct has_execution_space {
    using type = void;
    enum : bool { value = false };
  };

  template <typename T>
    requires requires { typename T::execution_space; }
  struct has_execution_space<T> {
    using type = typename T::execution_space;
    enum : bool { value = true };
  };

  using policy_has_space  = has_execution_space<Policy>;
  using functor_has_space = has_execution_space<Functor>;

  static_assert(!policy_has_space::value || !functor_has_space::value ||
                    std::is_same_v<typename policy_has_space::type,
                                   typename functor_has_space::type>,
                "Execution Policy and Functor execution space must match");

  //----------------------------------------
  // Check for Functor::value_type, which is either a simple type T or T[]

  // If the functor doesn't have a value_type alias, use OverrideValueType.
  template <typename F>
  struct has_value_type {
    using type = OverrideValueType;
  };

  template <typename F>
    requires requires { typename F::value_type; }
  struct has_value_type<F> {
    using type = typename F::value_type;

    static_assert(!std::is_reference_v<type> && std::rank_v<type> <= 1 &&
                      std::extent_v<type> == 0,
                  "Kokkos Functor::value_type is T or T[]");
  };

  //----------------------------------------
  // If Functor::value_type does not exist and OverrideValueType is void, then
  // evaluate operator(), depending upon the pattern and whether the policy has
  // a work tag, to determine the reduction or scan value_type.

  template <typename F, typename P = PatternInterface,
            typename V = typename has_value_type<F>::type,
            bool T     = std::is_void_v<Tag>>
  struct deduce_value_type {
    using type = V;
  };

  template <typename F>
  struct deduce_value_type<F, REDUCE, void, true> {
    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, M, A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, M, M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, M, M, M, M, M,
                                                             A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, M, M, M, M, M,
                                                             M, A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, M, M, M, M, M,
                                                             M, M, A&) const);

    using type = decltype(deduce(&F::operator()));
  };

  template <typename F>
  struct deduce_value_type<F, REDUCE, void, false> {
    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, M, M, M,
                                                             A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, M, M, M,
                                                             M, A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, M, M, M,
                                                             M, M, A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, M, M, M,
                                                             M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, M, M, M,
                                                             M, M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, M,
                                                             A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, M,
                                                             M, A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, M,
                                                             M, M, A&) const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, M,
                                                             M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, M,
                                                             M, M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, M,
                                                             M, M, M, M, M, A&)
                                               const);

    template <typename M, typename A>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, M,
                                                             M, M, M, M, M, M,
                                                             A&) const);

    using type = decltype(deduce(&F::operator()));
  };

  template <typename F>
  struct deduce_value_type<F, SCAN, void, true> {
    template <typename M, typename A, typename I>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(M, A&, I) const);

    using type = decltype(deduce(&F::operator()));
  };

  template <typename F>
  struct deduce_value_type<F, SCAN, void, false> {
    template <typename M, typename A, typename I>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag, M, A&, I)
                                               const);

    template <typename M, typename A, typename I>
    KOKKOS_INLINE_FUNCTION static A deduce(void (Functor::*)(WTag const&, M, A&,
                                                             I) const);

    using type = decltype(deduce(&F::operator()));
  };

  //----------------------------------------

  using candidate_type = typename deduce_value_type<Functor>::type;

  enum {
    candidate_is_void  = std::is_void_v<candidate_type>,
    candidate_is_array = std::rank_v<candidate_type> == 1
  };

  //----------------------------------------

 public:
  using execution_space =
      std::conditional_t<functor_has_space::value,
                         typename functor_has_space::type,
                         std::conditional_t<policy_has_space::value,
                                            typename policy_has_space::type,
                                            Kokkos::DefaultExecutionSpace>>;

  using value_type = std::remove_extent_t<candidate_type>;

  static_assert(!std::is_const_v<value_type>,
                "Kokkos functor operator reduce argument cannot be const");

 private:
  // Stub to avoid defining a type 'void &'
  using ValueType = std::conditional_t<candidate_is_void, void_tag, value_type>;

 public:
  using pointer_type = std::conditional_t<candidate_is_void, void, ValueType*>;

  using reference_type = std::conditional_t<
      candidate_is_array, ValueType*,
      std::conditional_t<!candidate_is_void, ValueType&, void>>;

 private:
  template <bool IsArray, class FF>
  KOKKOS_INLINE_FUNCTION static constexpr unsigned int get_length(FF const& f) {
    if constexpr (IsArray)
      return f.value_count;
    else
      return candidate_is_void ? 0 : 1;
  }

 public:
  enum {
    StaticValueSize =
        !candidate_is_void && !candidate_is_array ? sizeof(ValueType) : 0
  };

  KOKKOS_FORCEINLINE_FUNCTION static constexpr unsigned int value_count(
      const Functor& f) {
    return FunctorAnalysis::template get_length<candidate_is_array>(f);
  }

  KOKKOS_FORCEINLINE_FUNCTION static constexpr unsigned int value_size(
      const Functor& f) {
    return FunctorAnalysis::template get_length<candidate_is_array>(f) *
           sizeof(ValueType);
  }

  //----------------------------------------

  template <class Unknown>
  KOKKOS_FORCEINLINE_FUNCTION static constexpr unsigned int value_count(
      const Unknown&) {
    return candidate_is_void ? 0 : 1;
  }

  template <class Unknown>
  KOKKOS_FORCEINLINE_FUNCTION static constexpr unsigned int value_size(
      const Unknown&) {
    return candidate_is_void ? 0 : sizeof(ValueType);
  }

 private:
  using ref_type =
      std::conditional_t<candidate_is_array, ValueType*, ValueType&>;
  using cref_type = std::conditional_t<candidate_is_array, const ValueType*,
                                       const ValueType&>;

  //----------------------------------------
  // parallel_reduce join operator

  template <class F = Functor>
  struct DeduceJoinNoTag {
    enum : bool { value = false };

    KOKKOS_INLINE_FUNCTION static void join(F const* const f, ValueType* dst,
                                            ValueType const* src) {
      const int n = FunctorAnalysis::value_count(*f);
      for (int i = 0; i < n; ++i) dst[i] += src[i];
    }
  };

  template <has_join_no_tag_function<ref_type, cref_type> F>
    requires Kokkos::Reducer<F> || (!Kokkos::Reducer<F> && std::is_void_v<Tag>)
  struct DeduceJoinNoTag<F> {
    enum : bool { value = true };

    KOKKOS_INLINE_FUNCTION static void join(F const* const f, ValueType* dst,
                                            ValueType const* src) {
      if constexpr (candidate_is_array)
        f->join(dst, src);
      else
        f->join(*dst, *src);
    }
  };

  template <class F = Functor>
  struct DeduceJoin : public DeduceJoinNoTag<F> {};

  template <has_join_tag_function<WTag, ref_type, cref_type> F>
    requires(!Kokkos::Reducer<F>)
  struct DeduceJoin<F> {
    enum : bool { value = true };

    KOKKOS_INLINE_FUNCTION static void join(F const* const f, ValueType* dst,
                                            ValueType const* src) {
      if constexpr (candidate_is_array)
        f->join(WTag(), dst, src);
      else
        f->join(WTag(), *dst, *src);
    }
  };

  //----------------------------------------

  template <class F = Functor>
  struct DeduceInitNoTag {
    enum : bool { value = false };

    KOKKOS_INLINE_FUNCTION static void init(F const* const f, ValueType* dst) {
      const int n = FunctorAnalysis::value_count(*f);
      for (int i = 0; i < n; ++i) new (&dst[i]) ValueType();
    }
  };

  template <has_init_no_tag_function<ref_type> F>
    requires Kokkos::Reducer<F> || (!Kokkos::Reducer<F> && std::is_void_v<Tag>)
  struct DeduceInitNoTag<F> {
    enum : bool { value = true };

    KOKKOS_INLINE_FUNCTION static void init(F const* const f, ValueType* dst) {
      if constexpr (candidate_is_array)
        f->init(dst);
      else
        f->init(*dst);
    }
  };

  template <class F = Functor>
  struct DeduceInit : public DeduceInitNoTag<F> {};

  template <has_init_tag_function<WTag, ref_type> F>
    requires Kokkos::Reducer<F>
  struct DeduceInit<F> {
    enum : bool { value = true };

    KOKKOS_INLINE_FUNCTION static void init(F const* const f, ValueType* dst) {
      if constexpr (candidate_is_array)
        f->init(WTag(), dst);
      else
        f->init(WTag(), *dst);
    }
  };

  //----------------------------------------

  template <class F = Functor>
  struct DeduceFinalNoTag {
    enum : bool { value = false };

    KOKKOS_INLINE_FUNCTION
    static void final(F const* const, ValueType*) {}
  };

  template <has_final_no_tag_function<ref_type> F>
    requires Kokkos::Reducer<F> || (!Kokkos::Reducer<F> && std::is_void_v<Tag>)
  struct DeduceFinalNoTag<F> {
    enum : bool { value = true };

    KOKKOS_INLINE_FUNCTION static void final(F const* const f, ValueType* dst) {
      if constexpr (candidate_is_array)
        f->final(dst);
      else
        f->final(*dst);
    }
  };

  template <class F = Functor>
  struct DeduceFinal : public DeduceFinalNoTag<F> {};

  template <has_final_tag_function<WTag, ref_type> F>
    requires(!Kokkos::Reducer<F>)
  struct DeduceFinal<F> {
    enum : bool { value = true };

    KOKKOS_INLINE_FUNCTION static void final(F const* const f, ValueType* dst) {
      if constexpr (candidate_is_array)
        f->final(WTag(), dst);
      else
        f->final(WTag(), *dst);
    }
  };

 public:
  enum { has_join_member_function = DeduceJoin<>::value };
  enum { has_init_member_function = DeduceInit<>::value };
  enum { has_final_member_function = DeduceFinal<>::value };

  static_assert((Kokkos::is_reducer<Functor>::value &&
                 has_join_member_function) ||
                    !Kokkos::is_reducer<Functor>::value,
                "Reducer must have a join member function!");

  struct Reducer {
   private:
    Functor m_functor;

   public:
    using reducer        = Reducer;
    using value_type     = std::remove_const_t<FunctorAnalysis::value_type>;
    using pointer_type   = value_type*;
    using reference_type = FunctorAnalysis::reference_type;
    using functor_type   = Functor;  // Adapts a functor

    static constexpr bool has_join_member_function() {
      return DeduceJoin<>::value;
    }
    static constexpr bool has_init_member_function() {
      return DeduceInit<>::value;
    }
    static constexpr bool has_final_member_function() {
      return DeduceFinal<>::value;
    }

    KOKKOS_FUNCTION unsigned int value_size() const {
      return FunctorAnalysis::value_size(m_functor);
    }

    KOKKOS_FUNCTION unsigned int value_count() const {
      return FunctorAnalysis::value_count(m_functor);
    }

    KOKKOS_FUNCTION static constexpr unsigned int static_value_size() {
      return StaticValueSize;
    }

    KOKKOS_INLINE_FUNCTION static reference_type reference(
        ValueType* dst) noexcept {
      if constexpr (candidate_is_array)
        return dst;
      else
        return *dst;
    }

    KOKKOS_INLINE_FUNCTION constexpr int length() const noexcept {
      if constexpr (candidate_is_array)
        return m_functor.value_count;
      else
        return candidate_is_void ? 0 : 1;
    }

    KOKKOS_INLINE_FUNCTION
    void copy(ValueType* const dst, ValueType const* const src) const noexcept {
      for (int i = 0; i < length(); ++i) dst[i] = src[i];
    }

    KOKKOS_INLINE_FUNCTION
    void join(ValueType* dst, ValueType const* src) const noexcept {
      DeduceJoin<>::join(&m_functor, dst, src);
    }

    KOKKOS_INLINE_FUNCTION reference_type
    init(ValueType* const dst) const noexcept {
      DeduceInit<>::init(&m_functor, dst);
      return reference(dst);
    }

    KOKKOS_INLINE_FUNCTION
    void final(ValueType* dst) const noexcept {
      DeduceFinal<>::final(&m_functor, dst);
    }

    KOKKOS_INLINE_FUNCTION
    const Functor& get_functor() const { return m_functor; }

    Reducer(Reducer const&)            = default;
    Reducer(Reducer&&)                 = default;
    Reducer& operator=(Reducer const&) = delete;
    Reducer& operator=(Reducer&&)      = delete;
    ~Reducer()                         = default;

    KOKKOS_INLINE_FUNCTION explicit constexpr Reducer(
        Functor const& arg_functor) noexcept
        : m_functor(arg_functor) {}
  };
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_FUNCTORANALYSIS_HPP */
