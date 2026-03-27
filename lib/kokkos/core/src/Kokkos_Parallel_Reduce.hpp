// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_PARALLEL_REDUCE_HPP
#define KOKKOS_PARALLEL_REDUCE_HPP

#include <impl/Kokkos_BuiltinReducers.hpp>
#include <Kokkos_CheckUsage.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_View.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>
#include <impl/Kokkos_Tools_Generic.hpp>

#include <type_traits>

namespace Kokkos {
namespace Impl {

template <typename FunctorType, typename FunctorAnalysisReducerType,
          typename Enable>
class CombinedFunctorReducer {
 public:
  using functor_type = FunctorType;
  using reducer_type = FunctorAnalysisReducerType;
  CombinedFunctorReducer(const FunctorType& functor,
                         const FunctorAnalysisReducerType& reducer)
      : m_functor(functor), m_reducer(reducer) {}
  KOKKOS_FUNCTION const FunctorType& get_functor() const { return m_functor; }
  KOKKOS_FUNCTION const FunctorAnalysisReducerType& get_reducer() const {
    return m_reducer;
  }

 private:
  FunctorType m_functor;
  FunctorAnalysisReducerType m_reducer;
};
template <typename FunctorType, typename FunctorAnalysisReducerType>
class CombinedFunctorReducer<
    FunctorType, FunctorAnalysisReducerType,
    std::enable_if_t<std::is_same_v<
        FunctorType, typename FunctorAnalysisReducerType::functor_type>>> {
 public:
  using functor_type = FunctorType;
  using reducer_type = FunctorAnalysisReducerType;
  CombinedFunctorReducer(const FunctorType& functor,
                         const FunctorAnalysisReducerType&)
      : m_reducer(functor) {}
  KOKKOS_FUNCTION const FunctorType& get_functor() const {
    return m_reducer.get_functor();
  }
  KOKKOS_FUNCTION const FunctorAnalysisReducerType& get_reducer() const {
    return m_reducer;
  }

 private:
  FunctorAnalysisReducerType m_reducer;
};

template <class T, class ReturnType, class ValueTraits>
struct ParallelReduceReturnValue;

template <class ReturnType, class FunctorType>
struct ParallelReduceReturnValue<
    std::enable_if_t<Kokkos::is_view<ReturnType>::value>, ReturnType,
    FunctorType> {
  using return_type  = ReturnType;
  using reducer_type = InvalidType;

  using value_type_scalar = typename return_type::value_type;
  using value_type_array  = typename return_type::value_type* const;

  using value_type = std::conditional_t<return_type::rank == 0,
                                        value_type_scalar, value_type_array>;

  static return_type& return_value(ReturnType& return_val, const FunctorType&) {
    return return_val;  // NOLINT(bugprone-return-const-ref-from-parameter)
  }
};

template <class ReturnType, class FunctorType>
struct ParallelReduceReturnValue<
    std::enable_if_t<!Kokkos::is_view<ReturnType>::value &&
                     (!std::is_array_v<ReturnType> &&
                      !std::is_pointer_v<
                          ReturnType>)&&!Kokkos::is_reducer<ReturnType>::value>,
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
    std::enable_if_t<(std::is_array_v<ReturnType> ||
                      std::is_pointer_v<ReturnType>)>,
    ReturnType, FunctorType> {
  using return_type = Kokkos::View<std::remove_const_t<ReturnType>,
                                   Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

  using reducer_type = InvalidType;

  using value_type = typename return_type::value_type[];

  static return_type return_value(ReturnType& return_val,
                                  const FunctorType& functor) {
    if (std::is_array_v<ReturnType>)
      return return_type(return_val);
    else
      return return_type(return_val, functor.value_count);
  }
};

template <class ReturnType, class FunctorType>
struct ParallelReduceReturnValue<
    std::enable_if_t<Kokkos::is_reducer<ReturnType>::value>, ReturnType,
    FunctorType> {
  using return_type  = typename ReturnType::result_view_type;
  using reducer_type = ReturnType;
  using value_type   = typename return_type::value_type;

  static auto return_value(ReturnType& return_val, const FunctorType&) {
    return return_val.view();
  }
};

template <class T, class ReturnType, class FunctorType>
struct ParallelReducePolicyType;

template <class PolicyType, class FunctorType>
struct ParallelReducePolicyType<
    std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value>,
    PolicyType, FunctorType> {
  using policy_type = PolicyType;
  static PolicyType policy(const PolicyType& policy_) { return policy_; }
};

template <class PolicyType, class FunctorType>
struct ParallelReducePolicyType<
    std::enable_if_t<std::is_integral_v<PolicyType>>, PolicyType, FunctorType> {
  using execution_space =
      typename Impl::FunctorPolicyExecutionSpace<FunctorType,
                                                 void>::execution_space;

  using policy_type = Kokkos::RangePolicy<execution_space>;

  static policy_type policy(const PolicyType& policy_) {
    return policy_type(0, policy_);
  }
};

template <class PolicyType, class FunctorType, class ReturnType>
struct ParallelReduceAdaptor {
  using return_value_adapter =
      Impl::ParallelReduceReturnValue<void, ReturnType, FunctorType>;

  // Equivalent to std::get<I>(std::tuple) but callable on the device.
  template <bool B, class T1, class T2>
  static KOKKOS_FUNCTION std::conditional_t<B, T1&&, T2&&> forwarding_switch(
      T1&& v1, T2&& v2) {
    if constexpr (B)
      return static_cast<T1&&>(v1);
    else
      return static_cast<T2&&>(v2);
  }

  static inline void execute_impl(const std::string& label,
                                  const PolicyType& policy,
                                  const FunctorType& functor,
                                  ReturnType& return_value) {
    using PassedReducerType = typename return_value_adapter::reducer_type;
    uint64_t kpID           = 0;

    constexpr bool passed_reducer_type_is_invalid =
        std::is_same_v<InvalidType, PassedReducerType>;
    using TheReducerType = std::conditional_t<passed_reducer_type_is_invalid,
                                              FunctorType, PassedReducerType>;

    using Analysis = FunctorAnalysis<FunctorPatternInterface::REDUCE,
                                     PolicyType, TheReducerType,
                                     typename return_value_adapter::value_type>;
    using CombinedFunctorReducerType =
        CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>;

    CombinedFunctorReducerType functor_reducer(
        functor, typename Analysis::Reducer(
                     forwarding_switch<passed_reducer_type_is_invalid>(
                         functor, return_value)));
    const auto& response = Kokkos::Tools::Impl::begin_parallel_reduce<
        typename return_value_adapter::reducer_type>(policy, functor_reducer,
                                                     label, kpID);
    const auto& inner_policy = response.policy;

    auto closure = construct_with_shared_allocation_tracking_disabled<
        Impl::ParallelReduce<CombinedFunctorReducerType, PolicyType,
                             typename Impl::FunctorPolicyExecutionSpace<
                                 FunctorType, PolicyType>::execution_space>>(
        functor_reducer, inner_policy,
        return_value_adapter::return_value(return_value, functor));
    closure.execute();

    Kokkos::Tools::Impl::end_parallel_reduce<PassedReducerType>(
        inner_policy, functor, label, kpID);
  }

  static constexpr bool is_array_reduction =
      Impl::FunctorAnalysis<
          Impl::FunctorPatternInterface::REDUCE, PolicyType, FunctorType,
          typename return_value_adapter::value_type>::StaticValueSize == 0;

  template <typename Dummy = ReturnType>
  static inline std::enable_if_t<!(is_array_reduction &&
                                   std::is_pointer_v<Dummy>)>
  execute(const std::string& label, const PolicyType& policy,
          const FunctorType& functor, ReturnType& return_value) {
    execute_impl(label, policy, functor, return_value);
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
    value = std::is_same_v<std::true_type, decltype(test_func<T>(nullptr))>
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
  static void fence(const ExecutionSpace& ex, const std::string& name,
                    ArgsDeduced&&... args) {
    if (Impl::parallel_reduce_needs_fence(ex, (ArgsDeduced&&)args...)) {
      ex.fence(name);
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
 *    void join(       <podType> & update ,
 *               const <podType> & input ) const ;
 *
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
 *    void join(       <podType> update[] ,
 *               const <podType> input[] ) const ;
 *
 *    void final( <podType> update[] ) const ;
 *  };
 * \endcode
 */

// ReturnValue is scalar or array: take by reference
template <class PolicyType, class FunctorType, class ReturnType>
inline std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value &&
                        !(Kokkos::is_view<ReturnType>::value ||
                          Kokkos::is_reducer<ReturnType>::value ||
                          std::is_pointer_v<ReturnType>)>
parallel_reduce(const std::string& label, const PolicyType& policy,
                const FunctorType& functor, ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check(
      "parallel_reduce", policy, label.c_str());

  static_assert(
      !std::is_const_v<ReturnType>,
      "A const reduction result type is only allowed for a View, pointer or "
      "reducer return type!");

  Impl::ParallelReduceAdaptor<PolicyType, FunctorType, ReturnType>::execute(
      label, policy, functor, return_value);
  Impl::ParallelReduceFence<typename PolicyType::execution_space, ReturnType>::
      fence(
          policy.space(),
          "Kokkos::parallel_reduce: fence due to result being value, not view",
          return_value);
}

template <class PolicyType, class FunctorType, class ReturnType>
inline std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value &&
                        !(Kokkos::is_view<ReturnType>::value ||
                          Kokkos::is_reducer<ReturnType>::value ||
                          std::is_pointer_v<ReturnType>)>
parallel_reduce(const PolicyType& policy, const FunctorType& functor,
                ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check("parallel_reduce",
                                                              policy);

  parallel_reduce("", policy, functor, return_value);
}

template <class FunctorType, class ReturnType>
inline std::enable_if_t<!(Kokkos::is_view<ReturnType>::value ||
                          Kokkos::is_reducer<ReturnType>::value ||
                          std::is_pointer_v<ReturnType>)>
parallel_reduce(const std::string& label, const size_t& work_count,
                const FunctorType& functor, ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check(
      "parallel_reduce", work_count, label.c_str());

  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  parallel_reduce(label, policy_type(0, work_count), functor, return_value);
}

template <class FunctorType, class ReturnType>
inline std::enable_if_t<!(Kokkos::is_view<ReturnType>::value ||
                          Kokkos::is_reducer<ReturnType>::value ||
                          std::is_pointer_v<ReturnType>)>
parallel_reduce(const size_t& work_count, const FunctorType& functor,
                ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check("parallel_reduce",
                                                              work_count);

  parallel_reduce("", work_count, functor, return_value);
}

// ReturnValue as View or Reducer: take by copy to allow for inline construction
template <class PolicyType, class FunctorType, class ReturnType>
inline std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value &&
                        (Kokkos::is_view<ReturnType>::value ||
                         Kokkos::is_reducer<ReturnType>::value ||
                         std::is_pointer_v<ReturnType>)>
parallel_reduce(const std::string& label, const PolicyType& policy,
                const FunctorType& functor, const ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check(
      "parallel_reduce", policy, label.c_str());

  ReturnType return_value_impl = return_value;
  Impl::ParallelReduceAdaptor<PolicyType, FunctorType, ReturnType>::execute(
      label, policy, functor, return_value_impl);
  Impl::ParallelReduceFence<typename PolicyType::execution_space, ReturnType>::
      fence(policy.space(),
            "Kokkos::parallel_reduce: fence" /*FIXME: describe correct reason*/,
            return_value);
}

template <class PolicyType, class FunctorType, class ReturnType>
inline std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value &&
                        (Kokkos::is_view<ReturnType>::value ||
                         Kokkos::is_reducer<ReturnType>::value ||
                         std::is_pointer_v<ReturnType>)>
parallel_reduce(const PolicyType& policy, const FunctorType& functor,
                const ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check("parallel_reduce",
                                                              policy);

  parallel_reduce("", policy, functor, return_value);
}

template <class FunctorType, class ReturnType>
inline std::enable_if_t<Kokkos::is_view<ReturnType>::value ||
                        Kokkos::is_reducer<ReturnType>::value ||
                        std::is_pointer_v<ReturnType>>
parallel_reduce(const std::string& label, const size_t& work_count,
                const FunctorType& functor, const ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check(
      "parallel_reduce", work_count, label.c_str());

  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;
  parallel_reduce(label, policy_type(0, work_count), functor, return_value);
}

template <class FunctorType, class ReturnType>
inline std::enable_if_t<Kokkos::is_view<ReturnType>::value ||
                        Kokkos::is_reducer<ReturnType>::value ||
                        std::is_pointer_v<ReturnType>>
parallel_reduce(const size_t& work_count, const FunctorType& functor,
                const ReturnType& return_value) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check("parallel_reduce",
                                                              work_count);

  parallel_reduce("", work_count, functor, return_value);
}

// No Return Argument
template <class PolicyType, class FunctorType>
inline void parallel_reduce(
    const std::string& label, const PolicyType& policy,
    const FunctorType& functor,
    std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value>* =
        nullptr) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check(
      "parallel_reduce", policy, label.c_str());

  using FunctorAnalysis =
      Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE, PolicyType,
                            FunctorType, void>;
  using value_type = std::conditional_t<(FunctorAnalysis::StaticValueSize != 0),
                                        typename FunctorAnalysis::value_type,
                                        typename FunctorAnalysis::pointer_type>;

  static_assert(
      FunctorAnalysis::has_final_member_function,
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
    std::enable_if_t<Kokkos::is_execution_policy<PolicyType>::value>* =
        nullptr) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check("parallel_reduce",
                                                              policy);

  parallel_reduce("", policy, functor);
}

template <class FunctorType>
inline void parallel_reduce(const std::string& label, const size_t& work_count,
                            const FunctorType& functor) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check(
      "parallel_reduce", work_count, label.c_str());

  using policy_type =
      typename Impl::ParallelReducePolicyType<void, size_t,
                                              FunctorType>::policy_type;

  parallel_reduce(label, policy_type(0, work_count), functor);
}

template <class FunctorType>
inline void parallel_reduce(const size_t& work_count,
                            const FunctorType& functor) {
  /** Enforce correct use **/
  Impl::CheckUsage<Impl::UsageRequires::insideExecEnv>::check("parallel_reduce",
                                                              work_count);

  parallel_reduce("", work_count, functor);
}

}  // namespace Kokkos

#endif  // KOKKOS_PARALLEL_REDUCE_HPP
