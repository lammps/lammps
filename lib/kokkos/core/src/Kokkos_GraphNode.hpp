// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_KOKKOS_GRAPHNODE_HPP
#define KOKKOS_KOKKOS_GRAPHNODE_HPP

#include <Kokkos_Macros.hpp>

#include <impl/Kokkos_Error.hpp>  // contract macros

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Graph_fwd.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>
#include <Kokkos_Parallel_Reduce.hpp>
#include <impl/Kokkos_GraphImpl_Utilities.hpp>
#include <impl/Kokkos_GraphImpl.hpp>  // GraphAccess
#include <impl/Kokkos_GraphNodeThenPolicy.hpp>
#include <impl/Kokkos_GraphNodeCtorProps.hpp>

#include <memory>  // std::shared_ptr

namespace Kokkos {
namespace Experimental {

#define KOKKOS_IMPL_POLICY_ON_DEFAULT_EXEC(policy)                         \
  if (policy.space() != std::remove_cvref_t<decltype(policy.space())>{})   \
    Kokkos::abort(                                                         \
        "The execution space instance of the execution policy of a graph " \
        "node must be the default one.");

template <typename T, typename Exec>
concept ExecutionPolicyOn =
    ExecutionPolicy<T> && std::same_as<typename T::execution_space, Exec>;

template <Kokkos::ExecutionSpace ExecutionSpace,
          class Kernel /*= TypeErasedTag*/,
          class Predecessor /*= TypeErasedTag*/>
class GraphNodeRef {
  //----------------------------------------------------------------------------
  // <editor-fold desc="template parameter constraints"> {{{2

  // Note: because of these assertions, instantiating this class template is not
  //       intended to be SFINAE-safe, so do validation before you instantiate.

  static_assert(
      std::is_same_v<Predecessor, TypeErasedTag> ||
          Kokkos::Impl::is_specialization_of<Predecessor, GraphNodeRef>::value,
      "Invalid predecessor template parameter given to GraphNodeRef");

  static_assert(std::is_same_v<Predecessor, TypeErasedTag> ||
                    Kokkos::Impl::is_graph_kernel<Kernel>::value ||
                    Kokkos::Impl::is_graph_capture_v<Kernel> ||
                    Kokkos::Impl::is_graph_then_host_v<Kernel>,
                "Invalid kernel template parameter given to GraphNodeRef");

  static_assert(!Kokkos::Impl::is_more_type_erased<Kernel, Predecessor>::value,
                "The kernel of a graph node can't be more type-erased than the "
                "predecessor");

  // </editor-fold> end template parameter constraints }}}2
  //----------------------------------------------------------------------------

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="public member types"> {{{2

  using execution_space   = ExecutionSpace;
  using graph_kernel      = Kernel;
  using graph_predecessor = Predecessor;

  // </editor-fold> end public member types }}}2
  //----------------------------------------------------------------------------

 private:
  using device_handle_t = Kokkos::Impl::DeviceHandle<ExecutionSpace>;

  //----------------------------------------------------------------------------
  // <editor-fold desc="Friends"> {{{2

  template <Kokkos::ExecutionSpace, class, class>
  friend class GraphNodeRef;
  friend struct Kokkos::Impl::GraphAccess;
  friend struct Graph<execution_space>;

  // </editor-fold> end Friends }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="Private Data Members"> {{{2

  using graph_impl_t = Kokkos::Impl::GraphImpl<ExecutionSpace>;
  std::weak_ptr<graph_impl_t> m_graph_impl;

  // TODO @graphs figure out if we can get away with a weak reference here?
  //              GraphNodeRef instances shouldn't be stored by users outside
  //              of the create_graph closure, and so if we restructure things
  //              slightly, we could make it so that the graph owns the
  //              node_impl_t instance and this only holds a std::weak_ptr to
  //              it.
  using node_impl_t =
      Kokkos::Impl::GraphNodeImpl<ExecutionSpace, Kernel, Predecessor>;
  std::shared_ptr<node_impl_t> m_node_impl;

  // </editor-fold> end Private Data Members }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="Implementation detail accessors"> {{{2

  // Internally, use shallow constness
  node_impl_t& get_node_impl() const { return *m_node_impl.get(); }
  std::shared_ptr<node_impl_t> const& get_node_ptr() const& {
    return m_node_impl;
  }
  std::shared_ptr<node_impl_t> get_node_ptr() && {
    return std::move(m_node_impl);
  }
  std::weak_ptr<graph_impl_t> get_graph_weak_ptr() const {
    return m_graph_impl;
  }

  // </editor-fold> end Implementation detail accessors }}}2
  //----------------------------------------------------------------------------

  // TODO kernel name propagation and exposure

  template <class NextKernelDeduced>
  auto _then_kernel(NextKernelDeduced&& arg_kernel) const {
    static_assert(
        Kokkos::Impl::is_graph_kernel_v<std::remove_cvref_t<NextKernelDeduced>>,
        "Kokkos internal error");

    auto graph_ptr = m_graph_impl.lock();
    KOKKOS_EXPECTS(bool(graph_ptr))

    using next_kernel_t = std::remove_cvref_t<NextKernelDeduced>;

    using return_t = GraphNodeRef<ExecutionSpace, next_kernel_t, GraphNodeRef>;

    auto rv = Kokkos::Impl::GraphAccess::make_graph_node_ref(
        m_graph_impl,
        Kokkos::Impl::GraphAccess::make_node_shared_ptr<
            typename return_t::node_impl_t>(
            m_node_impl->get_device_handle(),
            Kokkos::Impl::_graph_node_kernel_ctor_tag{},
            (NextKernelDeduced&&)arg_kernel,
            // *this is the predecessor
            Kokkos::Impl::_graph_node_predecessor_ctor_tag{}, *this));

    // Add the node itself to the backend's graph data structure, now that
    // everything is set up.
    graph_ptr->add_node(rv.m_node_impl);
    // Add the predecessor we stored in the constructor above in the backend's
    // data structure, now that everything is set up.
    graph_ptr->add_predecessor(rv.m_node_impl, *this);
    KOKKOS_ENSURES(bool(rv.m_node_impl))
    return rv;
  }

  //----------------------------------------------------------------------------
  // <editor-fold desc="Private constructors"> {{{2

  GraphNodeRef(std::weak_ptr<graph_impl_t> arg_graph_impl,
               std::shared_ptr<node_impl_t> arg_node_impl)
      : m_graph_impl(std::move(arg_graph_impl)),
        m_node_impl(std::move(arg_node_impl)) {}

  // </editor-fold> end Private constructors }}}2
  //----------------------------------------------------------------------------

 public:
  //----------------------------------------------------------------------------
  // <editor-fold desc="Constructors, destructors, and assignment"> {{{2

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // <editor-fold desc="rule of 6 ctors"> {{{3

  // Copyable and movable (basically just shared_ptr semantics
  GraphNodeRef() noexcept = default;

  // </editor-fold> end rule of 6 ctors }}}3
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // <editor-fold desc="Type-erasing converting ctor and assignment"> {{{3

  template <class OtherKernel, class OtherPredecessor,
            std::enable_if_t<
                // Not a copy/move constructor
                !std::is_same_v<GraphNodeRef,
                                GraphNodeRef<execution_space, OtherKernel,
                                             OtherPredecessor>> &&
                    // must be an allowed type erasure of the kernel
                    Kokkos::Impl::is_compatible_type_erasure<
                        OtherKernel, graph_kernel>::value &&
                    // must be an allowed type erasure of the predecessor
                    Kokkos::Impl::is_compatible_type_erasure<
                        OtherPredecessor, graph_predecessor>::value,
                int> = 0>
  /* implicit */
  GraphNodeRef(
      GraphNodeRef<execution_space, OtherKernel, OtherPredecessor> const& other)
      : m_graph_impl(other.m_graph_impl), m_node_impl(other.m_node_impl) {}

  // Note: because this is an implicit conversion (as is supposed to be the
  //       case with most type-erasing wrappers like this), we don't also need
  //       a converting assignment operator.

  // </editor-fold> end Type-erasing converting ctor and assignment }}}3
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // </editor-fold> end Constructors, destructors, and assignment }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="then_parallel_for"> {{{2

  // TODO We should do better than a p-for (that uses registers, heavier).
  //      This should "just" launch the function on device with our driver.
  template <typename Props, typename Policy, typename Functor>
    requires(Kokkos::Impl::NodeProperties<std::remove_cvref_t<Props>> &&
             std::is_invocable_r_v<void, const std::remove_cvref_t<Functor>> &&
             Kokkos::Impl::is_specialization_of_v<Policy, ThenPolicy>)
  auto then(Props&& props, Policy&& policy, Functor&& functor) const {
    using next_kernel_t =
        Kokkos::Impl::GraphNodeThenImpl<ExecutionSpace,
                                        std::remove_cvref_t<Policy>,
                                        std::remove_cvref_t<Functor>>;
    auto graph_ptr = m_graph_impl.lock();
    KOKKOS_EXPECTS(bool(graph_ptr))
    auto full_props = Kokkos::Impl::with_properties_if_unset(
        std::forward<Props>(props), graph_ptr->get_device_handle(),
        "[unlabeled]");
    return this->_then_kernel(next_kernel_t{
        Kokkos::Impl::extract_property<std::string>(full_props),
        Kokkos::Impl::extract_property<device_handle_t>(full_props).m_exec,
        std::forward<Policy>(policy), std::forward<Functor>(functor)});
  }

  template <typename Props, typename Functor>
    requires Kokkos::Impl::NodeProperties<std::remove_cvref_t<Props>>
  auto then(Props&& props, Functor&& functor) const {
    return this->then(std::forward<Props>(props), ThenPolicy{},
                      std::forward<Functor>(functor));
  }

  // Overload for a label.
  template <typename Label, typename Functor>
    requires Kokkos::Impl::ViewLabel<std::remove_cvref_t<Label>>
  auto then(Label&& label, Functor&& functor) const {
    return this->then(node_props(std::forward<Label>(label)), ThenPolicy{},
                      std::forward<Functor>(functor));
  }

  // Overload for a policy.
  template <typename Policy, typename Functor,
            std::enable_if_t<Kokkos::Impl::is_specialization_of_v<
                                 std::remove_cvref_t<Policy>, ThenPolicy>,
                             int> = 0>
  auto then(Policy&& policy, Functor&& functor) const {
    return this->then(node_props(), std::forward<Policy>(policy),
                      std::forward<Functor>(functor));
  }

  // Overload for a label and a policy.
  template <typename Label, typename Policy, typename Functor,
            std::enable_if_t<
                Kokkos::Impl::is_specialization_of_v<
                    std::remove_cvref_t<Policy>, ThenPolicy> &&
                    Kokkos::Impl::is_view_label_v<std::remove_cvref_t<Label>>,
                int> = 0>
  auto then(Label&& label, Policy&& policy, Functor&& functor) const {
    return this->then(node_props(std::forward<Label>(label)),
                      std::forward<Policy>(policy),
                      std::forward<Functor>(functor));
  }

  // Overload without policy given.
  template <typename Functor>
  auto then(Functor&& functor) const {
    return this->then(node_props(), ThenPolicy{},
                      std::forward<Functor>(functor));
  }

  template <typename Label, typename Functor>
  auto then_host(Label&&, Functor&& functor) const {
    using host_t =
        Kokkos::Impl::GraphNodeThenHostImpl<ExecutionSpace,
                                            std::remove_cvref_t<Functor>>;
    using return_t = GraphNodeRef<ExecutionSpace, host_t, GraphNodeRef>;

    auto graph_ptr = m_graph_impl.lock();
    KOKKOS_EXPECTS(bool(graph_ptr))

    auto rv = Kokkos::Impl::GraphAccess::make_graph_node_ref(
        m_graph_impl,
        Kokkos::Impl::GraphAccess::make_node_shared_ptr<
            typename return_t::node_impl_t>(
            m_node_impl->get_device_handle(),
            Kokkos::Impl::_graph_node_host_ctor_tag{},
            std::forward<Functor>(functor),
            Kokkos::Impl::_graph_node_predecessor_ctor_tag{}, *this));

    // Add the node itself to the backend's graph data structure, now that
    // everything is set up.
    graph_ptr->add_node(rv.m_node_impl);
    // Add the predecessor we stored in the constructor above in the
    // backend's data structure, now that everything is set up.
    graph_ptr->add_predecessor(rv.m_node_impl, *this);
    KOKKOS_ENSURES(bool(rv.m_node_impl))
    return rv;
  }

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    (defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT))
  template <
      class Functor,
      typename = std::enable_if_t<std::is_invocable_r_v<
          void, const std::remove_cvref_t<Functor>, const ExecutionSpace&>>>
#if defined(KOKKOS_ENABLE_CUDA)
  auto cuda_capture(const ExecutionSpace& exec, Functor&& functor) const {
    if constexpr (std::is_same_v<ExecutionSpace, Kokkos::Cuda>) {
#elif defined(KOKKOS_ENABLE_HIP)
  auto hip_capture(const ExecutionSpace& exec, Functor&& functor) const {
    if constexpr (std::is_same_v<ExecutionSpace, Kokkos::HIP>) {
#elif defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_SYCL_GRAPH_SUPPORT)
  auto sycl_capture(const ExecutionSpace& exec, Functor&& functor) const {
    if constexpr (std::is_same_v<ExecutionSpace, Kokkos::SYCL>) {
#endif
      using capture_t =
          Kokkos::Impl::GraphNodeCaptureImpl<ExecutionSpace,
                                             std::remove_cvref_t<Functor>>;
      using return_t = GraphNodeRef<ExecutionSpace, capture_t, GraphNodeRef>;

      auto graph_ptr = m_graph_impl.lock();
      KOKKOS_EXPECTS(bool(graph_ptr))

      auto rv = Kokkos::Impl::GraphAccess::make_graph_node_ref(
          m_graph_impl,
          Kokkos::Impl::GraphAccess::make_node_shared_ptr<
              typename return_t::node_impl_t>(
              m_node_impl->get_device_handle(),
              Kokkos::Impl::_graph_node_capture_ctor_tag{},
              std::forward<Functor>(functor),
              Kokkos::Impl::_graph_node_predecessor_ctor_tag{}, *this));

      // Add the node itself to the backend's graph data structure, now that
      // everything is set up.
      graph_ptr->add_node(exec, rv.m_node_impl);
      // Add the predecessor we stored in the constructor above in the
      // backend's data structure, now that everything is set up.
      graph_ptr->add_predecessor(rv.m_node_impl, *this);
      KOKKOS_ENSURES(bool(rv.m_node_impl))
      return rv;
    }
  }
#endif

  template <typename Props, class Policy, class Functor>
    requires(Kokkos::Impl::NodeProperties<std::remove_cvref_t<Props>> &&
             ExecutionPolicyOn<std::remove_cvref_t<Policy>, ExecutionSpace>)
  auto then_parallel_for(Props&& props, Policy&& arg_policy,
                         Functor&& functor) const {
    //----------------------------------------
    KOKKOS_EXPECTS(!m_graph_impl.expired())
    KOKKOS_EXPECTS(bool(m_node_impl))

    KOKKOS_IMPL_POLICY_ON_DEFAULT_EXEC(arg_policy)

    // needs to static assert constraint: DataParallelFunctor<Functor>

    auto graph_ptr = m_graph_impl.lock();
    KOKKOS_EXPECTS(bool(graph_ptr))
    auto full_props =
        with_properties_if_unset(std::forward<Props>(props),
                                 graph_ptr->get_device_handle(), "[unlabeled]");

    auto policy = Experimental::require(
        Policy(Kokkos::Impl::PolicyUpdate{}, (Policy&&)arg_policy,
               Kokkos::Impl::get_property<device_handle_t>(full_props).m_exec),
        Kokkos::Impl::KernelInGraphProperty{});

    using next_policy_t = decltype(policy);
    using next_kernel_t =
        Kokkos::Impl::GraphNodeKernelImpl<ExecutionSpace, next_policy_t,
                                          std::decay_t<Functor>,
                                          Kokkos::ParallelForTag>;
    return this->_then_kernel(next_kernel_t{
        Kokkos::Impl::extract_property<std::string>(full_props),
        Kokkos::Impl::extract_property<device_handle_t>(full_props).m_exec,
        (Functor&&)functor, std::move(policy)});
  }

  template <class Policy, class Functor>
    requires ExecutionPolicyOn<std::remove_cvref_t<Policy>, ExecutionSpace>
  auto then_parallel_for(Policy&& policy, Functor&& functor) const {
    // needs to static assert constraint: DataParallelFunctor<Functor>
    return this->then_parallel_for(node_props(), (Policy&&)policy,
                                   (Functor&&)functor);
  }

  template <class Label, class Policy, class Functor>
    requires(ExecutionPolicyOn<std::remove_cvref_t<Policy>, ExecutionSpace> &&
             Kokkos::Impl::is_view_label_v<std::remove_cvref_t<Label>>)
  auto then_parallel_for(Label&& label, Policy&& policy,
                         Functor&& functor) const {
    // needs to static assert constraint: DataParallelFunctor<Functor>
    return this->then_parallel_for(node_props(std::forward<Label>(label)),
                                   (Policy&&)policy, (Functor&&)functor);
  }

  template <class Functor>
  auto then_parallel_for(std::string name, std::size_t n,
                         Functor&& functor) const {
    // needs to static assert constraint: DataParallelFunctor<Functor>
    return this->then_parallel_for(node_props(std::move(name)),
                                   Kokkos::RangePolicy<execution_space>(0, n),
                                   (Functor&&)functor);
  }

  template <class Functor>
  auto then_parallel_for(std::size_t n, Functor&& functor) const {
    // needs to static assert constraint: DataParallelFunctor<Functor>
    return this->then_parallel_for(node_props(),
                                   Kokkos::RangePolicy<execution_space>(0, n),
                                   (Functor&&)functor);
  }

  // </editor-fold> end then_parallel_for }}}2
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // <editor-fold desc="then_parallel_reduce"> {{{2

  // Equivalent to std::get<I>(std::tuple) but callable on the device.
  template <bool B, class T1, class T2>
  static KOKKOS_FUNCTION std::conditional_t<B, T1&&, T2&&>
  impl_forwarding_switch(T1&& v1, T2&& v2) {
    if constexpr (B)
      return static_cast<T1&&>(v1);
    else
      return static_cast<T2&&>(v2);
  }

  template <typename Props, class Policy, class Functor, class ReturnType>
    requires(Kokkos::Impl::NodeProperties<std::remove_cvref_t<Props>> &&
             ExecutionPolicyOn<std::remove_cvref_t<Policy>, ExecutionSpace>)
  auto then_parallel_reduce(Props&& props, Policy&& arg_policy,
                            Functor&& functor,
                            ReturnType&& return_value) const {
    auto graph_impl_ptr = m_graph_impl.lock();
    KOKKOS_EXPECTS(bool(graph_impl_ptr))
    KOKKOS_EXPECTS(bool(m_node_impl))

    KOKKOS_IMPL_POLICY_ON_DEFAULT_EXEC(arg_policy)

    // needs static assertion of constraint:
    //   DataParallelReductionFunctor<Functor, ReturnType>

    auto full_props = with_properties_if_unset(
        std::forward<Props>(props), graph_impl_ptr->get_device_handle(),
        "[unlabeled]");

    // This is also just an expectation, but it's one that we expect the user
    // to interact with (even in release mode), so we should throw an exception
    // with an explanation rather than just doing a contract assertion.
    // We can't static_assert this because of the way that Reducers store
    // whether or not they point to a View as a runtime boolean rather than part
    // of the type.
    if (Kokkos::Impl::parallel_reduce_needs_fence(
            Kokkos::Impl::get_property<device_handle_t>(full_props).m_exec,
            return_value)) {
      Kokkos::Impl::throw_runtime_exception(
          "Parallel reductions in graphs can't operate on Reducers that "
          "reference a scalar because they can't complete synchronously. Use a "
          "Kokkos::View instead and keep in mind the result will only be "
          "available once the graph is submitted (or in tasks that depend on "
          "this one).");
    }

    //----------------------------------------
    // This is a disaster, but I guess it's not a my disaster to fix right now
    using return_type_remove_cvref =
        std::remove_cv_t<std::remove_reference_t<ReturnType>>;
    static_assert(Kokkos::is_view<return_type_remove_cvref>::value ||
                      Kokkos::is_reducer<return_type_remove_cvref>::value,
                  "Output argument to parallel reduce in a graph must be a "
                  "View or a Reducer");

    if constexpr (Kokkos::is_reducer_v<return_type_remove_cvref>) {
      static_assert(
          Kokkos::SpaceAccessibility<
              ExecutionSpace, typename return_type_remove_cvref::
                                  result_view_type::memory_space>::accessible,
          "The reduction target must be accessible by the graph execution "
          "space.");
    } else {
      static_assert(
          Kokkos::SpaceAccessibility<
              ExecutionSpace,
              typename return_type_remove_cvref::memory_space>::accessible,
          "The reduction target must be accessible by the graph execution "
          "space.");
    }

    using return_type =
        // Yes, you do really have to do this...
        std::conditional_t<Kokkos::is_reducer<return_type_remove_cvref>::value,
                           return_type_remove_cvref,
                           const return_type_remove_cvref>;
    using functor_type = std::remove_cvref_t<Functor>;
    // see Kokkos_Parallel_Reduce.hpp for how these details are used there;
    // we're just doing the same thing here
    using return_value_adapter =
        Kokkos::Impl::ParallelReduceReturnValue<void, return_type,
                                                functor_type>;
    // End of Kokkos reducer disaster
    //----------------------------------------

    auto policy = Experimental::require(
        Policy(Kokkos::Impl::PolicyUpdate{}, (Policy&&)arg_policy,
               Kokkos::Impl::get_property<device_handle_t>(full_props).m_exec),
        Kokkos::Impl::KernelInGraphProperty{});

    using passed_reducer_type = typename return_value_adapter::reducer_type;

    constexpr bool passed_reducer_type_is_invalid =
        std::is_same_v<InvalidType, passed_reducer_type>;
    using TheReducerType =
        std::conditional_t<passed_reducer_type_is_invalid, functor_type,
                           passed_reducer_type>;

    using analysis = Kokkos::Impl::FunctorAnalysis<
        Kokkos::Impl::FunctorPatternInterface::REDUCE, Policy, TheReducerType,
        typename return_value_adapter::value_type>;
    typename analysis::Reducer final_reducer(
        impl_forwarding_switch<passed_reducer_type_is_invalid>(functor,
                                                               return_value));
    Kokkos::Impl::CombinedFunctorReducer<functor_type,
                                         typename analysis::Reducer>
        functor_reducer(functor, final_reducer);

    using next_policy_t = decltype(policy);
    using next_kernel_t =
        Kokkos::Impl::GraphNodeKernelImpl<ExecutionSpace, next_policy_t,
                                          decltype(functor_reducer),
                                          Kokkos::ParallelReduceTag>;

    return this->_then_kernel(next_kernel_t{
        Kokkos::Impl::extract_property<std::string>(full_props),
        Kokkos::Impl::extract_property<device_handle_t>(full_props).m_exec,
        std::move(functor_reducer), std::move(policy),
        return_value_adapter::return_value(return_value,
                                           std::forward<Functor>(functor))});
  }

  template <class Policy, class Functor, class ReturnType>
    requires ExecutionPolicyOn<std::remove_cvref_t<Policy>, ExecutionSpace>
  auto then_parallel_reduce(Policy&& arg_policy, Functor&& functor,
                            ReturnType&& return_value) const {
    return this->then_parallel_reduce(node_props(), (Policy&&)arg_policy,
                                      (Functor&&)functor,
                                      (ReturnType&&)return_value);
  }

  template <class Functor, class ReturnType>
  auto then_parallel_reduce(std::string label,
                            typename execution_space::size_type idx_end,
                            Functor&& functor,
                            ReturnType&& return_value) const {
    return this->then_parallel_reduce(
        node_props(std::move(label)),
        Kokkos::RangePolicy<execution_space>{0, idx_end}, (Functor&&)functor,
        (ReturnType&&)return_value);
  }

  template <class Functor, class ReturnType>
  auto then_parallel_reduce(typename execution_space::size_type idx_end,
                            Functor&& functor,
                            ReturnType&& return_value) const {
    return this->then_parallel_reduce(
        node_props(), Kokkos::RangePolicy<execution_space>{0, idx_end},
        (Functor&&)functor, (ReturnType&&)return_value);
  }

  // </editor-fold> end then_parallel_reduce }}}2
  //----------------------------------------------------------------------------

  // TODO @graph parallel scan, deep copy, etc.
};

}  // end namespace Experimental
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_GRAPHNODE_HPP
