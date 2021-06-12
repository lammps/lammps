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

#ifndef KOKKOS_KOKKOS_HOST_GRAPHNODEKERNEL_HPP
#define KOKKOS_KOKKOS_HOST_GRAPHNODEKERNEL_HPP

#include <Kokkos_Macros.hpp>

#include <impl/Kokkos_Default_Graph_fwd.hpp>

#include <Kokkos_Graph.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Parallel_Reduce.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="GraphNodeKernelImpl"> {{{1

template <class ExecutionSpace>
struct GraphNodeKernelDefaultImpl {
  // TODO @graphs decide if this should use vtable or intrusive erasure via
  //      function pointers like in the rest of the graph interface
  virtual void execute_kernel() = 0;
};

// TODO Indicate that this kernel specialization is only for the Host somehow?
template <class ExecutionSpace, class PolicyType, class Functor,
          class PatternTag, class... Args>
class GraphNodeKernelImpl
    : public PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                              Args..., ExecutionSpace>::type,
      public GraphNodeKernelDefaultImpl<ExecutionSpace> {
 public:
  using base_t =
      typename PatternImplSpecializationFromTag<PatternTag, Functor, PolicyType,
                                                Args..., ExecutionSpace>::type;
  using execute_kernel_vtable_base_t =
      GraphNodeKernelDefaultImpl<ExecutionSpace>;
  // We have to use this name here because that's how it was done way back when
  // then implementations of Impl::Parallel*<> were written
  using Policy       = PolicyType;
  using graph_kernel = GraphNodeKernelImpl;

  // TODO @graph kernel name info propagation
  template <class PolicyDeduced, class... ArgsDeduced>
  GraphNodeKernelImpl(std::string const&, ExecutionSpace const&,
                      Functor arg_functor, PolicyDeduced&& arg_policy,
                      ArgsDeduced&&... args)
      : base_t(std::move(arg_functor), (PolicyDeduced &&) arg_policy,
               (ArgsDeduced &&) args...),
        execute_kernel_vtable_base_t() {}

  // FIXME @graph Forward through the instance once that works in the backends
  template <class PolicyDeduced, class... ArgsDeduced>
  GraphNodeKernelImpl(ExecutionSpace const& ex, Functor arg_functor,
                      PolicyDeduced&& arg_policy, ArgsDeduced&&... args)
      : GraphNodeKernelImpl("", ex, std::move(arg_functor),
                            (PolicyDeduced &&) arg_policy,
                            (ArgsDeduced &&) args...) {}

  void execute_kernel() final { this->base_t::execute(); }
};

// </editor-fold> end GraphNodeKernelImpl }}}1
//==============================================================================

template <class ExecutionSpace>
struct GraphNodeAggregateKernelDefaultImpl
    : GraphNodeKernelDefaultImpl<ExecutionSpace> {
  // Aggregates don't need a policy, but for the purposes of checking the static
  // assertions about graph kernels,
  struct Policy {
    using is_graph_kernel = std::true_type;
  };
  using graph_kernel = GraphNodeAggregateKernelDefaultImpl;
  void execute_kernel() final {}
};

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_HOST_GRAPHNODEKERNEL_HPP
