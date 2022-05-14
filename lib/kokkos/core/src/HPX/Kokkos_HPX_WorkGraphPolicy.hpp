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

#ifndef KOKKOS_HPX_WORKGRAPHPOLICY_HPP
#define KOKKOS_HPX_WORKGRAPHPOLICY_HPP

#include <Kokkos_HPX.hpp>

#include <hpx/local/algorithm.hpp>
#include <hpx/local/execution.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::WorkGraphPolicy<Traits...>,
                  Kokkos::Experimental::HPX> {
 private:
  using Policy  = Kokkos::WorkGraphPolicy<Traits...>;
  using WorkTag = typename Policy::work_tag;

  Policy m_policy;
  FunctorType m_functor;

  template <class TagType>
  typename std::enable_if<std::is_same<TagType, void>::value>::type
  execute_functor(const std::int32_t w) const noexcept {
    m_functor(w);
  }

  template <class TagType>
  typename std::enable_if<!std::is_same<TagType, void>::value>::type
  execute_functor(const std::int32_t w) const noexcept {
    const TagType t{};
    m_functor(t, w);
  }

 public:
  void execute() const {
    dispatch_execute_task(this, m_policy.space());
    m_policy.space().fence(
        "Kokkos::Experimental::Impl::HPX::ParallelFor<WorkGraphPolicy>: fence "
        "after kernel execution");
  }

  void execute_task() const {
    // See [note 1] in Kokkos_HPX.hpp for an explanation. The work graph policy
    // does not store an execution space instance, so we only need to reset the
    // parallel region count here.
    Kokkos::Experimental::HPX::reset_count_on_exit_parallel reset_count_on_exit;

    const int num_worker_threads = Kokkos::Experimental::HPX::concurrency();

    using hpx::for_loop;
    using hpx::execution::par;
    using hpx::execution::static_chunk_size;

    auto exec = Kokkos::Experimental::HPX::impl_get_executor();

    for_loop(par.on(exec).with(static_chunk_size(1)), 0, num_worker_threads,
             [this](int) {
               std::int32_t w = m_policy.pop_work();
               while (w != Policy::COMPLETED_TOKEN) {
                 if (w != Policy::END_TOKEN) {
                   execute_functor<WorkTag>(w);
                   m_policy.completed_work(w);
                 }

                 w = m_policy.pop_work();
               }
             });
  }

  inline ParallelFor(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_policy(arg_policy), m_functor(arg_functor) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #define KOKKOS_HPX_WORKGRAPHPOLICY_HPP */
