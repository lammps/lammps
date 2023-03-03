//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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
  std::enable_if_t<std::is_void<TagType>::value> execute_functor(
      const std::int32_t w) const noexcept {
    m_functor(w);
  }

  template <class TagType>
  std::enable_if_t<!std::is_void<TagType>::value> execute_functor(
      const std::int32_t w) const noexcept {
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
