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

#include <HPX/Kokkos_HPX.hpp>

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

 public:
  void execute_range(int) const {
    std::int32_t w = m_policy.pop_work();
    while (w != Policy::COMPLETED_TOKEN) {
      if (w != Policy::END_TOKEN) {
        if constexpr (std::is_same_v<WorkTag, void>) {
          m_functor(w);
        } else {
          m_functor(WorkTag{}, w);
        }
        m_policy.completed_work(w);
      }

      w = m_policy.pop_work();
    }
  }

  void execute() const {
    const int num_worker_threads = Kokkos::Experimental::HPX().concurrency();
    Kokkos::Experimental::HPX().impl_bulk_plain(
        true, is_light_weight_policy<Policy>(), *this, num_worker_threads,
        hpx::threads::thread_stacksize::nostack);
  }

  inline ParallelFor(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_policy(arg_policy), m_functor(arg_functor) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #define KOKKOS_HPX_WORKGRAPHPOLICY_HPP */
