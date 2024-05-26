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

#ifndef KOKKOS_THREADS_PARALLEL_FOR_MDRANGE_HPP
#define KOKKOS_THREADS_PARALLEL_FOR_MDRANGE_HPP

#include <Kokkos_Parallel.hpp>

#include <KokkosExp_MDRangePolicy.hpp>
namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Threads> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;

  using WorkTag = typename MDRangePolicy::work_tag;

  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  using iterate_type = typename Kokkos::Impl::HostIterateTile<
      MDRangePolicy, FunctorType, typename MDRangePolicy::work_tag, void>;

  const iterate_type m_iter;

  inline void exec_range(const Member ibeg, const Member iend) const {
    for (Member i = ibeg; i < iend; ++i) {
      m_iter(i);
    }
  }

  static void exec(ThreadsInternal &instance, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(instance, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsInternal &instance, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    auto const num_tiles = self.m_iter.m_rp.m_num_tiles;
    WorkRange range(Policy(0, num_tiles).set_chunk_size(1),
                    instance.pool_rank(), instance.pool_size());

    self.exec_range(range.begin(), range.end());

    instance.fan_in();
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsInternal &instance, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    auto const num_tiles = self.m_iter.m_rp.m_num_tiles;
    WorkRange range(Policy(0, num_tiles).set_chunk_size(1),
                    instance.pool_rank(), instance.pool_size());

    instance.set_work_range(range.begin(), range.end(), 1);
    instance.reset_steal_target();
    instance.barrier();

    long work_index = instance.get_work_index();

    while (work_index != -1) {
      const Member begin = static_cast<Member>(work_index);
      const Member end   = begin + 1 < num_tiles ? begin + 1 : num_tiles;

      self.exec_range(begin, end);
      work_index = instance.get_work_index();
    }

    instance.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsInternal::start(&ParallelFor::exec, this);
    ThreadsInternal::fence();
  }

  ParallelFor(const FunctorType &arg_functor, const MDRangePolicy &arg_policy)
      : m_iter(arg_policy, arg_functor) {}

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy &, const Functor &) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
};

}  // namespace Impl
}  // namespace Kokkos
#endif
