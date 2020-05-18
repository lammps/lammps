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

#ifndef KOKKOS_HPX_CHUNKEDROUNDROBINEXECUTOR_HPP
#define KOKKOS_HPX_CHUNKEDROUNDROBINEXECUTOR_HPP

#include <hpx/config.hpp>
#include <hpx/async_launch_policy_dispatch.hpp>
#include <hpx/lcos/local/latch.hpp>
#include <hpx/parallel/executors/execution.hpp>
#include <hpx/parallel/executors/post_policy_dispatch.hpp>
#include <hpx/runtime/get_os_thread_count.hpp>
#include <hpx/runtime/threads/thread_helpers.hpp>
#include <hpx/traits/is_executor.hpp>
#include <hpx/traits/is_launch_policy.hpp>
#include <hpx/util/deferred_call.hpp>

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace Kokkos {
namespace Impl {

///////////////////////////////////////////////////////////////////////////
/// A \a ChunkedRoundRobinExecutor creates groups of parallel execution
/// agents which execute in threads implicitly created by the executor. This
/// executor uses the scheduling hint to spawn threads with the first grouped on
/// the first core, the second group getting the next consecutive threads, etc.
/// For example, if 10 tasks are spawned (num_tasks is set to 10) and num_cores
/// is set to 2 the executor will schedule the tasks in the following order:
///
/// worker thread | 1 | 2
/// --------------+---+---
/// tasks         | 1 | 6
///               | 2 | 7
///               | 3 | 8
///               | 4 | 9
///               | 5 | 10
///
/// rather than the typical round robin:
///
/// worker thread | 1 | 2
/// --------------+---+---
/// tasks         | 1 | 2
///               | 3 | 4
///               | 5 | 6
///               | 7 | 8
///               | 9 | 10
struct ChunkedRoundRobinExecutor {
  using execution_category = hpx::parallel::execution::parallel_execution_tag;

  HPX_CONSTEXPR explicit ChunkedRoundRobinExecutor(
      std::size_t num_tasks = std::size_t(-1), std::size_t core_offset = 0,
      std::size_t num_cores = hpx::get_os_thread_count())
      : num_tasks_(num_tasks),
        core_offset_(core_offset),
        num_cores_(num_cores),
        num_tasks_per_core_(double(num_tasks_) / num_cores_),
        num_tasks_spawned_(0) {}

  bool operator==(ChunkedRoundRobinExecutor const &rhs) const noexcept {
    return num_cores_ == rhs.num_cores_ && num_tasks_ == rhs.num_tasks_;
  }

  bool operator!=(ChunkedRoundRobinExecutor const &rhs) const noexcept {
    return !(*this == rhs);
  }

  ChunkedRoundRobinExecutor const &context() const noexcept { return *this; }

  template <typename F, typename... Ts>
  hpx::future<
      typename hpx::util::detail::invoke_deferred_result<F, Ts...>::type>
  async_execute(F &&f, Ts &&... ts) const {
    return hpx::detail::async_launch_policy_dispatch<hpx::launch>::call(
        hpx::launch::async_policy{}, std::forward<F>(f),
        std::forward<Ts>(ts)...);
  }

  template <typename F, typename... Ts>
  void post(F &&f, Ts &&... ts) const {
    hpx::util::thread_description const desc(
        f, "Kokkos::Impl::ChunkedRoundRobinExecutor::async_execute");
    hpx::threads::thread_schedule_hint const hint(
        hpx::threads::thread_schedule_hint_mode_thread,
        core_offset_ + std::floor(double(num_tasks_spawned_ % num_tasks_) /
                                  num_tasks_per_core_));

    hpx::threads::register_thread_nullary(
        hpx::util::deferred_call(std::forward<F>(f), std::forward<Ts>(ts)...),
        desc, hpx::threads::pending, false,
        hpx::threads::thread_priority_normal, hint,
        hpx::threads::thread_stacksize_default);

    ++num_tasks_spawned_;
  }

  template <typename F, typename Shape, typename... Ts>
  std::vector<hpx::future<typename hpx::parallel::execution::detail::
                              bulk_function_result<F, Shape, Ts...>::type>>
  bulk_async_execute(F &&f, Shape const &shape, Ts &&... ts) {
    hpx::util::thread_description desc(
        f, "Kokkos::Impl::ChunkedRoundRobinExecutor::bulk_sync_execute");

    hpx::lcos::local::latch l(hpx::util::size(shape));
    // Keep a separate counter for bulk launch
    std::size_t num_tasks_spawned = 0;

    for (auto const &s : shape) {
      hpx::threads::thread_schedule_hint const hint(
          hpx::threads::thread_schedule_hint_mode_thread,
          core_offset_ + std::floor(double(num_tasks_spawned % num_tasks_) /
                                    num_tasks_per_core_));

      hpx::threads::register_thread_nullary(
          [&, s]() {
            hpx::util::invoke(f, s, ts...);
            l.count_down(1);
          },
          desc, hpx::threads::pending, false,
          hpx::threads::thread_priority_normal, hint,
          hpx::threads::thread_stacksize_default);

      ++num_tasks_spawned;
    }

    // NOTE: We block here to avoid extra synchronization. Since this executor
    // is only used in the HPX backend we get away with this.
    l.wait();

    return {};
  }

 private:
  std::size_t num_tasks_;
  std::size_t core_offset_;
  std::size_t num_cores_;
  double num_tasks_per_core_;
  mutable std::size_t num_tasks_spawned_;
};

}  // namespace Impl
}  // namespace Kokkos

namespace hpx {
namespace parallel {
namespace execution {

template <>
struct is_one_way_executor<Kokkos::Impl::ChunkedRoundRobinExecutor>
    : std::true_type {};

template <>
struct is_two_way_executor<Kokkos::Impl::ChunkedRoundRobinExecutor>
    : std::true_type {};

template <>
struct is_bulk_two_way_executor<Kokkos::Impl::ChunkedRoundRobinExecutor>
    : std::true_type {};

}  // namespace execution
}  // namespace parallel
}  // namespace hpx

#endif
