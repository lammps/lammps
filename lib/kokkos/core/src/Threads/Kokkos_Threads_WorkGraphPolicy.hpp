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

#ifndef KOKKOS_THREADS_WORKGRAPHPOLICY_HPP
#define KOKKOS_THREADS_WORKGRAPHPOLICY_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Threads/Kokkos_ThreadsExec.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::WorkGraphPolicy<Traits...>,
                  Kokkos::Threads> {
 private:
  using Policy = Kokkos::WorkGraphPolicy<Traits...>;

  using Self = ParallelFor<FunctorType, Kokkos::WorkGraphPolicy<Traits...>,
                           Kokkos::Threads>;

  Policy m_policy;
  FunctorType m_functor;

  template <class TagType>
  std::enable_if_t<std::is_void<TagType>::value> exec_one(
      const std::int32_t w) const noexcept {
    m_functor(w);
  }

  template <class TagType>
  std::enable_if_t<!std::is_void<TagType>::value> exec_one(
      const std::int32_t w) const noexcept {
    const TagType t{};
    m_functor(t, w);
  }

  inline void exec_one_thread() const noexcept {
    // Spin until COMPLETED_TOKEN.
    // END_TOKEN indicates no work is currently available.

    for (std::int32_t w = Policy::END_TOKEN;
         Policy::COMPLETED_TOKEN != (w = m_policy.pop_work());) {
      if (Policy::END_TOKEN != w) {
        exec_one<typename Policy::work_tag>(w);
        m_policy.completed_work(w);
      }
    }
  }

  static inline void thread_main(ThreadsExec& exec, const void* arg) noexcept {
    const Self& self = *(static_cast<const Self*>(arg));
    self.exec_one_thread();
    exec.fan_in();
  }

 public:
  inline void execute() {
    ThreadsExec::start(&Self::thread_main, this);
    ThreadsExec::fence();
  }

  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_policy(arg_policy), m_functor(arg_functor) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #define KOKKOS_THREADS_WORKGRAPHPOLICY_HPP */
