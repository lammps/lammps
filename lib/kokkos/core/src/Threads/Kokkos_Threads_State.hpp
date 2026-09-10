// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_THREADS_STATE_HPP
#define KOKKOS_THREADS_STATE_HPP

namespace Kokkos {
namespace Impl {
/** \brief States of a worker thread */
enum class ThreadState {
  Terminating  ///<  Termination in progress
  ,
  Inactive  ///<  Exists, waiting for work
  ,
  Active  ///<  Exists, performing work
  ,
  Rendezvous  ///<  Exists, waiting in a barrier or reduce
  ,
  ScanCompleted,
  ScanAvailable,
  ReductionAvailable
};
}  // namespace Impl
}  // namespace Kokkos

#endif
