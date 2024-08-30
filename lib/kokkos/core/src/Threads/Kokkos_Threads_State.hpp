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
