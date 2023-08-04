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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Core.hpp>

#include <Serial/Kokkos_Serial_Task.hpp>
#include <impl/Kokkos_TaskQueue_impl.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue<Kokkos::Serial, typename Kokkos::Serial::memory_space>;

}
}  // namespace Kokkos

#else
void KOKKOS_CORE_SRC_IMPL_SERIAL_TASK_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
