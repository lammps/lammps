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
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_TaskQueue_impl.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue<
    Kokkos::Cuda,
    Impl::default_tasking_memory_space_for_execution_space_t<Kokkos::Cuda> >;
template class TaskQueueMultiple<
    Kokkos::Cuda,
    Impl::default_tasking_memory_space_for_execution_space_t<Kokkos::Cuda> >;

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
#else
void KOKKOS_CORE_SRC_CUDA_KOKKOS_CUDA_TASK_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_CUDA ) && defined( KOKKOS_ENABLE_TASKDAG \
          ) */
