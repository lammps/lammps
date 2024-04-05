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

#ifndef KOKKOS_HIP_PARALLEL_FOR_MDRANGE_HPP
#define KOKKOS_HIP_PARALLEL_FOR_MDRANGE_HPP

#include <Kokkos_Parallel.hpp>

#include <HIP/Kokkos_HIP_BlockSize_Deduction.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <KokkosExp_MDRangePolicy.hpp>
#include <impl/KokkosExp_IterateTileGPU.hpp>

namespace Kokkos {
namespace Impl {

// ParallelFor
template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>, HIP> {
 public:
  using Policy       = Kokkos::MDRangePolicy<Traits...>;
  using functor_type = FunctorType;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;
  using LaunchBounds     = typename Policy::launch_bounds;

  const FunctorType m_functor;
  const Policy m_policy;

 public:
  ParallelFor()                   = delete;
  ParallelFor(ParallelFor const&) = default;
  ParallelFor& operator=(ParallelFor const&) = delete;

  inline __device__ void operator()() const {
    Kokkos::Impl::DeviceIterateTile<Policy::rank, Policy, FunctorType,
                                    typename Policy::work_tag>(m_policy,
                                                               m_functor)
        .exec_range();
  }

  inline void execute() const {
    using ClosureType = ParallelFor<FunctorType, Policy, HIP>;
    if (m_policy.m_num_tiles == 0) return;
    auto const maxblocks = hip_internal_maximum_grid_count();
    if (Policy::rank == 2) {
      dim3 const block(m_policy.m_tile[0], m_policy.m_tile[1], 1);
      dim3 const grid(
          std::min<array_index_type>(
              (m_policy.m_upper[0] - m_policy.m_lower[0] + block.x - 1) /
                  block.x,
              maxblocks[0]),
          std::min<array_index_type>(
              (m_policy.m_upper[1] - m_policy.m_lower[1] + block.y - 1) /
                  block.y,
              maxblocks[1]),
          1);
      hip_parallel_launch<ClosureType, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 3) {
      dim3 const block(m_policy.m_tile[0], m_policy.m_tile[1],
                       m_policy.m_tile[2]);
      dim3 const grid(
          std::min<array_index_type>(
              (m_policy.m_upper[0] - m_policy.m_lower[0] + block.x - 1) /
                  block.x,
              maxblocks[0]),
          std::min<array_index_type>(
              (m_policy.m_upper[1] - m_policy.m_lower[1] + block.y - 1) /
                  block.y,
              maxblocks[1]),
          std::min<array_index_type>(
              (m_policy.m_upper[2] - m_policy.m_lower[2] + block.z - 1) /
                  block.z,
              maxblocks[2]));
      hip_parallel_launch<ClosureType, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 4) {
      // id0,id1 encoded within threadIdx.x; id2 to threadIdx.y; id3 to
      // threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2], m_policy.m_tile[3]);
      dim3 const grid(
          std::min<array_index_type>(
              m_policy.m_tile_end[0] * m_policy.m_tile_end[1], maxblocks[0]),
          std::min<array_index_type>(
              (m_policy.m_upper[2] - m_policy.m_lower[2] + block.y - 1) /
                  block.y,
              maxblocks[1]),
          std::min<array_index_type>(
              (m_policy.m_upper[3] - m_policy.m_lower[3] + block.z - 1) /
                  block.z,
              maxblocks[2]));
      hip_parallel_launch<ClosureType, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 5) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4
      // to threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2] * m_policy.m_tile[3],
                       m_policy.m_tile[4]);
      dim3 const grid(
          std::min<array_index_type>(
              m_policy.m_tile_end[0] * m_policy.m_tile_end[1], maxblocks[0]),
          std::min<array_index_type>(
              m_policy.m_tile_end[2] * m_policy.m_tile_end[3], maxblocks[1]),
          std::min<array_index_type>(
              (m_policy.m_upper[4] - m_policy.m_lower[4] + block.z - 1) /
                  block.z,
              maxblocks[2]));
      hip_parallel_launch<ClosureType, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else if (Policy::rank == 6) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y;
      // id4,id5 to threadIdx.z
      dim3 const block(m_policy.m_tile[0] * m_policy.m_tile[1],
                       m_policy.m_tile[2] * m_policy.m_tile[3],
                       m_policy.m_tile[4] * m_policy.m_tile[5]);
      dim3 const grid(
          std::min<array_index_type>(
              m_policy.m_tile_end[0] * m_policy.m_tile_end[1], maxblocks[0]),
          std::min<array_index_type>(
              m_policy.m_tile_end[2] * m_policy.m_tile_end[3], maxblocks[1]),
          std::min<array_index_type>(
              m_policy.m_tile_end[4] * m_policy.m_tile_end[5], maxblocks[2]));
      hip_parallel_launch<ClosureType, LaunchBounds>(
          *this, grid, block, 0,
          m_policy.space().impl_internal_space_instance(), false);
    } else {
      Kokkos::abort("Kokkos::MDRange Error: Exceeded rank bounds with HIP\n");
    }

  }  // end execute

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    using closure_type =
        ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>, HIP>;
    unsigned block_size = hip_get_max_blocksize<closure_type, LaunchBounds>();
    if (block_size == 0)
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelFor< HIP > could not find a valid "
                      "tile size."));
    return block_size;
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
