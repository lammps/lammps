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

#ifndef KOKKOS_SYCL_WORKGROUP_REDUCTION_HPP
#define KOKKOS_SYCL_WORKGROUP_REDUCTION_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos::Impl::SYCLReduction {

// FIXME_SYCL It appears that using shuffles is slower than going through local
// memory.
template <class ReducerType>
inline constexpr bool use_shuffle_based_algorithm = false;
// std::is_reference_v<typename ReducerType::reference_type>;

template <typename ValueType, typename ReducerType, int dim>
std::enable_if_t<!use_shuffle_based_algorithm<ReducerType>> workgroup_reduction(
    sycl::nd_item<dim>& item, sycl::local_accessor<ValueType> local_mem,
    sycl::device_ptr<ValueType> results_ptr,
    sycl::global_ptr<ValueType> device_accessible_result_ptr,
    const unsigned int value_count_, const ReducerType& final_reducer,
    bool final, unsigned int max_size) {
  const unsigned int value_count =
      std::is_reference_v<typename ReducerType::reference_type> ? 1
                                                                : value_count_;
  const int local_id = item.get_local_linear_id();

  // Perform the actual workgroup reduction in each subgroup
  // separately.
  auto sg            = item.get_sub_group();
  auto* result       = &local_mem[local_id * value_count];
  const int id_in_sg = sg.get_local_id()[0];
  const auto local_range =
      std::min<unsigned int>(sg.get_local_range()[0], max_size);
  const auto upper_stride_bound =
      std::min<unsigned int>(local_range - id_in_sg, max_size - local_id);
  for (unsigned int stride = 1; stride < local_range; stride <<= 1) {
    if (stride < upper_stride_bound)
      final_reducer.join(result, &local_mem[(local_id + stride) * value_count]);
    sycl::group_barrier(sg);
  }
  sycl::group_barrier(item.get_group());

  // Do the final reduction only using the first subgroup.
  if (sg.get_group_id()[0] == 0) {
    const unsigned int n_subgroups = sg.get_group_range()[0];
    const int max_subgroup_size    = sg.get_max_local_range()[0];
    auto* result_ = &local_mem[id_in_sg * max_subgroup_size * value_count];
    // In case the number of subgroups is larger than the range of
    // the first subgroup, we first combine the items with a higher
    // index.
    for (unsigned int offset = local_range; offset < n_subgroups;
         offset += local_range)
      if (id_in_sg + offset < n_subgroups)
        final_reducer.join(
            result_,
            &local_mem[(id_in_sg + offset) * max_subgroup_size * value_count]);
    sycl::group_barrier(sg);

    // Then, we proceed as before.
    for (unsigned int stride = 1; stride < local_range; stride <<= 1) {
      if (id_in_sg + stride < n_subgroups)
        final_reducer.join(
            result_,
            &local_mem[(id_in_sg + stride) * max_subgroup_size * value_count]);
      sycl::group_barrier(sg);
    }

    // Finally, we copy the workgroup results back to global memory
    // to be used in the next iteration. If this is the last
    // iteration, i.e., there is only one workgroup also call
    // final() if necessary.
    if (id_in_sg == 0) {
      if (final) {
        final_reducer.final(&local_mem[0]);
        if (device_accessible_result_ptr != nullptr)
          final_reducer.copy(&device_accessible_result_ptr[0], &local_mem[0]);
        else
          final_reducer.copy(&results_ptr[0], &local_mem[0]);
      } else
        final_reducer.copy(
            &results_ptr[(item.get_group_linear_id()) * value_count],
            &local_mem[0]);
    }
  }
}

template <typename ValueType, typename ReducerType, int dim>
std::enable_if_t<use_shuffle_based_algorithm<ReducerType>> workgroup_reduction(
    sycl::nd_item<dim>& item, sycl::local_accessor<ValueType> local_mem,
    ValueType local_value, sycl::device_ptr<ValueType> results_ptr,
    sycl::global_ptr<ValueType> device_accessible_result_ptr,
    const ReducerType& final_reducer, bool final, unsigned int max_size) {
  const auto local_id = item.get_local_linear_id();

  // Perform the actual workgroup reduction in each subgroup
  // separately.
  auto sg            = item.get_sub_group();
  const int id_in_sg = sg.get_local_id()[0];
  const auto local_range =
      std::min<unsigned int>(sg.get_local_range()[0], max_size);

  const auto upper_stride_bound =
      std::min<unsigned int>(local_range - id_in_sg, max_size - local_id);
  for (unsigned int stride = 1; stride < local_range; stride <<= 1) {
    auto tmp = sg.shuffle_down(local_value, stride);
    if (stride < upper_stride_bound) final_reducer.join(&local_value, &tmp);
  }

  // Copy the subgroup results into the first positions of the
  // reduction array.
  const int max_subgroup_size = sg.get_max_local_range()[0];
  const int n_active_subgroups =
      (max_size + max_subgroup_size - 1) / max_subgroup_size;
  const int sg_group_id = sg.get_group_id()[0];
  if (id_in_sg == 0 && sg_group_id <= n_active_subgroups)
    local_mem[sg_group_id] = local_value;

  item.barrier(sycl::access::fence_space::local_space);

  // Do the final reduction only using the first subgroup.
  if (sg.get_group_id()[0] == 0) {
    auto sg_value = local_mem[id_in_sg < n_active_subgroups ? id_in_sg : 0];

    // In case the number of subgroups is larger than the range of
    // the first subgroup, we first combine the items with a higher
    // index.
    if (n_active_subgroups > local_range) {
      for (unsigned int offset = local_range; offset < n_active_subgroups;
           offset += local_range)
        if (id_in_sg + offset < n_active_subgroups) {
          final_reducer.join(&sg_value, &local_mem[(id_in_sg + offset)]);
        }
      sg.barrier();
    }

    // Then, we proceed as before.
    for (unsigned int stride = 1; stride < local_range; stride <<= 1) {
      auto tmp = sg.shuffle_down(sg_value, stride);
      if (id_in_sg + stride < n_active_subgroups)
        final_reducer.join(&sg_value, &tmp);
    }

    // Finally, we copy the workgroup results back to global memory
    // to be used in the next iteration. If this is the last
    // iteration, i.e., there is only one workgroup also call
    // final() if necessary.
    if (id_in_sg == 0) {
      if (final) {
        final_reducer.final(&sg_value);
        if (device_accessible_result_ptr != nullptr)
          device_accessible_result_ptr[0] = sg_value;
        else
          results_ptr[0] = sg_value;
      } else
        results_ptr[(item.get_group_linear_id())] = sg_value;
    }
  }
}

}  // namespace Kokkos::Impl::SYCLReduction

#endif /* KOKKOS_SYCL_WORKGROUP_REDUCTION_HPP */
