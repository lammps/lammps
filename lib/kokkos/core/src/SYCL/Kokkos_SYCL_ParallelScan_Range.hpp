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

#ifndef KOKKOS_SYCL_PARALLEL_SCAN_RANGE_HPP
#define KOKKOS_SYCL_PARALLEL_SCAN_RANGE_HPP

#include <Kokkos_Macros.hpp>
#include <SYCL/Kokkos_SYCL_WorkgroupReduction.hpp>
#include <memory>
#include <vector>

namespace Kokkos::Impl {

// Perform a scan over a workgroup.
// At the end of this function, the subgroup scans are stored in the local array
// such that the last value (at position n_active_subgroups-1) contains the
// total sum.
template <int dim, typename ValueType, typename FunctorType>
void workgroup_scan(sycl::nd_item<dim> item, const FunctorType& final_reducer,
                    sycl::local_accessor<ValueType> local_mem,
                    ValueType& local_value, int global_range) {
  // subgroup scans
  auto sg               = item.get_sub_group();
  const int sg_group_id = sg.get_group_id()[0];
  const int id_in_sg    = sg.get_local_id()[0];
  const int local_range = std::min<int>(sg.get_local_range()[0], global_range);

#if defined(KOKKOS_ARCH_INTEL_GPU) || defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  auto shuffle_combine = [&](int stride) {
    if (stride < local_range) {
      auto tmp = Kokkos::Impl::SYCLReduction::shift_group_right(sg, local_value,
                                                                stride);
      if (id_in_sg >= stride) final_reducer.join(&local_value, &tmp);
    }
  };
  shuffle_combine(1);
  shuffle_combine(2);
  shuffle_combine(4);
  shuffle_combine(8);
  shuffle_combine(16);
  KOKKOS_ASSERT(local_range <= 32);
#else
  for (int stride = 1; stride < local_range; stride <<= 1) {
    auto tmp =
        Kokkos::Impl::SYCLReduction::shift_group_right(sg, local_value, stride);
    if (id_in_sg >= stride) final_reducer.join(&local_value, &tmp);
  }
#endif

  const int max_subgroup_size = sg.get_max_local_range()[0];
  const int n_active_subgroups =
      (global_range + max_subgroup_size - 1) / max_subgroup_size;

  if (id_in_sg == local_range - 1 && sg_group_id < n_active_subgroups)
    local_mem[sg_group_id] = local_value;
  local_value =
      Kokkos::Impl::SYCLReduction::shift_group_right(sg, local_value, 1);
  if (id_in_sg == 0) final_reducer.init(&local_value);
  sycl::group_barrier(item.get_group());

  // scan subgroup results using the first subgroup
  if (n_active_subgroups > 1) {
    if (sg_group_id == 0) {
      const int n_rounds = (n_active_subgroups + local_range - 1) / local_range;
      for (int round = 0; round < n_rounds; ++round) {
        const int idx = id_in_sg + round * local_range;
        const auto upper_bound =
            std::min(local_range, n_active_subgroups - round * local_range);
        auto local_sg_value = local_mem[idx < n_active_subgroups ? idx : 0];
#if defined(KOKKOS_ARCH_INTEL_GPU) || defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
        auto shuffle_combine_sg = [&](int stride) {
          if (stride < upper_bound) {
            auto tmp = Kokkos::Impl::SYCLReduction::shift_group_right(
                sg, local_sg_value, stride);
            if (id_in_sg >= stride) {
              if (idx < n_active_subgroups)
                final_reducer.join(&local_sg_value, &tmp);
              else
                local_sg_value = tmp;
            }
          }
        };
        shuffle_combine_sg(1);
        shuffle_combine_sg(2);
        shuffle_combine_sg(4);
        shuffle_combine_sg(8);
        shuffle_combine_sg(16);
        KOKKOS_ASSERT(upper_bound <= 32);
#else
        for (int stride = 1; stride < upper_bound; stride <<= 1) {
          auto tmp = Kokkos::Impl::SYCLReduction::shift_group_right(
              sg, local_sg_value, stride);
          if (id_in_sg >= stride) {
            if (idx < n_active_subgroups)
              final_reducer.join(&local_sg_value, &tmp);
            else
              local_sg_value = tmp;
          }
        }
#endif
        if (idx < n_active_subgroups) {
          local_mem[idx] = local_sg_value;
          if (round > 0)
            final_reducer.join(&local_mem[idx],
                               &local_mem[round * local_range - 1]);
        }
        if (round + 1 < n_rounds) sycl::group_barrier(sg);
      }
    }
    sycl::group_barrier(item.get_group());
  }

  // add results to all subgroups
  if (sg_group_id > 0)
    final_reducer.join(&local_value, &local_mem[sg_group_id - 1]);
}

template <class FunctorType, class ValueType, class... Traits>
class ParallelScanSYCLBase {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 protected:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

 public:
  using Analysis       = FunctorAnalysis<FunctorPatternInterface::SCAN, Policy,
                                   FunctorType, ValueType>;
  using pointer_type   = typename Analysis::pointer_type;
  using value_type     = typename Analysis::value_type;
  using reference_type = typename Analysis::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Kokkos::Experimental::SYCL::size_type;
  using index_type     = typename Policy::index_type;

 protected:
  const CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>
      m_functor_reducer;
  const Policy m_policy;
  sycl_host_ptr<value_type> m_scratch_host = nullptr;
  pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;

 private:
  template <typename FunctorWrapper>
  sycl::event sycl_direct_launch(const FunctorWrapper& functor_wrapper,
                                 sycl::event memcpy_event) {
    // Convenience references
    const Kokkos::Experimental::SYCL& space = m_policy.space();
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *space.impl_internal_space_instance();
    sycl::queue& q = space.sycl_queue();

    const auto size = m_policy.end() - m_policy.begin();

    auto scratch_flags = static_cast<sycl_device_ptr<unsigned int>>(
        instance.scratch_flags(sizeof(unsigned int)));

    const auto begin = m_policy.begin();

    // Initialize global memory
    auto scan_lambda_factory = [&](sycl::local_accessor<value_type> local_mem,
                                   sycl::local_accessor<unsigned int>
                                       num_teams_done,
                                   sycl_device_ptr<value_type> global_mem_,
                                   sycl_device_ptr<value_type> group_results_) {
      auto lambda = [=](sycl::nd_item<1> item) {
        auto global_mem    = global_mem_;
        auto group_results = group_results_;

        const CombinedFunctorReducer<FunctorType, typename Analysis::Reducer>&
            functor_reducer        = functor_wrapper.get_functor();
        const FunctorType& functor = functor_reducer.get_functor();
        const typename Analysis::Reducer& reducer =
            functor_reducer.get_reducer();

        const auto n_wgroups  = item.get_group_range()[0];
        const int wgroup_size = item.get_local_range()[0];

        const int local_id         = item.get_local_linear_id();
        const index_type global_id = item.get_global_linear_id();

        // Initialize local memory
        value_type local_value;
        reducer.init(&local_value);
        if (global_id < size) {
          if constexpr (std::is_void<WorkTag>::value)
            functor(global_id + begin, local_value, false);
          else
            functor(WorkTag(), global_id + begin, local_value, false);
        }

        workgroup_scan<>(item, reducer, local_mem, local_value, wgroup_size);

        // Write results to global memory
        if (global_id < size) global_mem[global_id] = local_value;

        if (local_id == wgroup_size - 1) {
          group_results[item.get_group_linear_id()] =
              local_mem[item.get_sub_group().get_group_range()[0] - 1];

          sycl::atomic_ref<unsigned, sycl::memory_order::acq_rel,
                           sycl::memory_scope::device,
                           sycl::access::address_space::global_space>
              scratch_flags_ref(*scratch_flags);
          num_teams_done[0] = ++scratch_flags_ref;
        }
        item.barrier(sycl::access::fence_space::global_space);
        if (num_teams_done[0] == n_wgroups) {
          if (local_id == 0) *scratch_flags = 0;
          value_type total;
          reducer.init(&total);

          for (unsigned int offset = 0; offset < n_wgroups;
               offset += wgroup_size) {
            index_type id = local_id + offset;
            if (id < static_cast<index_type>(n_wgroups))
              local_value = group_results[id];
            else
              reducer.init(&local_value);
            workgroup_scan<>(
                item, reducer, local_mem, local_value,
                std::min<index_type>(n_wgroups - offset, wgroup_size));
            if (id < static_cast<index_type>(n_wgroups)) {
              reducer.join(&local_value, &total);
              group_results[id] = local_value;
            }
            reducer.join(
                &total,
                &local_mem[item.get_sub_group().get_group_range()[0] - 1]);
            if (offset + wgroup_size < n_wgroups)
              item.barrier(sycl::access::fence_space::global_space);
          }
        }
      };
      return lambda;
    };

    size_t wgroup_size;
    size_t n_wgroups;
    sycl_device_ptr<value_type> global_mem;
    sycl_device_ptr<value_type> group_results;

    desul::ensure_sycl_lock_arrays_on_device(q);

    auto perform_work_group_scans = q.submit([&](sycl::handler& cgh) {
      sycl::local_accessor<unsigned int> num_teams_done(1, cgh);

      auto dummy_scan_lambda =
          scan_lambda_factory({1, cgh}, num_teams_done, nullptr, nullptr);

      static sycl::kernel kernel = [&] {
        sycl::kernel_id functor_kernel_id =
            sycl::get_kernel_id<decltype(dummy_scan_lambda)>();
        auto kernel_bundle =
            sycl::get_kernel_bundle<sycl::bundle_state::executable>(
                q.get_context(), std::vector{functor_kernel_id});
        return kernel_bundle.get_kernel(functor_kernel_id);
      }();
      auto multiple = kernel.get_info<sycl::info::kernel_device_specific::
                                          preferred_work_group_size_multiple>(
          q.get_device());
      auto max =
          kernel.get_info<sycl::info::kernel_device_specific::work_group_size>(
              q.get_device());

      wgroup_size = static_cast<size_t>(max / multiple) * multiple;
      n_wgroups   = (size + wgroup_size - 1) / wgroup_size;

      // Compute the total amount of memory we will need.
      // We need to allocate memory for the whole range (rounded towards the
      // next multiple of the workgroup size) and for one element per workgroup
      // that will contain the sum of the previous workgroups totals.
      // FIXME_SYCL consider only storing one value per block and recreate
      // initial results in the end before doing the final pass
      global_mem =
          static_cast<sycl_device_ptr<value_type>>(instance.scratch_space(
              n_wgroups * (wgroup_size + 1) * sizeof(value_type)));
      m_scratch_host = static_cast<sycl_host_ptr<value_type>>(
          instance.scratch_host(sizeof(value_type)));

      group_results = global_mem + n_wgroups * wgroup_size;

      // Store subgroup totals in local space
      const auto min_subgroup_size =
          q.get_device()
              .template get_info<sycl::info::device::sub_group_sizes>()
              .front();
      sycl::local_accessor<value_type> local_mem(
          sycl::range<1>((wgroup_size + min_subgroup_size - 1) /
                         min_subgroup_size),
          cgh);

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      cgh.depends_on(memcpy_event);
#else
      (void)memcpy_event;
#endif

      auto scan_lambda = scan_lambda_factory(local_mem, num_teams_done,
                                             global_mem, group_results);
      cgh.parallel_for(sycl::nd_range<1>(n_wgroups * wgroup_size, wgroup_size),
                       scan_lambda);
    });

    // Write results to global memory
    auto update_global_results = q.submit([&](sycl::handler& cgh) {
      // The compiler failed with CL_INVALID_ARG_VALUE if using m_result_ptr
      // directly.
      pointer_type result_ptr = m_result_ptr_device_accessible
                                    ? m_result_ptr
                                    : static_cast<pointer_type>(m_scratch_host);

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      cgh.depends_on(perform_work_group_scans);
#endif

      cgh.parallel_for(
          sycl::nd_range<1>(n_wgroups * wgroup_size, wgroup_size),
          [=](sycl::nd_item<1> item) {
            const index_type global_id = item.get_global_linear_id();
            const CombinedFunctorReducer<
                FunctorType, typename Analysis::Reducer>& functor_reducer =
                functor_wrapper.get_functor();
            const FunctorType& functor = functor_reducer.get_functor();
            const typename Analysis::Reducer& reducer =
                functor_reducer.get_reducer();

            if (global_id < size) {
              value_type update = global_mem[global_id];

              reducer.join(&update, &group_results[item.get_group_linear_id()]);

              if constexpr (std::is_void<WorkTag>::value)
                functor(global_id + begin, update, true);
              else
                functor(WorkTag(), global_id + begin, update, true);

              if (global_id == size - 1) *result_ptr = update;
            }
          });
    });
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
    q.ext_oneapi_submit_barrier(
        std::vector<sycl::event>{update_global_results});
#endif
    return update_global_results;
  }

 public:
  template <typename PostFunctor>
  void impl_execute(const PostFunctor& post_functor) {
    if (m_policy.begin() == m_policy.end()) return;

    auto& instance = *m_policy.space().impl_internal_space_instance();

    // Only let one instance at a time resize the instance's scratch memory
    // allocations.
    std::scoped_lock<std::mutex> scratch_buffers_lock(
        instance.m_mutexScratchSpace);

    Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem&
        indirectKernelMem = instance.get_indirect_kernel_mem();

    auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor_reducer, indirectKernelMem);

    sycl::event event =
        sycl_direct_launch(functor_wrapper, functor_wrapper.get_copy_event());
    functor_wrapper.register_event(event);
    post_functor();
  }

  ParallelScanSYCLBase(const FunctorType& arg_functor, const Policy& arg_policy,
                       pointer_type arg_result_ptr,
                       bool arg_result_ptr_device_accessible)
      : m_functor_reducer(arg_functor, typename Analysis::Reducer{arg_functor}),
        m_policy(arg_policy),
        m_result_ptr(arg_result_ptr),
        m_result_ptr_device_accessible(arg_result_ptr_device_accessible) {}
};

}  // namespace Kokkos::Impl

template <class FunctorType, class... Traits>
class Kokkos::Impl::ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                                 Kokkos::Experimental::SYCL>
    : private ParallelScanSYCLBase<FunctorType, void, Traits...> {
 public:
  using Base = ParallelScanSYCLBase<FunctorType, void, Traits...>;

  inline void execute() {
    Base::impl_execute([]() {});
  }

  ParallelScan(const FunctorType& arg_functor,
               const typename Base::Policy& arg_policy)
      : Base(arg_functor, arg_policy, nullptr, false) {}
};

//----------------------------------------------------------------------------

template <class FunctorType, class ReturnType, class... Traits>
class Kokkos::Impl::ParallelScanWithTotal<
    FunctorType, Kokkos::RangePolicy<Traits...>, ReturnType,
    Kokkos::Experimental::SYCL>
    : public ParallelScanSYCLBase<FunctorType, ReturnType, Traits...> {
 public:
  using Base = ParallelScanSYCLBase<FunctorType, ReturnType, Traits...>;

  const Kokkos::Experimental::SYCL& m_exec;

  inline void execute() {
    Base::impl_execute([&]() {
      const long long nwork = Base::m_policy.end() - Base::m_policy.begin();
      if (nwork > 0 && !Base::m_result_ptr_device_accessible) {
        // Using DeepCopy instead of fence+memcpy turned out to be up to 2x
        // slower.
        m_exec.fence(
            "Kokkos::Impl::ParallelReduce<SYCL, MDRangePolicy>::execute: "
            "result not device-accessible");
        const int size = Base::m_functor_reducer.get_reducer().value_size();
        std::memcpy(Base::m_result_ptr, Base::m_scratch_host, size);
      }
    });
  }

  template <class ViewType>
  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const typename Base::Policy& arg_policy,
                        const ViewType& arg_result_view)
      : Base(arg_functor, arg_policy, arg_result_view.data(),
             MemorySpaceAccess<Experimental::SYCLDeviceUSMSpace,
                               typename ViewType::memory_space>::accessible),
        m_exec(arg_policy.space()) {}
};

#endif
