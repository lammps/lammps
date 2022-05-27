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

#ifndef KOKKO_SYCL_PARALLEL_SCAN_HPP
#define KOKKO_SYCL_PARALLEL_SCAN_HPP

#include <Kokkos_Macros.hpp>
#include <memory>
#include <vector>
#if defined(KOKKOS_ENABLE_SYCL)

namespace Kokkos {
namespace Impl {

// Perform a scan over a workgroup.
// At the end of this function, the subgroup scans are stored in the local array
// such that the last value (at position n_active_subgroups-1) contains the
// total sum.
template <class ValueJoin, class ValueInit, int dim, typename ValueType,
          typename FunctorType>
void workgroup_scan(sycl::nd_item<dim> item, const FunctorType& functor,
                    sycl::local_ptr<ValueType> local_mem,
                    ValueType& local_value, unsigned int global_range) {
  // subgroup scans
  auto sg                = item.get_sub_group();
  const auto sg_group_id = sg.get_group_id()[0];
  const auto id_in_sg    = sg.get_local_id()[0];
  for (unsigned int stride = 1; stride < global_range; stride <<= 1) {
    auto tmp = sg.shuffle_up(local_value, stride);
    if (id_in_sg >= stride) ValueJoin::join(functor, &local_value, &tmp);
  }

  const auto max_subgroup_size = sg.get_max_local_range()[0];
  const auto n_active_subgroups =
      (global_range + max_subgroup_size - 1) / max_subgroup_size;

  const auto local_range = sg.get_local_range()[0];
  if (id_in_sg == local_range - 1 && sg_group_id < n_active_subgroups)
    local_mem[sg_group_id] = local_value;
  local_value = sg.shuffle_up(local_value, 1);
  if (id_in_sg == 0) ValueInit::init(functor, &local_value);
  sycl::group_barrier(item.get_group());

  // scan subgroup results using the first subgroup
  if (n_active_subgroups > 1) {
    if (sg_group_id == 0) {
      const auto n_rounds =
          (n_active_subgroups + local_range - 1) / local_range;
      for (unsigned int round = 0; round < n_rounds; ++round) {
        const unsigned int idx = id_in_sg + round * local_range;
        const auto upper_bound =
            std::min(local_range, n_active_subgroups - round * local_range);
        auto local_sg_value = local_mem[idx < n_active_subgroups ? idx : 0];
        for (unsigned int stride = 1; stride < upper_bound; stride <<= 1) {
          auto tmp = sg.shuffle_up(local_sg_value, stride);
          if (id_in_sg >= stride) {
            if (idx < n_active_subgroups)
              ValueJoin::join(functor, &local_sg_value, &tmp);
            else
              local_sg_value = tmp;
          }
        }
        if (idx < n_active_subgroups) {
          local_mem[idx] = local_sg_value;
          if (round > 0)
            ValueJoin::join(functor, &local_mem[idx],
                            &local_mem[round * local_range - 1]);
        }
        if (round + 1 < n_rounds) sycl::group_barrier(sg);
      }
    }
    sycl::group_barrier(item.get_group());
  }

  // add results to all subgroups
  if (sg_group_id > 0)
    ValueJoin::join(functor, &local_value, &local_mem[sg_group_id - 1]);
}

template <class FunctorType, class... Traits>
class ParallelScanSYCLBase {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 protected:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using WorkRange    = typename Policy::WorkRange;
  using LaunchBounds = typename Policy::launch_bounds;

  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, WorkTag>;
  using ValueInit   = Kokkos::Impl::FunctorValueInit<FunctorType, WorkTag>;
  using ValueJoin   = Kokkos::Impl::FunctorValueJoin<FunctorType, WorkTag>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using value_type     = typename ValueTraits::value_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Kokkos::Experimental::SYCL::size_type;
  using index_type     = typename Policy::index_type;

 protected:
  const FunctorType m_functor;
  const Policy m_policy;
  pointer_type m_scratch_space = nullptr;

  // Only let one Parallel/Scan modify the shared memory. The
  // constructor acquires the mutex which is released in the destructor.
  std::scoped_lock<std::mutex> m_shared_memory_lock;

 private:
  template <typename FunctorWrapper>
  void scan_internal(sycl::queue& q, const FunctorWrapper& functor_wrapper,
                     pointer_type global_mem, std::size_t size) const {
    // FIXME_SYCL optimize
    constexpr size_t wgroup_size = 128;
    auto n_wgroups               = (size + wgroup_size - 1) / wgroup_size;
    pointer_type group_results   = global_mem + n_wgroups * wgroup_size;

    auto local_scans = q.submit([&](sycl::handler& cgh) {
      // Store subgroup totals
      const auto min_subgroup_size =
          q.get_device()
              .template get_info<sycl::info::device::sub_group_sizes>()
              .front();
      sycl::accessor<value_type, 1, sycl::access::mode::read_write,
                     sycl::access::target::local>
          local_mem(sycl::range<1>((wgroup_size + min_subgroup_size - 1) /
                                   min_subgroup_size),
                    cgh);

      cgh.parallel_for(
          sycl::nd_range<1>(n_wgroups * wgroup_size, wgroup_size),
          [=](sycl::nd_item<1> item) {
            const auto local_id  = item.get_local_linear_id();
            const auto global_id = item.get_global_linear_id();

            // Initialize local memory
            value_type local_value;
            if (global_id < size)
              local_value = global_mem[global_id];
            else
              ValueInit::init(functor_wrapper.get_functor(), &local_value);

            workgroup_scan<ValueJoin, ValueInit>(
                item, functor_wrapper.get_functor(), local_mem.get_pointer(),
                local_value, wgroup_size);

            if (n_wgroups > 1 && local_id == wgroup_size - 1)
              group_results[item.get_group_linear_id()] =
                  local_mem[item.get_sub_group().get_group_range()[0] - 1];

            // Write results to global memory
            if (global_id < size) global_mem[global_id] = local_value;
          });
    });
    q.ext_oneapi_submit_barrier(std::vector<sycl::event>{local_scans});

    if (n_wgroups > 1) {
      scan_internal(q, functor_wrapper, group_results, n_wgroups);
      auto update_with_group_results = q.submit([&](sycl::handler& cgh) {
        cgh.parallel_for(
            sycl::nd_range<1>(n_wgroups * wgroup_size, wgroup_size),
            [=](sycl::nd_item<1> item) {
              const auto global_id = item.get_global_linear_id();
              if (global_id < size)
                ValueJoin::join(functor_wrapper.get_functor(),
                                &global_mem[global_id],
                                &group_results[item.get_group_linear_id()]);
            });
      });
      q.ext_oneapi_submit_barrier(
          std::vector<sycl::event>{update_with_group_results});
    }
  }

  template <typename FunctorWrapper>
  sycl::event sycl_direct_launch(const FunctorWrapper& functor_wrapper,
                                 sycl::event memcpy_event) const {
    // Convenience references
    const Kokkos::Experimental::SYCL& space = m_policy.space();
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *space.impl_internal_space_instance();
    sycl::queue& q = *instance.m_queue;

    const std::size_t len = m_policy.end() - m_policy.begin();

    // Initialize global memory
    auto initialize_global_memory = q.submit([&](sycl::handler& cgh) {
      auto global_mem = m_scratch_space;
      auto begin      = m_policy.begin();

      cgh.depends_on(memcpy_event);
      cgh.parallel_for(sycl::range<1>(len), [=](sycl::item<1> item) {
        const typename Policy::index_type id =
            static_cast<typename Policy::index_type>(item.get_id()) + begin;
        value_type update{};
        ValueInit::init(functor_wrapper.get_functor(), &update);
        if constexpr (std::is_same<WorkTag, void>::value)
          functor_wrapper.get_functor()(id, update, false);
        else
          functor_wrapper.get_functor()(WorkTag(), id, update, false);
        global_mem[id] = update;
      });
    });
    q.ext_oneapi_submit_barrier(
        std::vector<sycl::event>{initialize_global_memory});

    // Perform the actual exclusive scan
    scan_internal(q, functor_wrapper, m_scratch_space, len);

    // Write results to global memory
    auto update_global_results = q.submit([&](sycl::handler& cgh) {
      auto global_mem = m_scratch_space;
      cgh.parallel_for(sycl::range<1>(len), [=](sycl::item<1> item) {
        auto global_id = item.get_id(0);

        value_type update = global_mem[global_id];
        if constexpr (std::is_same<WorkTag, void>::value)
          functor_wrapper.get_functor()(global_id, update, true);
        else
          functor_wrapper.get_functor()(WorkTag(), global_id, update, true);
        global_mem[global_id] = update;
      });
    });
    q.ext_oneapi_submit_barrier(
        std::vector<sycl::event>{update_global_results});
    return update_global_results;
  }

 public:
  template <typename PostFunctor>
  void impl_execute(const PostFunctor& post_functor) {
    if (m_policy.begin() == m_policy.end()) return;

    auto& instance        = *m_policy.space().impl_internal_space_instance();
    const std::size_t len = m_policy.end() - m_policy.begin();

    // Compute the total amount of memory we will need. We emulate the recursive
    // structure that is used to do the actual scan. Essentially, we need to
    // allocate memory for the whole range and then recursively for the reduced
    // group results until only one group is left.
    std::size_t total_memory = 0;
    {
      size_t wgroup_size   = 128;
      size_t n_nested_size = len;
      size_t n_nested_wgroups;
      do {
        n_nested_wgroups = (n_nested_size + wgroup_size - 1) / wgroup_size;
        n_nested_size    = n_nested_wgroups;
        total_memory += sizeof(value_type) * n_nested_wgroups * wgroup_size;
      } while (n_nested_wgroups > 1);
      total_memory += sizeof(value_type) * wgroup_size;
    }

    // FIXME_SYCL consider only storing one value per block and recreate initial
    // results in the end before doing the final pass
    m_scratch_space =
        static_cast<pointer_type>(instance.scratch_space(total_memory));

    Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem&
        indirectKernelMem = instance.get_indirect_kernel_mem();

    auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);

    sycl::event event =
        sycl_direct_launch(functor_wrapper, functor_wrapper.get_copy_event());
    functor_wrapper.register_event(event);
    post_functor();
  }

  ParallelScanSYCLBase(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_shared_memory_lock(m_policy.space()
                                 .impl_internal_space_instance()
                                 ->m_mutexScratchSpace) {}
};

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Experimental::SYCL>
    : private ParallelScanSYCLBase<FunctorType, Traits...> {
 public:
  using Base = ParallelScanSYCLBase<FunctorType, Traits...>;

  inline void execute() {
    Base::impl_execute([]() {});
  }

  ParallelScan(const FunctorType& arg_functor,
               const typename Base::Policy& arg_policy)
      : Base(arg_functor, arg_policy) {}
};

//----------------------------------------------------------------------------

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Experimental::SYCL>
    : private ParallelScanSYCLBase<FunctorType, Traits...> {
 public:
  using Base = ParallelScanSYCLBase<FunctorType, Traits...>;

  ReturnType& m_returnvalue;

  inline void execute() {
    Base::impl_execute([&]() {
      const long long nwork = Base::m_policy.end() - Base::m_policy.begin();
      if (nwork > 0) {
        const int size = Base::ValueTraits::value_size(Base::m_functor);
        DeepCopy<HostSpace, Kokkos::Experimental::SYCLDeviceUSMSpace>(
            &m_returnvalue, Base::m_scratch_space + nwork - 1, size);
      }
    });
  }

  ParallelScanWithTotal(const FunctorType& arg_functor,
                        const typename Base::Policy& arg_policy,
                        ReturnType& arg_returnvalue)
      : Base(arg_functor, arg_policy), m_returnvalue(arg_returnvalue) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif

#endif
