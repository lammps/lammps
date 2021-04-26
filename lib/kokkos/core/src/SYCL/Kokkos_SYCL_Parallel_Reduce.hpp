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

#ifndef KOKKOS_SYCL_PARALLEL_REDUCE_HPP
#define KOKKOS_SYCL_PARALLEL_REDUCE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_SYCL)

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::SYCL> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::REDUCE, Policy, FunctorType>;
  using execution_space = typename Analysis::execution_space;
  using value_type      = typename Analysis::value_type;
  using pointer_type    = typename Analysis::pointer_type;
  using reference_type  = typename Analysis::reference_type;

  using WorkTag = typename Policy::work_tag;

 public:
  // V - View
  template <typename V>
  ParallelReduce(
      const FunctorType& f, const Policy& p, const V& v,
      typename std::enable_if<Kokkos::is_view<V>::value, void*>::type = nullptr)
      : m_functor(f), m_policy(p), m_result_ptr(v.data()) {}

  ParallelReduce(const FunctorType& f, const Policy& p,
                 const ReducerType& reducer)
      : m_functor(f),
        m_policy(p),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {}

 private:
  template <typename PolicyType, typename Functor, typename Reducer>
  void sycl_direct_launch(const PolicyType& policy, const Functor& functor,
                          const Reducer& reducer) const {
    using ReducerConditional =
        Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                           FunctorType, ReducerType>;
    using ReducerTypeFwd = typename ReducerConditional::type;
    using WorkTagFwd =
        std::conditional_t<std::is_same<InvalidType, ReducerType>::value,
                           WorkTag, void>;
    using ValueInit =
        Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
    using ValueJoin =
        Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;
    using ValueOps = Kokkos::Impl::FunctorValueOps<FunctorType, WorkTag>;

    auto selected_reducer = ReducerConditional::select(functor, reducer);

    // Convenience references
    const Kokkos::Experimental::SYCL& space = policy.space();
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *space.impl_internal_space_instance();
    sycl::queue& q = *instance.m_queue;

    // FIXME_SYCL optimize
    constexpr size_t wgroup_size       = 128;
    constexpr size_t values_per_thread = 2;
    std::size_t size                   = policy.end() - policy.begin();
    const auto init_size               = std::max<std::size_t>(
        ((size + values_per_thread - 1) / values_per_thread + wgroup_size - 1) /
            wgroup_size,
        1);
    const unsigned int value_count =
        FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>::value_count(
            selected_reducer);
    // FIXME_SYCL only use the first half
    const auto results_ptr = static_cast<pointer_type>(instance.scratch_space(
        sizeof(value_type) * std::max(value_count, 1u) * init_size * 2));
    // FIXME_SYCL without this we are running into a race condition
    const auto results_ptr2 =
        results_ptr + std::max(value_count, 1u) * init_size;

    // If size<=1 we only call init(), the functor and possibly final once
    // working with the global scratch memory but don't copy back to
    // m_result_ptr yet.
    if (size <= 1) {
      q.submit([&](sycl::handler& cgh) {
        const auto begin = policy.begin();
        cgh.single_task([=]() {
          const auto& selected_reducer = ReducerConditional::select(
              static_cast<const FunctorType&>(functor),
              static_cast<const ReducerType&>(reducer));
          reference_type update =
              ValueInit::init(selected_reducer, results_ptr);
          if (size == 1) {
            if constexpr (std::is_same<WorkTag, void>::value)
              functor(begin, update);
            else
              functor(WorkTag(), begin, update);
          }
          if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
            FunctorFinal<FunctorType, WorkTag>::final(
                static_cast<const FunctorType&>(functor), results_ptr);
        });
      });
      space.fence();
    }

    // Otherwise, we perform a reduction on the values in all workgroups
    // separately, write the workgroup results back to global memory and recurse
    // until only one workgroup does the reduction and thus gets the final
    // value.
    bool first_run = true;
    while (size > 1) {
      auto n_wgroups = ((size + values_per_thread - 1) / values_per_thread +
                        wgroup_size - 1) /
                       wgroup_size;
      q.submit([&](sycl::handler& cgh) {
        sycl::accessor<value_type, 1, sycl::access::mode::read_write,
                       sycl::access::target::local>
            local_mem(sycl::range<1>(wgroup_size) * std::max(value_count, 1u),
                      cgh);
        const auto begin = policy.begin();

        cgh.parallel_for(
            sycl::nd_range<1>(n_wgroups * wgroup_size, wgroup_size),
            [=](sycl::nd_item<1> item) {
              const auto local_id = item.get_local_linear_id();
              const auto global_id =
                  wgroup_size * item.get_group_linear_id() * values_per_thread +
                  local_id;
              const auto& selected_reducer = ReducerConditional::select(
                  static_cast<const FunctorType&>(functor),
                  static_cast<const ReducerType&>(reducer));

              // In the first iteration, we call functor to initialize the local
              // memory. Otherwise, the local memory is initialized with the
              // results from the previous iteration that are stored in global
              // memory. Note that we load values_per_thread values per thread
              // and immediately combine them to avoid too many threads being
              // idle in the actual workgroup reduction.
              using index_type       = typename Policy::index_type;
              const auto upper_bound = std::min<index_type>(
                  global_id + values_per_thread * wgroup_size, size);
              if (first_run) {
                reference_type update = ValueInit::init(
                    selected_reducer, &local_mem[local_id * value_count]);
                for (index_type id = global_id; id < upper_bound;
                     id += wgroup_size) {
                  if constexpr (std::is_same<WorkTag, void>::value)
                    functor(id + begin, update);
                  else
                    functor(WorkTag(), id + begin, update);
                }
              } else {
                if (global_id >= size)
                  ValueInit::init(selected_reducer,
                                  &local_mem[local_id * value_count]);
                else {
                  ValueOps::copy(functor, &local_mem[local_id * value_count],
                                 &results_ptr[global_id * value_count]);
                  for (index_type id = global_id + wgroup_size;
                       id < upper_bound; id += wgroup_size) {
                    ValueJoin::join(selected_reducer,
                                    &local_mem[local_id * value_count],
                                    &results_ptr[id * value_count]);
                  }
                }
              }
              item.barrier(sycl::access::fence_space::local_space);

              // Perform the actual workgroup reduction. To achieve a better
              // memory access pattern, we use sequential addressing and a
              // reversed loop. If the workgroup size is 8, the first element
              // contains all the values with index%4==0, after the second one
              // the values with index%2==0 and after the third one index%1==0,
              // i.e., all values.
              for (unsigned int stride = wgroup_size / 2; stride > 0;
                   stride >>= 1) {
                const auto idx = local_id;
                if (idx < stride) {
                  ValueJoin::join(selected_reducer,
                                  &local_mem[idx * value_count],
                                  &local_mem[(idx + stride) * value_count]);
                }
                item.barrier(sycl::access::fence_space::local_space);
              }

              // Finally, we copy the workgroup results back to global memory to
              // be used in the next iteration. If this is the last iteration,
              // i.e., there is only one workgroup also call final() if
              // necessary.
              if (local_id == 0) {
                ValueOps::copy(
                    functor,
                    &results_ptr2[(item.get_group_linear_id()) * value_count],
                    &local_mem[0]);
                if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
                  if (n_wgroups <= 1)
                    FunctorFinal<FunctorType, WorkTag>::final(
                        static_cast<const FunctorType&>(functor),
                        &results_ptr2[(item.get_group_linear_id()) *
                                      value_count]);
              }
            });
      });
      space.fence();

      // FIXME_SYCL this is likely not necessary, see above
      Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                             Kokkos::Experimental::SYCLDeviceUSMSpace>(
          space, results_ptr, results_ptr2,
          sizeof(*m_result_ptr) * value_count * n_wgroups);
      space.fence();

      first_run = false;
      size      = n_wgroups;
    }

    // At this point, the reduced value is written to the entry in results_ptr
    // and all that is left is to copy it back to the given result pointer if
    // necessary.
    if (m_result_ptr) {
      Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                             Kokkos::Experimental::SYCLDeviceUSMSpace>(
          space, m_result_ptr, results_ptr,
          sizeof(*m_result_ptr) * value_count);
      space.fence();
    }
  }

 public:
  void execute() const {
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_policy.space().impl_internal_space_instance();
    using IndirectKernelMem =
        Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem;
    IndirectKernelMem& indirectKernelMem  = instance.m_indirectKernelMem;
    IndirectKernelMem& indirectReducerMem = instance.m_indirectReducerMem;

    const auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);
    const auto reducer_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_reducer, indirectReducerMem);

    sycl_direct_launch(m_policy, functor_wrapper.get_functor(),
                       reducer_wrapper.get_functor());
  }

 private:
  FunctorType m_functor;
  Policy m_policy;
  ReducerType m_reducer;
  pointer_type m_result_ptr;
};

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::SYCL> {
 public:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

 private:
  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::REDUCE, Policy, FunctorType>;
  using execution_space = typename Analysis::execution_space;
  using value_type      = typename Analysis::value_type;
  using pointer_type    = typename Analysis::pointer_type;
  using reference_type  = typename Analysis::reference_type;

  using WorkTag = typename Policy::work_tag;

  // MDRangePolicy is not trivially copyable. Hence, replicate the data we
  // really need in DeviceIterateTile in a trivially copyable struct.
  struct BarePolicy {
    using index_type = typename Policy::index_type;

    BarePolicy(const Policy& policy)
        : m_lower(policy.m_lower),
          m_upper(policy.m_upper),
          m_tile(policy.m_tile),
          m_tile_end(policy.m_tile_end),
          m_num_tiles(policy.m_num_tiles),
          m_prod_tile_dims(policy.m_prod_tile_dims) {}

    const typename Policy::point_type m_lower;
    const typename Policy::point_type m_upper;
    const typename Policy::tile_type m_tile;
    const typename Policy::point_type m_tile_end;
    const typename Policy::index_type m_num_tiles;
    const typename Policy::index_type m_prod_tile_dims;
    static constexpr Iterate inner_direction = Policy::inner_direction;
    static constexpr int rank                = Policy::rank;
  };

 public:
  // V - View
  template <typename V>
  ParallelReduce(
      const FunctorType& f, const Policy& p, const V& v,
      typename std::enable_if<Kokkos::is_view<V>::value, void*>::type = nullptr)
      : m_functor(f), m_policy(p), m_space(p.space()), m_result_ptr(v.data()) {}

  ParallelReduce(const FunctorType& f, const Policy& p,
                 const ReducerType& reducer)
      : m_functor(f),
        m_policy(p),
        m_space(p.space()),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {}

 private:
  template <typename PolicyType, typename Functor, typename Reducer>
  void sycl_direct_launch(const PolicyType& policy, const Functor& functor,
                          const Reducer& reducer) const {
    using ReducerConditional =
        Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                           FunctorType, ReducerType>;
    using ReducerTypeFwd = typename ReducerConditional::type;
    using WorkTagFwd =
        std::conditional_t<std::is_same<InvalidType, ReducerType>::value,
                           WorkTag, void>;
    using ValueInit =
        Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
    using ValueJoin =
        Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;
    using ValueOps = Kokkos::Impl::FunctorValueOps<FunctorType, WorkTag>;

    // Convenience references
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_space.impl_internal_space_instance();
    sycl::queue& q = *instance.m_queue;

    const int nwork = m_policy.m_num_tiles;
    const int block_size =
        std::pow(2, std::ceil(std::log2(m_policy.m_prod_tile_dims)));

    const sycl::range<1> local_range(block_size);
    // REMEMBER swap local x<->y to be conforming with Cuda/HIP implementation
    const sycl::range<1> global_range(nwork * block_size);
    const sycl::nd_range<1> range{global_range, local_range};

    const size_t wgroup_size = range.get_local_range().size();
    size_t size              = range.get_global_range().size();
    const auto init_size =
        std::max<std::size_t>((size + wgroup_size - 1) / wgroup_size, 1);
    const auto& selected_reducer = ReducerConditional::select(functor, reducer);
    const unsigned int value_count =
        FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>::value_count(
            selected_reducer);
    // FIXME_SYCL only use the first half
    const auto results_ptr = static_cast<pointer_type>(instance.scratch_space(
        sizeof(value_type) * std::max(value_count, 1u) * init_size * 2));
    // FIXME_SYCL without this we are running into a race condition
    const auto results_ptr2 =
        results_ptr + std::max(value_count, 1u) * init_size;

    // If size<=1 we only call init(), the functor and possibly final once
    // working with the global scratch memory but don't copy back to
    // m_result_ptr yet.
    if (size <= 1) {
      q.submit([&](sycl::handler& cgh) {
        cgh.single_task([=]() {
          const auto& selected_reducer = ReducerConditional::select(
              static_cast<const FunctorType&>(functor),
              static_cast<const ReducerType&>(reducer));
          reference_type update =
              ValueInit::init(selected_reducer, results_ptr);
          if (size == 1) {
            Kokkos::Impl::Reduce::DeviceIterateTile<
                Policy::rank, BarePolicy, Functor, typename Policy::work_tag,
                reference_type>(policy, functor, update, {1, 1, 1}, {0, 0, 0},
                                {0, 0, 0})
                .exec_range();
          }
          if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
            FunctorFinal<FunctorType, WorkTag>::final(
                static_cast<const FunctorType&>(functor), results_ptr);
        });
      });
      m_space.fence();
    }

    // Otherwise, we perform a reduction on the values in all workgroups
    // separately, write the workgroup results back to global memory and recurse
    // until only one workgroup does the reduction and thus gets the final
    // value.
    bool first_run = true;
    while (size > 1) {
      auto n_wgroups = (size + wgroup_size - 1) / wgroup_size;
      q.submit([&](sycl::handler& cgh) {
        sycl::accessor<value_type, 1, sycl::access::mode::read_write,
                       sycl::access::target::local>
            local_mem(sycl::range<1>(wgroup_size) * std::max(value_count, 1u),
                      cgh);

        const BarePolicy bare_policy = m_policy;

        cgh.parallel_for(range, [=](sycl::nd_item<1> item) {
          const auto local_id = item.get_local_linear_id();
          const auto global_id =
              wgroup_size * item.get_group_linear_id() + local_id;
          const auto& selected_reducer = ReducerConditional::select(
              static_cast<const FunctorType&>(functor),
              static_cast<const ReducerType&>(reducer));

          // In the first iteration, we call functor to initialize the local
          // memory. Otherwise, the local memory is initialized with the
          // results from the previous iteration that are stored in global
          // memory.
          using index_type = typename Policy::index_type;
          const auto upper_bound =
              std::min<index_type>(global_id + wgroup_size, size);
          if (first_run) {
            reference_type update = ValueInit::init(
                selected_reducer, &local_mem[local_id * value_count]);

            // SWAPPED here to be conforming with CUDA implementation
            const index_type local_x    = 0;
            const index_type local_y    = item.get_local_id(0);
            const index_type local_z    = 0;
            const index_type global_x   = item.get_group(0);
            const index_type global_y   = 0;
            const index_type global_z   = 0;
            const index_type n_global_x = item.get_group_range(0);
            const index_type n_global_y = 1;
            const index_type n_global_z = 1;

            Kokkos::Impl::Reduce::DeviceIterateTile<
                Policy::rank, BarePolicy, Functor, typename Policy::work_tag,
                reference_type>(bare_policy, functor, update,
                                {n_global_x, n_global_y, n_global_z},
                                {global_x, global_y, global_z},
                                {local_x, local_y, local_z})
                .exec_range();
          } else {
            if (global_id >= size)
              ValueInit::init(selected_reducer,
                              &local_mem[local_id * value_count]);
            else {
              ValueOps::copy(functor, &local_mem[local_id * value_count],
                             &results_ptr[global_id * value_count]);
              for (index_type id = global_id + wgroup_size; id < upper_bound;
                   id += wgroup_size) {
                ValueJoin::join(selected_reducer,
                                &local_mem[local_id * value_count],
                                &results_ptr[id * value_count]);
              }
            }
          }
          item.barrier(sycl::access::fence_space::local_space);

          // Perform the actual workgroup reduction. To achieve a better
          // memory access pattern, we use sequential addressing and a
          // reversed loop. If the workgroup size is 8, the first element
          // contains all the values with index%4==0, after the second one
          // the values with index%2==0 and after the third one index%1==0,
          // i.e., all values.
          for (unsigned int stride = wgroup_size / 2; stride > 0;
               stride >>= 1) {
            const auto idx = local_id;
            if (idx < stride) {
              ValueJoin::join(selected_reducer, &local_mem[idx * value_count],
                              &local_mem[(idx + stride) * value_count]);
            }
            item.barrier(sycl::access::fence_space::local_space);
          }

          // Finally, we copy the workgroup results back to global memory to
          // be used in the next iteration. If this is the last iteration,
          // i.e., there is only one workgroup also call final() if
          // necessary.
          if (local_id == 0) {
            ValueOps::copy(
                functor,
                &results_ptr2[(item.get_group_linear_id()) * value_count],
                &local_mem[0]);
            if constexpr (ReduceFunctorHasFinal<FunctorType>::value)
              if (n_wgroups <= 1)
                FunctorFinal<FunctorType, WorkTag>::final(
                    static_cast<const FunctorType&>(functor),
                    &results_ptr2[(item.get_group_linear_id()) * value_count]);
          }
        });
      });
      m_space.fence();

      // FIXME_SYCL this is likely not necessary, see above
      Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                             Kokkos::Experimental::SYCLDeviceUSMSpace>(
          m_space, results_ptr, results_ptr2,
          sizeof(*m_result_ptr) * value_count * n_wgroups);
      m_space.fence();

      first_run = false;
      size      = n_wgroups;
    }

    // At this point, the reduced value is written to the entry in results_ptr
    // and all that is left is to copy it back to the given result pointer if
    // necessary.
    if (m_result_ptr) {
      Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                             Kokkos::Experimental::SYCLDeviceUSMSpace>(
          m_space, m_result_ptr, results_ptr,
          sizeof(*m_result_ptr) * value_count);
      m_space.fence();
    }
  }

 public:
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& policy, const Functor&) {
    return policy.space().impl_internal_space_instance()->m_maxThreadsPerSM;
  }

  void execute() const {
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_space.impl_internal_space_instance();
    using IndirectKernelMem =
        Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem;
    IndirectKernelMem& indirectKernelMem  = instance.m_indirectKernelMem;
    IndirectKernelMem& indirectReducerMem = instance.m_indirectReducerMem;

    const auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);
    const auto reducer_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_reducer, indirectReducerMem);

    sycl_direct_launch(m_policy, functor_wrapper.get_functor(),
                       reducer_wrapper.get_functor());
  }

 private:
  FunctorType m_functor;
  BarePolicy m_policy;
  const Kokkos::Experimental::SYCL& m_space;
  ReducerType m_reducer;
  pointer_type m_result_ptr;
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
#endif /* KOKKOS_SYCL_PARALLEL_REDUCE_HPP */
