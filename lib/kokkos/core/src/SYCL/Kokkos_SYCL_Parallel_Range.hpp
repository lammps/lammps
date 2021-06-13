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

#ifndef KOKKOS_SYCL_PARALLEL_RANGE_HPP_
#define KOKKOS_SYCL_PARALLEL_RANGE_HPP_

#include <impl/KokkosExp_IterateTileGPU.hpp>

template <class FunctorType, class ExecPolicy>
class Kokkos::Impl::ParallelFor<FunctorType, ExecPolicy,
                                Kokkos::Experimental::SYCL> {
 public:
  using Policy = ExecPolicy;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

  const FunctorType m_functor;
  const Policy m_policy;

  template <typename Functor>
  static void sycl_direct_launch(const Policy& policy, const Functor& functor) {
    // Convenience references
    const Kokkos::Experimental::SYCL& space = policy.space();
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *space.impl_internal_space_instance();
    sycl::queue& q = *instance.m_queue;

    space.fence();

    q.submit([functor, policy](sycl::handler& cgh) {
      sycl::range<1> range(policy.end() - policy.begin());
      const auto begin = policy.begin();

      cgh.parallel_for(range, [=](sycl::item<1> item) {
        const typename Policy::index_type id = item.get_linear_id() + begin;
        if constexpr (std::is_same<WorkTag, void>::value)
          functor(id);
        else
          functor(WorkTag(), id);
      });
    });

    space.fence();
  }

 public:
  using functor_type = FunctorType;

  void execute() const {
    if (m_policy.begin() == m_policy.end()) return;

    Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem&
        indirectKernelMem = m_policy.space()
                                .impl_internal_space_instance()
                                ->m_indirectKernelMem;

    const auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);
    sycl_direct_launch(m_policy, functor_wrapper.get_functor());
  }

  ParallelFor(const ParallelFor&) = delete;
  ParallelFor(ParallelFor&&)      = delete;
  ParallelFor& operator=(const ParallelFor&) = delete;
  ParallelFor& operator=(ParallelFor&&) = delete;
  ~ParallelFor()                        = default;

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

// ParallelFor
template <class FunctorType, class... Traits>
class Kokkos::Impl::ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                                Kokkos::Experimental::SYCL> {
 public:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;
  using LaunchBounds     = typename Policy::launch_bounds;
  using WorkTag          = typename Policy::work_tag;

  const FunctorType m_functor;
  // MDRangePolicy is not trivially copyable. Hence, replicate the data we
  // really need in DeviceIterateTile in a trivially copyable struct.
  const struct BarePolicy {
    using index_type = typename Policy::index_type;

    BarePolicy(const Policy& policy)
        : m_lower(policy.m_lower),
          m_upper(policy.m_upper),
          m_tile(policy.m_tile),
          m_tile_end(policy.m_tile_end),
          m_num_tiles(policy.m_num_tiles) {}

    const typename Policy::point_type m_lower;
    const typename Policy::point_type m_upper;
    const typename Policy::tile_type m_tile;
    const typename Policy::point_type m_tile_end;
    const typename Policy::index_type m_num_tiles;
    static constexpr Iterate inner_direction = Policy::inner_direction;
  } m_policy;
  const Kokkos::Experimental::SYCL& m_space;

  sycl::nd_range<3> compute_ranges() const {
    const auto& m_tile     = m_policy.m_tile;
    const auto& m_tile_end = m_policy.m_tile_end;

    if constexpr (Policy::rank == 2) {
      sycl::range<3> local_sizes(m_tile[0], m_tile[1], 1);
      sycl::range<3> global_sizes(m_tile_end[0] * m_tile[0],
                                  m_tile_end[1] * m_tile[1], 1);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 3) {
      sycl::range<3> local_sizes(m_tile[0], m_tile[1], m_tile[2]);
      sycl::range<3> global_sizes(m_tile_end[0] * m_tile[0],
                                  m_tile_end[1] * m_tile[1],
                                  m_tile_end[2] * m_tile[2]);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 4) {
      // id0,id1 encoded within first index; id2 to second index; id3 to third
      // index
      sycl::range<3> local_sizes(m_tile[0] * m_tile[1], m_tile[2], m_tile[3]);
      sycl::range<3> global_sizes(
          m_tile_end[0] * m_tile[0] * m_tile_end[1] * m_tile[1],
          m_tile_end[2] * m_tile[2], m_tile_end[3] * m_tile[3]);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 5) {
      // id0,id1 encoded within first index; id2,id3 to second index; id4 to
      // third index
      sycl::range<3> local_sizes(m_tile[0] * m_tile[1], m_tile[2] * m_tile[3],
                                 m_tile[4]);
      sycl::range<3> global_sizes(
          m_tile_end[0] * m_tile[0] * m_tile_end[1] * m_tile[1],
          m_tile_end[2] * m_tile[2] * m_tile_end[3] * m_tile[3],
          m_tile_end[4] * m_tile[4]);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 6) {
      // id0,id1 encoded within first index; id2,id3 to second index; id4,id5 to
      // third index
      sycl::range<3> local_sizes(m_tile[0] * m_tile[1], m_tile[2] * m_tile[3],
                                 m_tile[4] * m_tile[5]);
      sycl::range<3> global_sizes(
          m_tile_end[0] * m_tile[0] * m_tile_end[1] * m_tile[1],
          m_tile_end[2] * m_tile[2] * m_tile_end[3] * m_tile[3],
          m_tile_end[4] * m_tile[4] * m_tile_end[5] * m_tile[5]);
      return {global_sizes, local_sizes};
    }
    static_assert(Policy::rank > 1 && Policy::rank < 7,
                  "Kokkos::MDRange Error: Exceeded rank bounds with SYCL\n");
  }

  template <typename Functor>
  void sycl_direct_launch(const Functor& functor) const {
    // Convenience references
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_space.impl_internal_space_instance();
    sycl::queue& q = *instance.m_queue;

    m_space.fence();

    if (m_policy.m_num_tiles == 0) return;

    const BarePolicy bare_policy(m_policy);

    q.submit([functor, this, bare_policy](sycl::handler& cgh) {
      const auto range = compute_ranges();

      cgh.parallel_for(range, [functor, bare_policy](sycl::nd_item<3> item) {
        const index_type local_x    = item.get_local_id(0);
        const index_type local_y    = item.get_local_id(1);
        const index_type local_z    = item.get_local_id(2);
        const index_type global_x   = item.get_group(0);
        const index_type global_y   = item.get_group(1);
        const index_type global_z   = item.get_group(2);
        const index_type n_global_x = item.get_group_range(0);
        const index_type n_global_y = item.get_group_range(1);
        const index_type n_global_z = item.get_group_range(2);

        Kokkos::Impl::DeviceIterateTile<Policy::rank, BarePolicy, Functor,
                                        typename Policy::work_tag>(
            bare_policy, functor, {n_global_x, n_global_y, n_global_z},
            {global_x, global_y, global_z}, {local_x, local_y, local_z})
            .exec_range();
      });
    });

    m_space.fence();
  }

 public:
  using functor_type = FunctorType;

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& policy, const Functor&) {
    return policy.space().impl_internal_space_instance()->m_maxWorkgroupSize;
  }

  void execute() const {
    Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem&
        indirectKernelMem =
            m_space.impl_internal_space_instance()->m_indirectKernelMem;

    const auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);
    sycl_direct_launch(functor_wrapper.get_functor());
  }

  ParallelFor(const ParallelFor&) = delete;
  ParallelFor(ParallelFor&&)      = delete;
  ParallelFor& operator=(const ParallelFor&) = delete;
  ParallelFor& operator=(ParallelFor&&) = delete;
  ~ParallelFor()                        = default;

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_space(arg_policy.space()) {}
};

#endif  // KOKKOS_SYCL_PARALLEL_RANGE_HPP_
