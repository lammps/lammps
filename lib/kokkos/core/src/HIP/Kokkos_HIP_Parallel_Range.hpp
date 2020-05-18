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

#ifndef KOKKO_HIP_PARALLEL_RANGE_HPP
#define KOKKO_HIP_PARALLEL_RANGE_HPP

#include <Kokkos_Parallel.hpp>

#if defined(__HIPCC__)

#include <HIP/Kokkos_HIP_BlockSize_Deduction.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_ReduceScan.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                  Kokkos::Experimental::HIP> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using LaunchBounds = typename Policy::launch_bounds;

  const FunctorType m_functor;
  const Policy m_policy;

  ParallelFor()        = delete;
  ParallelFor& operator=(const ParallelFor&) = delete;

  template <class TagType>
  inline __device__
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_range(const Member i) const {
    m_functor(i);
  }

  template <class TagType>
  inline __device__
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_range(const Member i) const {
    m_functor(TagType(), i);
  }

 public:
  using functor_type = FunctorType;

  inline __device__ void operator()(void) const {
    const Member work_stride = hipBlockDim_y * hipGridDim_x;
    const Member work_end    = m_policy.end();

    for (Member iwork =
             m_policy.begin() + hipThreadIdx_y + hipBlockDim_y * hipBlockIdx_x;
         iwork < work_end;
         iwork = iwork < work_end - work_stride ? iwork + work_stride
                                                : work_end) {
      this->template exec_range<WorkTag>(iwork);
    }
  }

  inline void execute() const {
    const typename Policy::index_type nwork = m_policy.end() - m_policy.begin();

    const int block_size = 256;  // FIXME_HIP Choose block_size better
    const dim3 block(1, block_size, 1);
    const dim3 grid(
        typename Policy::index_type((nwork + block.y - 1) / block.y), 1, 1);

    Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelFor, LaunchBounds>(
        *this, grid, block, 0, m_policy.space().impl_internal_space_instance(),
        false);
  }

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::HIP> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using WorkRange    = typename Policy::WorkRange;
  using WorkTag      = typename Policy::work_tag;
  using Member       = typename Policy::member_type;
  using LaunchBounds = typename Policy::launch_bounds;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type;

  using ValueTraits =
      Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>;
  using ValueInit = Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd>;
  using ValueJoin = Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using value_type     = typename ValueTraits::value_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Kokkos::Experimental::HIP::size_type;
  using index_type     = typename Policy::index_type;

  // Algorithmic constraints: blockSize is a power of two AND hipBlockDim_y ==
  // hipBlockDim_z == 1

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  size_type* m_scratch_space = nullptr;
  size_type* m_scratch_flags = nullptr;

  // Shall we use the shfl based reduction or not (only use it for static sized
  // types of more than 128bit)
  enum {
    UseShflReduction = false
  };  //((sizeof(value_type)>2*sizeof(double)) && ValueTraits::StaticValueSize)
      //};
      // Some crutch to do function overloading
 private:
  using DummyShflReductionType  = double;
  using DummySHMEMReductionType = int;

 public:
  // Make the exec_range calls call to Reduce::DeviceIterateTile
  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update) const {
    m_functor(i, update);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update) const {
    m_functor(TagType(), i, update);
  }

  __device__ inline void operator()() const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(
                       ReducerConditional::select(m_functor, m_reducer)) /
                   sizeof(size_type));

    {
      reference_type value = ValueInit::init(
          ReducerConditional::select(m_functor, m_reducer),
          ::Kokkos::Experimental::kokkos_impl_hip_shared_memory<size_type>() +
              hipThreadIdx_y * word_count.value);

      // Number of blocks is bounded so that the reduction can be limited to two
      // passes. Each thread block is given an approximately equal amount of
      // work to perform. Accumulate the values for this block. The accumulation
      // ordering does not match the final pass, but is arithmetically
      // equivalent.

      const WorkRange range(m_policy, hipBlockIdx_x, hipGridDim_x);

      for (Member iwork     = range.begin() + hipThreadIdx_y,
                  iwork_end = range.end();
           iwork < iwork_end; iwork += hipBlockDim_y) {
        this->template exec_range<WorkTag>(iwork, value);
      }
    }

    // Reduce with final value at hipBlockDim_y - 1 location.
    if (hip_single_inter_block_reduce_scan<false, ReducerTypeFwd, WorkTagFwd>(
            ReducerConditional::select(m_functor, m_reducer), hipBlockIdx_x,
            hipGridDim_x,
            ::Kokkos::Experimental::kokkos_impl_hip_shared_memory<size_type>(),
            m_scratch_space, m_scratch_flags)) {
      // This is the final block with the final result at the final threads'
      // location

      size_type* const shared =
          ::Kokkos::Experimental::kokkos_impl_hip_shared_memory<size_type>() +
          (hipBlockDim_y - 1) * word_count.value;
      size_type* const global = m_result_ptr_device_accessible
                                    ? reinterpret_cast<size_type*>(m_result_ptr)
                                    : m_scratch_space;

      if (hipThreadIdx_y == 0) {
        Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
            ReducerConditional::select(m_functor, m_reducer), shared);
      }

      if (::Kokkos::Experimental::Impl::HIPTraits::WarpSize <
          word_count.value) {
        __syncthreads();
      }

      for (unsigned i = hipThreadIdx_y; i < word_count.value;
           i += hipBlockDim_y) {
        global[i] = shared[i];
      }
    }
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    // FIXME_HIP I don't know where 8 comes from
    unsigned int n = ::Kokkos::Experimental::Impl::HIPTraits::WarpSize * 8;
    int shmem_size =
        hip_single_inter_block_reduce_scan_shmem<false, FunctorType, WorkTag>(
            f, n);
    while (
        (n &&
         (m_policy.space().impl_internal_space_instance()->m_maxShmemPerBlock <
          shmem_size)) ||
        (n > static_cast<unsigned int>(
                 Kokkos::Experimental::Impl::hip_get_max_block_size<
                     ParallelReduce, LaunchBounds>(f, 1, shmem_size, 0)))) {
      n >>= 1;
      shmem_size =
          hip_single_inter_block_reduce_scan_shmem<false, FunctorType, WorkTag>(
              f, n);
    }
    return n;
  }

  inline void execute() {
    const index_type nwork = m_policy.end() - m_policy.begin();
    if (nwork) {
      const int block_size = local_block_size(m_functor);

      m_scratch_space =
          ::Kokkos::Experimental::Impl::hip_internal_scratch_space(
              ValueTraits::value_size(
                  ReducerConditional::select(m_functor, m_reducer)) *
              block_size /* block_size == max block_count */);
      m_scratch_flags =
          ::Kokkos::Experimental::Impl::hip_internal_scratch_flags(
              sizeof(size_type));

      // REQUIRED ( 1 , N , 1 )
      const dim3 block(1, block_size, 1);
      // Required grid.x <= block.y
      const dim3 grid(
          std::min(int(block.y), int((nwork + block.y - 1) / block.y)), 1, 1);

      const int shmem =
          UseShflReduction
              ? 0
              : hip_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                         WorkTag>(m_functor,
                                                                  block.y);

      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelReduce,
                                                    LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      if (!m_result_ptr_device_accessible) {
        ::Kokkos::Experimental::HIP().fence();

        if (m_result_ptr) {
          const int size = ValueTraits::value_size(
              ReducerConditional::select(m_functor, m_reducer));
          DeepCopy<HostSpace, ::Kokkos::Experimental::HIPSpace>(
              m_result_ptr, m_scratch_space, size);
        }
      }
    } else {
      if (m_result_ptr) {
        ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                        m_result_ptr);
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ViewType& arg_result,
                 typename std::enable_if<Kokkos::is_view<ViewType>::value,
                                         void*>::type = NULL)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                              typename ViewType::memory_space>::accessible) {}

  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::Experimental::HIPSpace,
                              typename ReducerType::result_view_type::
                                  memory_space>::accessible) {}
};

template <class FunctorType, class... Traits>
class ParallelScanHIPBase {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 protected:
  using Member       = typename Policy::member_type;
  using WorkTag      = typename Policy::work_tag;
  using WorkRange    = typename Policy::WorkRange;
  using LaunchBounds = typename Policy::launch_bounds;

  using ValueTraits = Kokkos::Impl::FunctorValueTraits<FunctorType, WorkTag>;
  using ValueInit   = Kokkos::Impl::FunctorValueInit<FunctorType, WorkTag>;
  using ValueOps    = Kokkos::Impl::FunctorValueOps<FunctorType, WorkTag>;

 public:
  using pointer_type   = typename ValueTraits::pointer_type;
  using reference_type = typename ValueTraits::reference_type;
  using functor_type   = FunctorType;
  using size_type      = Kokkos::Experimental::HIP::size_type;
  using index_type     = typename Policy::index_type;

 protected:
  // Algorithmic constraints:
  //  (a) hipBlockDim_y is a power of two
  //  (b) hipBlockDim_x == hipBlockDim_z == 1
  //  (c) hipGridDim_x  <= hipBlockDim_y * hipBlockDim_y
  //  (d) hipGridDim_y  == hipGridDim_z == 1

  const FunctorType m_functor;
  const Policy m_policy;
  size_type* m_scratch_space = nullptr;
  size_type* m_scratch_flags = nullptr;
  size_type m_final          = false;
  int m_grid_x               = 0;

 private:
  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update,
                 const bool final_result) const {
    m_functor(i, update, final_result);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_range(const Member& i, reference_type update,
                 const bool final_result) const {
    m_functor(TagType(), i, update, final_result);
  }

  //----------------------------------------

  __device__ inline void initial(void) const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(m_functor) / sizeof(size_type));

    size_type* const shared_value =
        Kokkos::Experimental::kokkos_impl_hip_shared_memory<size_type>() +
        word_count.value * hipThreadIdx_y;

    ValueInit::init(m_functor, shared_value);

    // Number of blocks is bounded so that the reduction can be limited to two
    // passes. Each thread block is given an approximately equal amount of work
    // to perform. Accumulate the values for this block. The accumulation
    // ordering does not match the final pass, but is arithmetically equivalent.

    const WorkRange range(m_policy, hipBlockIdx_x, hipGridDim_x);

    for (Member iwork = range.begin() + hipThreadIdx_y, iwork_end = range.end();
         iwork < iwork_end; iwork += hipBlockDim_y) {
      this->template exec_range<WorkTag>(
          iwork, ValueOps::reference(shared_value), false);
    }

    // Reduce and scan, writing out scan of blocks' totals and block-groups'
    // totals. Blocks' scan values are written to 'hipBlockIdx_x' location.
    // Block-groups' scan values are at: i = ( j * hipBlockDim_y - 1 ) for i <
    // hipGridDim_x
    hip_single_inter_block_reduce_scan<true, FunctorType, WorkTag>(
        m_functor, hipBlockIdx_x, hipGridDim_x,
        Kokkos::Experimental::kokkos_impl_hip_shared_memory<size_type>(),
        m_scratch_space, m_scratch_flags);
  }

  //----------------------------------------

  __device__ inline void final(void) const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(size_type)>
        word_count(ValueTraits::value_size(m_functor) / sizeof(size_type));

    // Use shared memory as an exclusive scan: { 0 , value[0] , value[1] ,
    // value[2] , ... }
    size_type* const shared_data =
        Kokkos::Experimental::kokkos_impl_hip_shared_memory<size_type>();
    size_type* const shared_prefix =
        shared_data + word_count.value * hipThreadIdx_y;
    size_type* const shared_accum =
        shared_data + word_count.value * (hipBlockDim_y + 1);

    // Starting value for this thread block is the previous block's total.
    if (hipBlockIdx_x) {
      size_type* const block_total =
          m_scratch_space + word_count.value * (hipBlockIdx_x - 1);
      for (unsigned i = hipThreadIdx_y; i < word_count.value; ++i) {
        shared_accum[i] = block_total[i];
      }
    } else if (0 == hipThreadIdx_y) {
      ValueInit::init(m_functor, shared_accum);
    }

    const WorkRange range(m_policy, hipBlockIdx_x, hipGridDim_x);

    for (typename Policy::member_type iwork_base = range.begin();
         iwork_base < range.end(); iwork_base += hipBlockDim_y) {
      const typename Policy::member_type iwork = iwork_base + hipThreadIdx_y;

      __syncthreads();  // Don't overwrite previous iteration values until they
                        // are used

      ValueInit::init(m_functor, shared_prefix + word_count.value);

      // Copy previous block's accumulation total into thread[0] prefix and
      // inclusive scan value of this block
      for (unsigned i = hipThreadIdx_y; i < word_count.value; ++i) {
        shared_data[i + word_count.value] = shared_data[i] = shared_accum[i];
      }

      // Make sure the write is seen by all threads
      __threadfence_block();

      // Call functor to accumulate inclusive scan value for this work item
      const bool doWork = (iwork < range.end());
      if (doWork) {
        this->template exec_range<WorkTag>(
            iwork, ValueOps::reference(shared_prefix + word_count.value),
            false);
      }

      // Scan block values into locations shared_data[1..hipBlockDim_y]
      hip_intra_block_reduce_scan<true, FunctorType, WorkTag>(
          m_functor,
          typename ValueTraits::pointer_type(shared_data + word_count.value));

      {
        size_type* const block_total =
            shared_data + word_count.value * hipBlockDim_y;
        for (unsigned i = hipThreadIdx_y; i < word_count.value; ++i) {
          shared_accum[i] = block_total[i];
        }
      }

      // Call functor with exclusive scan value
      if (doWork) {
        this->template exec_range<WorkTag>(
            iwork, ValueOps::reference(shared_prefix), true);
      }
    }
  }

 public:
  //----------------------------------------

  __device__ inline void operator()(void) const {
    if (!m_final) {
      initial();
    } else {
      final();
    }
  }

  // Determine block size constrained by shared memory:
  inline unsigned local_block_size(const FunctorType& f) {
    // hipBlockDim_y must be power of two = 128 (2 warps) or 256 (4 warps) or
    // 512 (8 warps) hipGridDim_x <= hipBlockDim_y * hipBlockDim_y
    //
    // TODO check best option

    unsigned n = Experimental::Impl::HIPTraits::WarpSize * 4;
    while (n && static_cast<unsigned>(m_policy.space()
                                          .impl_internal_space_instance()
                                          ->m_maxShmemPerBlock) <
                    hip_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                             WorkTag>(f, n)) {
      n >>= 1;
    }
    return n;
  }

  inline void impl_execute() {
    const index_type nwork = m_policy.end() - m_policy.begin();
    if (nwork) {
      // FIXME_HIP we cannot choose it larger for large work sizes to work
      // correctly, the unit tests fail with wrong results
      const int gridMaxComputeCapability_2x = 0x01fff;

      // FIXME_HIP block sizes greater than 256 don't work correctly,
      // the unit tests fail with wrong results
      const int block_size =
          std::min(static_cast<int>(local_block_size(m_functor)), 256);

      const int grid_max =
          std::min(block_size * block_size, gridMaxComputeCapability_2x);

      // At most 'max_grid' blocks:
      const int max_grid =
          std::min<int>(grid_max, (nwork + block_size - 1) / block_size);

      // How much work per block:
      const int work_per_block = (nwork + max_grid - 1) / max_grid;

      // How many block are really needed for this much work:
      m_grid_x = (nwork + work_per_block - 1) / work_per_block;

      m_scratch_space = Kokkos::Experimental::Impl::hip_internal_scratch_space(
          ValueTraits::value_size(m_functor) * m_grid_x);
      m_scratch_flags = Kokkos::Experimental::Impl::hip_internal_scratch_flags(
          sizeof(size_type) * 1);

      dim3 grid(m_grid_x, 1, 1);
      dim3 block(1, block_size, 1);  // REQUIRED DIMENSIONS ( 1 , N , 1 )
      const int shmem = ValueTraits::value_size(m_functor) * (block_size + 2);

      m_final = false;
      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelScanHIPBase,
                                                    LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute

      m_final = true;
      Kokkos::Experimental::Impl::HIPParallelLaunch<ParallelScanHIPBase,
                                                    LaunchBounds>(
          *this, grid, block, shmem,
          m_policy.space().impl_internal_space_instance(),
          false);  // copy to device and execute
    }
  }

  ParallelScanHIPBase(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Experimental::HIP>
    : private ParallelScanHIPBase<FunctorType, Traits...> {
 public:
  using Base = ParallelScanHIPBase<FunctorType, Traits...>;
  using Base::operator();

  inline void execute() { Base::impl_execute(); }

  ParallelScan(const FunctorType& arg_functor,
               const typename Base::Policy& arg_policy)
      : Base(arg_functor, arg_policy) {}
};

//----------------------------------------------------------------------------

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Experimental::HIP>
    : private ParallelScanHIPBase<FunctorType, Traits...> {
 public:
  using Base = ParallelScanHIPBase<FunctorType, Traits...>;
  using Base::operator();

  ReturnType& m_returnvalue;

  inline void execute() {
    Base::impl_execute();

    const auto nwork = Base::m_policy.end() - Base::m_policy.begin();
    if (nwork) {
      const int size = Base::ValueTraits::value_size(Base::m_functor);
      DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace>(
          &m_returnvalue,
          Base::m_scratch_space + (Base::m_grid_x - 1) * size / sizeof(int),
          size);
    }
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
