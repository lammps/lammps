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

#ifndef KOKKOS_CUDA_WORKGRAPHPOLICY_HPP
#define KOKKOS_CUDA_WORKGRAPHPOLICY_HPP

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::WorkGraphPolicy<Traits...>,
                  Kokkos::Cuda> {
 public:
  typedef Kokkos::WorkGraphPolicy<Traits...> Policy;
  typedef ParallelFor<FunctorType, Policy, Kokkos::Cuda> Self;

 private:
  Policy m_policy;
  FunctorType m_functor;

  template <class TagType>
  __device__ inline
      typename std::enable_if<std::is_same<TagType, void>::value>::type
      exec_one(const std::int32_t w) const noexcept {
    m_functor(w);
  }

  template <class TagType>
  __device__ inline
      typename std::enable_if<!std::is_same<TagType, void>::value>::type
      exec_one(const std::int32_t w) const noexcept {
    const TagType t{};
    m_functor(t, w);
  }

 public:
  __device__ inline void operator()() const noexcept {
    if (0 == (threadIdx.y % 16)) {
      // Spin until COMPLETED_TOKEN.
      // END_TOKEN indicates no work is currently available.

      for (std::int32_t w = Policy::END_TOKEN;
           Policy::COMPLETED_TOKEN != (w = m_policy.pop_work());) {
        if (Policy::END_TOKEN != w) {
          exec_one<typename Policy::work_tag>(w);
          m_policy.completed_work(w);
        }
      }
    }
  }

  inline void execute() {
    const int warps_per_block = 4;
    const dim3 grid(Kokkos::Impl::cuda_internal_multiprocessor_count(), 1, 1);
    const dim3 block(1, Kokkos::Impl::CudaTraits::WarpSize, warps_per_block);
    const int shared = 0;

    Kokkos::Impl::CudaParallelLaunch<Self>(
        *this, grid, block, shared, Cuda().impl_internal_space_instance(),
        false);
  }

  inline ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_policy(arg_policy), m_functor(arg_functor) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #define KOKKOS_CUDA_WORKGRAPHPOLICY_HPP */
