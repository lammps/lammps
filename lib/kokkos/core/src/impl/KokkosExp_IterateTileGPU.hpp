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

#ifndef KOKKOS_EXP_ITERATE_TILE_GPU_HPP
#define KOKKOS_EXP_ITERATE_TILE_GPU_HPP

#include <Kokkos_Macros.hpp>

#include <algorithm>

#include <utility>

#include <impl/Kokkos_Profiling_Interface.hpp>
#include <typeinfo>

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_SYCL
template <typename index_type>
struct EmulateCUDADim3 {
  index_type x;
  index_type y;
  index_type z;
};
#endif

template <class Tag, class Functor, class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION std::enable_if_t<std::is_void<Tag>::value>
_tag_invoke(Functor const& f, Args&&... args) {
  f((Args &&) args...);
}

template <class Tag, class Functor, class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION std::enable_if_t<!std::is_void<Tag>::value>
_tag_invoke(Functor const& f, Args&&... args) {
  f(Tag{}, (Args &&) args...);
}

template <class Tag, class Functor, class T, size_t N, size_t... Idxs,
          class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION void _tag_invoke_array_helper(
    Functor const& f, T (&vals)[N], std::integer_sequence<size_t, Idxs...>,
    Args&&... args) {
  _tag_invoke<Tag>(f, vals[Idxs]..., (Args &&) args...);
}

template <class Tag, class Functor, class T, size_t N, class... Args>
KOKKOS_IMPL_FORCEINLINE_FUNCTION void _tag_invoke_array(Functor const& f,
                                                        T (&vals)[N],
                                                        Args&&... args) {
  _tag_invoke_array_helper<Tag>(f, vals, std::make_index_sequence<N>{},
                                (Args &&) args...);
}

// ------------------------------------------------------------------ //
// ParallelFor iteration pattern
template <int N, typename PolicyType, typename Functor, typename Tag>
struct DeviceIterateTile;

// Rank 2
template <typename PolicyType, typename Functor, typename Tag>
struct DeviceIterateTile<2, PolicyType, Functor, Tag> {
  using index_type = typename PolicyType::index_type;

#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(
      const PolicyType& policy_, const Functor& f_,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_policy(policy_),
        m_func(f_),
        gridDim(gridDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(const PolicyType& policy_,
                                                const Functor& f_)
      : m_policy(policy_), m_func(f_) {}
#endif

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const {
    if (PolicyType::inner_direction == Iterate::Left) {
      // Loop over size maxnumblocks until full range covered
      for (index_type tile_id1 = static_cast<index_type>(blockIdx.y);
           tile_id1 < m_policy.m_tile_end[1]; tile_id1 += gridDim.y) {
        const index_type offset_1 =
            tile_id1 * m_policy.m_tile[1] +
            static_cast<index_type>(threadIdx.y) +
            static_cast<index_type>(m_policy.m_lower[1]);
        if (offset_1 < m_policy.m_upper[1] &&
            static_cast<index_type>(threadIdx.y) < m_policy.m_tile[1]) {
          for (index_type tile_id0 = static_cast<index_type>(blockIdx.x);
               tile_id0 < m_policy.m_tile_end[0]; tile_id0 += gridDim.x) {
            const index_type offset_0 =
                tile_id0 * m_policy.m_tile[0] +
                static_cast<index_type>(threadIdx.x) +
                static_cast<index_type>(m_policy.m_lower[0]);
            if (offset_0 < m_policy.m_upper[0] &&
                static_cast<index_type>(threadIdx.x) < m_policy.m_tile[0]) {
              Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1);
            }
          }
        }
      }
    } else {
      for (index_type tile_id0 = static_cast<index_type>(blockIdx.x);
           tile_id0 < m_policy.m_tile_end[0]; tile_id0 += gridDim.x) {
        const index_type offset_0 =
            tile_id0 * m_policy.m_tile[0] +
            static_cast<index_type>(threadIdx.x) +
            static_cast<index_type>(m_policy.m_lower[0]);
        if (offset_0 < m_policy.m_upper[0] &&
            static_cast<index_type>(threadIdx.x) < m_policy.m_tile[0]) {
          for (index_type tile_id1 = static_cast<index_type>(blockIdx.y);
               tile_id1 < m_policy.m_tile_end[1]; tile_id1 += gridDim.y) {
            const index_type offset_1 =
                tile_id1 * m_policy.m_tile[1] +
                static_cast<index_type>(threadIdx.y) +
                static_cast<index_type>(m_policy.m_lower[1]);
            if (offset_1 < m_policy.m_upper[1] &&
                static_cast<index_type>(threadIdx.y) < m_policy.m_tile[1]) {
              Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1);
            }
          }
        }
      }
    }
  }  // end exec_range

 private:
  const PolicyType& m_policy;
  const Functor& m_func;
#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif
};

// Rank 3
template <typename PolicyType, typename Functor, typename Tag>
struct DeviceIterateTile<3, PolicyType, Functor, Tag> {
  using index_type = typename PolicyType::index_type;

#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(
      const PolicyType& policy_, const Functor& f_,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_policy(policy_),
        m_func(f_),
        gridDim(gridDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(const PolicyType& policy_,
                                                const Functor& f_)
      : m_policy(policy_), m_func(f_) {}
#endif

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const {
    if (PolicyType::inner_direction == Iterate::Left) {
      for (index_type tile_id2 = static_cast<index_type>(blockIdx.z);
           tile_id2 < m_policy.m_tile_end[2]; tile_id2 += gridDim.z) {
        const index_type offset_2 =
            tile_id2 * m_policy.m_tile[2] +
            static_cast<index_type>(threadIdx.z) +
            static_cast<index_type>(m_policy.m_lower[2]);
        if (offset_2 < m_policy.m_upper[2] &&
            static_cast<index_type>(threadIdx.z) < m_policy.m_tile[2]) {
          for (index_type tile_id1 = static_cast<index_type>(blockIdx.y);
               tile_id1 < m_policy.m_tile_end[1]; tile_id1 += gridDim.y) {
            const index_type offset_1 =
                tile_id1 * m_policy.m_tile[1] +
                static_cast<index_type>(threadIdx.y) +
                static_cast<index_type>(m_policy.m_lower[1]);
            if (offset_1 < m_policy.m_upper[1] &&
                static_cast<index_type>(threadIdx.y) < m_policy.m_tile[1]) {
              for (index_type tile_id0 = static_cast<index_type>(blockIdx.x);
                   tile_id0 < m_policy.m_tile_end[0]; tile_id0 += gridDim.x) {
                const index_type offset_0 =
                    tile_id0 * m_policy.m_tile[0] +
                    static_cast<index_type>(threadIdx.x) +
                    static_cast<index_type>(m_policy.m_lower[0]);
                if (offset_0 < m_policy.m_upper[0] &&
                    static_cast<index_type>(threadIdx.x) < m_policy.m_tile[0]) {
                  Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1, offset_2);
                }
              }
            }
          }
        }
      }
    } else {
      for (index_type tile_id0 = static_cast<index_type>(blockIdx.x);
           tile_id0 < m_policy.m_tile_end[0]; tile_id0 += gridDim.x) {
        const index_type offset_0 =
            tile_id0 * m_policy.m_tile[0] +
            static_cast<index_type>(threadIdx.x) +
            static_cast<index_type>(m_policy.m_lower[0]);
        if (offset_0 < m_policy.m_upper[0] &&
            static_cast<index_type>(threadIdx.x) < m_policy.m_tile[0]) {
          for (index_type tile_id1 = static_cast<index_type>(blockIdx.y);
               tile_id1 < m_policy.m_tile_end[1]; tile_id1 += gridDim.y) {
            const index_type offset_1 =
                tile_id1 * m_policy.m_tile[1] +
                static_cast<index_type>(threadIdx.y) +
                static_cast<index_type>(m_policy.m_lower[1]);
            if (offset_1 < m_policy.m_upper[1] &&
                static_cast<index_type>(threadIdx.y) < m_policy.m_tile[1]) {
              for (index_type tile_id2 = static_cast<index_type>(blockIdx.z);
                   tile_id2 < m_policy.m_tile_end[2]; tile_id2 += gridDim.z) {
                const index_type offset_2 =
                    tile_id2 * m_policy.m_tile[2] +
                    static_cast<index_type>(threadIdx.z) +
                    static_cast<index_type>(m_policy.m_lower[2]);
                if (offset_2 < m_policy.m_upper[2] &&
                    static_cast<index_type>(threadIdx.z) < m_policy.m_tile[2]) {
                  Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1, offset_2);
                }
              }
            }
          }
        }
      }
    }
  }  // end exec_range

 private:
  const PolicyType& m_policy;
  const Functor& m_func;
#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif
};

// Rank 4
template <typename PolicyType, typename Functor, typename Tag>
struct DeviceIterateTile<4, PolicyType, Functor, Tag> {
  using index_type = typename PolicyType::index_type;

#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(
      const PolicyType& policy_, const Functor& f_,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_policy(policy_),
        m_func(f_),
        gridDim(gridDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(const PolicyType& policy_,
                                                const Functor& f_)
      : m_policy(policy_), m_func(f_) {}
#endif

  static constexpr index_type max_blocks = 65535;

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const {
    if (PolicyType::inner_direction == Iterate::Left) {
      const index_type temp0  = m_policy.m_tile_end[0];
      const index_type temp1  = m_policy.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = static_cast<index_type>(blockIdx.x) % numbl0;
      const index_type tile_id1 = static_cast<index_type>(blockIdx.x) / numbl0;
      const index_type thr_id0 =
          static_cast<index_type>(threadIdx.x) % m_policy.m_tile[0];
      const index_type thr_id1 =
          static_cast<index_type>(threadIdx.x) / m_policy.m_tile[0];

      for (index_type tile_id3 = static_cast<index_type>(blockIdx.z);
           tile_id3 < m_policy.m_tile_end[3]; tile_id3 += gridDim.z) {
        const index_type offset_3 =
            tile_id3 * m_policy.m_tile[3] +
            static_cast<index_type>(threadIdx.z) +
            static_cast<index_type>(m_policy.m_lower[3]);
        if (offset_3 < m_policy.m_upper[3] &&
            static_cast<index_type>(threadIdx.z) < m_policy.m_tile[3]) {
          for (index_type tile_id2 = static_cast<index_type>(blockIdx.y);
               tile_id2 < m_policy.m_tile_end[2]; tile_id2 += gridDim.y) {
            const index_type offset_2 =
                tile_id2 * m_policy.m_tile[2] +
                static_cast<index_type>(threadIdx.y) +
                static_cast<index_type>(m_policy.m_lower[2]);
            if (offset_2 < m_policy.m_upper[2] &&
                static_cast<index_type>(threadIdx.y) < m_policy.m_tile[2]) {
              for (index_type j = tile_id1; j < m_policy.m_tile_end[1];
                   j += numbl1) {
                const index_type offset_1 =
                    j * m_policy.m_tile[1] + thr_id1 +
                    static_cast<index_type>(m_policy.m_lower[1]);
                if (offset_1 < m_policy.m_upper[1] &&
                    thr_id1 < m_policy.m_tile[1]) {
                  for (index_type i = tile_id0; i < m_policy.m_tile_end[0];
                       i += numbl0) {
                    const index_type offset_0 =
                        i * m_policy.m_tile[0] + thr_id0 +
                        static_cast<index_type>(m_policy.m_lower[0]);
                    if (offset_0 < m_policy.m_upper[0] &&
                        thr_id0 < m_policy.m_tile[0]) {
                      Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1,
                                             offset_2, offset_3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      const index_type temp0  = m_policy.m_tile_end[0];
      const index_type temp1  = m_policy.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = static_cast<index_type>(blockIdx.x) / numbl1;
      const index_type tile_id1 = static_cast<index_type>(blockIdx.x) % numbl1;
      const index_type thr_id0 =
          static_cast<index_type>(threadIdx.x) / m_policy.m_tile[1];
      const index_type thr_id1 =
          static_cast<index_type>(threadIdx.x) % m_policy.m_tile[1];

      for (index_type i = tile_id0; i < m_policy.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_policy.m_tile[0] + thr_id0 +
            static_cast<index_type>(m_policy.m_lower[0]);
        if (offset_0 < m_policy.m_upper[0] && thr_id0 < m_policy.m_tile[0]) {
          for (index_type j = tile_id1; j < m_policy.m_tile_end[1];
               j += numbl1) {
            const index_type offset_1 =
                j * m_policy.m_tile[1] + thr_id1 +
                static_cast<index_type>(m_policy.m_lower[1]);
            if (offset_1 < m_policy.m_upper[1] &&
                thr_id1 < m_policy.m_tile[1]) {
              for (index_type tile_id2 = static_cast<index_type>(blockIdx.y);
                   tile_id2 < m_policy.m_tile_end[2]; tile_id2 += gridDim.y) {
                const index_type offset_2 =
                    tile_id2 * m_policy.m_tile[2] +
                    static_cast<index_type>(threadIdx.y) +
                    static_cast<index_type>(m_policy.m_lower[2]);
                if (offset_2 < m_policy.m_upper[2] &&
                    static_cast<index_type>(threadIdx.y) < m_policy.m_tile[2]) {
                  for (index_type tile_id3 =
                           static_cast<index_type>(blockIdx.z);
                       tile_id3 < m_policy.m_tile_end[3];
                       tile_id3 += gridDim.z) {
                    const index_type offset_3 =
                        tile_id3 * m_policy.m_tile[3] +
                        static_cast<index_type>(threadIdx.z) +
                        static_cast<index_type>(m_policy.m_lower[3]);
                    if (offset_3 < m_policy.m_upper[3] &&
                        static_cast<index_type>(threadIdx.z) <
                            m_policy.m_tile[3]) {
                      Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1,
                                             offset_2, offset_3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  // end exec_range

 private:
  const PolicyType& m_policy;
  const Functor& m_func;
#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif
};

// Rank 5
template <typename PolicyType, typename Functor, typename Tag>
struct DeviceIterateTile<5, PolicyType, Functor, Tag> {
  using index_type = typename PolicyType::index_type;

#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(
      const PolicyType& policy_, const Functor& f_,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_policy(policy_),
        m_func(f_),
        gridDim(gridDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(const PolicyType& policy_,
                                                const Functor& f_)
      : m_policy(policy_), m_func(f_) {}
#endif

  static constexpr index_type max_blocks = 65535;

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const {
    // LL
    if (PolicyType::inner_direction == Iterate::Left) {
      index_type temp0        = m_policy.m_tile_end[0];
      index_type temp1        = m_policy.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = static_cast<index_type>(blockIdx.x) % numbl0;
      const index_type tile_id1 = static_cast<index_type>(blockIdx.x) / numbl0;
      const index_type thr_id0 =
          static_cast<index_type>(threadIdx.x) % m_policy.m_tile[0];
      const index_type thr_id1 =
          static_cast<index_type>(threadIdx.x) / m_policy.m_tile[0];

      temp0                   = m_policy.m_tile_end[2];
      temp1                   = m_policy.m_tile_end[3];
      const index_type numbl2 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl3 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl2)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id2 = static_cast<index_type>(blockIdx.y) % numbl2;
      const index_type tile_id3 = static_cast<index_type>(blockIdx.y) / numbl2;
      const index_type thr_id2 =
          static_cast<index_type>(threadIdx.y) % m_policy.m_tile[2];
      const index_type thr_id3 =
          static_cast<index_type>(threadIdx.y) / m_policy.m_tile[2];

      for (index_type tile_id4 = static_cast<index_type>(blockIdx.z);
           tile_id4 < m_policy.m_tile_end[4]; tile_id4 += gridDim.z) {
        const index_type offset_4 =
            tile_id4 * m_policy.m_tile[4] +
            static_cast<index_type>(threadIdx.z) +
            static_cast<index_type>(m_policy.m_lower[4]);
        if (offset_4 < m_policy.m_upper[4] &&
            static_cast<index_type>(threadIdx.z) < m_policy.m_tile[4]) {
          for (index_type l = tile_id3; l < m_policy.m_tile_end[3];
               l += numbl3) {
            const index_type offset_3 =
                l * m_policy.m_tile[3] + thr_id3 +
                static_cast<index_type>(m_policy.m_lower[3]);
            if (offset_3 < m_policy.m_upper[3] &&
                thr_id3 < m_policy.m_tile[3]) {
              for (index_type k = tile_id2; k < m_policy.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_policy.m_tile[2] + thr_id2 +
                    static_cast<index_type>(m_policy.m_lower[2]);
                if (offset_2 < m_policy.m_upper[2] &&
                    thr_id2 < m_policy.m_tile[2]) {
                  for (index_type j = tile_id1; j < m_policy.m_tile_end[1];
                       j += numbl1) {
                    const index_type offset_1 =
                        j * m_policy.m_tile[1] + thr_id1 +
                        static_cast<index_type>(m_policy.m_lower[1]);
                    if (offset_1 < m_policy.m_upper[1] &&
                        thr_id1 < m_policy.m_tile[1]) {
                      for (index_type i = tile_id0; i < m_policy.m_tile_end[0];
                           i += numbl0) {
                        const index_type offset_0 =
                            i * m_policy.m_tile[0] + thr_id0 +
                            static_cast<index_type>(m_policy.m_lower[0]);
                        if (offset_0 < m_policy.m_upper[0] &&
                            thr_id0 < m_policy.m_tile[0]) {
                          Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1,
                                                 offset_2, offset_3, offset_4);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      index_type temp0        = m_policy.m_tile_end[0];
      index_type temp1        = m_policy.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = static_cast<index_type>(blockIdx.x) / numbl1;
      const index_type tile_id1 = static_cast<index_type>(blockIdx.x) % numbl1;
      const index_type thr_id0 =
          static_cast<index_type>(threadIdx.x) / m_policy.m_tile[1];
      const index_type thr_id1 =
          static_cast<index_type>(threadIdx.x) % m_policy.m_tile[1];

      temp0                   = m_policy.m_tile_end[2];
      temp1                   = m_policy.m_tile_end[3];
      const index_type numbl3 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl2 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl3)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id2 = static_cast<index_type>(blockIdx.y) / numbl3;
      const index_type tile_id3 = static_cast<index_type>(blockIdx.y) % numbl3;
      const index_type thr_id2 =
          static_cast<index_type>(threadIdx.y) / m_policy.m_tile[3];
      const index_type thr_id3 =
          static_cast<index_type>(threadIdx.y) % m_policy.m_tile[3];

      for (index_type i = tile_id0; i < m_policy.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_policy.m_tile[0] + thr_id0 +
            static_cast<index_type>(m_policy.m_lower[0]);
        if (offset_0 < m_policy.m_upper[0] && thr_id0 < m_policy.m_tile[0]) {
          for (index_type j = tile_id1; j < m_policy.m_tile_end[1];
               j += numbl1) {
            const index_type offset_1 =
                j * m_policy.m_tile[1] + thr_id1 +
                static_cast<index_type>(m_policy.m_lower[1]);
            if (offset_1 < m_policy.m_upper[1] &&
                thr_id1 < m_policy.m_tile[1]) {
              for (index_type k = tile_id2; k < m_policy.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_policy.m_tile[2] + thr_id2 +
                    static_cast<index_type>(m_policy.m_lower[2]);
                if (offset_2 < m_policy.m_upper[2] &&
                    thr_id2 < m_policy.m_tile[2]) {
                  for (index_type l = tile_id3; l < m_policy.m_tile_end[3];
                       l += numbl3) {
                    const index_type offset_3 =
                        l * m_policy.m_tile[3] + thr_id3 +
                        static_cast<index_type>(m_policy.m_lower[3]);
                    if (offset_3 < m_policy.m_upper[3] &&
                        thr_id3 < m_policy.m_tile[3]) {
                      for (index_type tile_id4 =
                               static_cast<index_type>(blockIdx.z);
                           tile_id4 < m_policy.m_tile_end[4];
                           tile_id4 += gridDim.z) {
                        const index_type offset_4 =
                            tile_id4 * m_policy.m_tile[4] +
                            static_cast<index_type>(threadIdx.z) +
                            static_cast<index_type>(m_policy.m_lower[4]);
                        if (offset_4 < m_policy.m_upper[4] &&
                            static_cast<index_type>(threadIdx.z) <
                                m_policy.m_tile[4]) {
                          Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1,
                                                 offset_2, offset_3, offset_4);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  // end exec_range

 private:
  const PolicyType& m_policy;
  const Functor& m_func;
#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif
};

// Rank 6
template <typename PolicyType, typename Functor, typename Tag>
struct DeviceIterateTile<6, PolicyType, Functor, Tag> {
  using index_type = typename PolicyType::index_type;

#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(
      const PolicyType& policy_, const Functor& f_,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_policy(policy_),
        m_func(f_),
        gridDim(gridDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(const PolicyType& policy_,
                                                const Functor& f_)
      : m_policy(policy_), m_func(f_) {}
#endif

  static constexpr index_type max_blocks = 65535;

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const {
    // LL
    if (PolicyType::inner_direction == Iterate::Left) {
      index_type temp0        = m_policy.m_tile_end[0];
      index_type temp1        = m_policy.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = static_cast<index_type>(blockIdx.x) % numbl0;
      const index_type tile_id1 = static_cast<index_type>(blockIdx.x) / numbl0;
      const index_type thr_id0 =
          static_cast<index_type>(threadIdx.x) % m_policy.m_tile[0];
      const index_type thr_id1 =
          static_cast<index_type>(threadIdx.x) / m_policy.m_tile[0];

      temp0                   = m_policy.m_tile_end[2];
      temp1                   = m_policy.m_tile_end[3];
      const index_type numbl2 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl3 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl2)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id2 = static_cast<index_type>(blockIdx.y) % numbl2;
      const index_type tile_id3 = static_cast<index_type>(blockIdx.y) / numbl2;
      const index_type thr_id2 =
          static_cast<index_type>(threadIdx.y) % m_policy.m_tile[2];
      const index_type thr_id3 =
          static_cast<index_type>(threadIdx.y) / m_policy.m_tile[2];

      temp0                   = m_policy.m_tile_end[4];
      temp1                   = m_policy.m_tile_end[5];
      const index_type numbl4 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl5 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl4)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id4 = static_cast<index_type>(blockIdx.z) % numbl4;
      const index_type tile_id5 = static_cast<index_type>(blockIdx.z) / numbl4;
      const index_type thr_id4 =
          static_cast<index_type>(threadIdx.z) % m_policy.m_tile[4];
      const index_type thr_id5 =
          static_cast<index_type>(threadIdx.z) / m_policy.m_tile[4];

      for (index_type n = tile_id5; n < m_policy.m_tile_end[5]; n += numbl5) {
        const index_type offset_5 =
            n * m_policy.m_tile[5] + thr_id5 +
            static_cast<index_type>(m_policy.m_lower[5]);
        if (offset_5 < m_policy.m_upper[5] && thr_id5 < m_policy.m_tile[5]) {
          for (index_type m = tile_id4; m < m_policy.m_tile_end[4];
               m += numbl4) {
            const index_type offset_4 =
                m * m_policy.m_tile[4] + thr_id4 +
                static_cast<index_type>(m_policy.m_lower[4]);
            if (offset_4 < m_policy.m_upper[4] &&
                thr_id4 < m_policy.m_tile[4]) {
              for (index_type l = tile_id3; l < m_policy.m_tile_end[3];
                   l += numbl3) {
                const index_type offset_3 =
                    l * m_policy.m_tile[3] + thr_id3 +
                    static_cast<index_type>(m_policy.m_lower[3]);
                if (offset_3 < m_policy.m_upper[3] &&
                    thr_id3 < m_policy.m_tile[3]) {
                  for (index_type k = tile_id2; k < m_policy.m_tile_end[2];
                       k += numbl2) {
                    const index_type offset_2 =
                        k * m_policy.m_tile[2] + thr_id2 +
                        static_cast<index_type>(m_policy.m_lower[2]);
                    if (offset_2 < m_policy.m_upper[2] &&
                        thr_id2 < m_policy.m_tile[2]) {
                      for (index_type j = tile_id1; j < m_policy.m_tile_end[1];
                           j += numbl1) {
                        const index_type offset_1 =
                            j * m_policy.m_tile[1] + thr_id1 +
                            static_cast<index_type>(m_policy.m_lower[1]);
                        if (offset_1 < m_policy.m_upper[1] &&
                            thr_id1 < m_policy.m_tile[1]) {
                          for (index_type i = tile_id0;
                               i < m_policy.m_tile_end[0]; i += numbl0) {
                            const index_type offset_0 =
                                i * m_policy.m_tile[0] + thr_id0 +
                                static_cast<index_type>(m_policy.m_lower[0]);
                            if (offset_0 < m_policy.m_upper[0] &&
                                thr_id0 < m_policy.m_tile[0]) {
                              Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1,
                                                     offset_2, offset_3,
                                                     offset_4, offset_5);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      index_type temp0        = m_policy.m_tile_end[0];
      index_type temp1        = m_policy.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = static_cast<index_type>(blockIdx.x) / numbl1;
      const index_type tile_id1 = static_cast<index_type>(blockIdx.x) % numbl1;
      const index_type thr_id0 =
          static_cast<index_type>(threadIdx.x) / m_policy.m_tile[1];
      const index_type thr_id1 =
          static_cast<index_type>(threadIdx.x) % m_policy.m_tile[1];

      temp0                   = m_policy.m_tile_end[2];
      temp1                   = m_policy.m_tile_end[3];
      const index_type numbl3 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl2 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl3)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id2 = static_cast<index_type>(blockIdx.y) / numbl3;
      const index_type tile_id3 = static_cast<index_type>(blockIdx.y) % numbl3;
      const index_type thr_id2 =
          static_cast<index_type>(threadIdx.y) / m_policy.m_tile[3];
      const index_type thr_id3 =
          static_cast<index_type>(threadIdx.y) % m_policy.m_tile[3];

      temp0                   = m_policy.m_tile_end[4];
      temp1                   = m_policy.m_tile_end[5];
      const index_type numbl5 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl4 =
          (temp0 * temp1 > max_blocks
               ? static_cast<index_type>(max_blocks / numbl5)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id4 = static_cast<index_type>(blockIdx.z) / numbl5;
      const index_type tile_id5 = static_cast<index_type>(blockIdx.z) % numbl5;
      const index_type thr_id4 =
          static_cast<index_type>(threadIdx.z) / m_policy.m_tile[5];
      const index_type thr_id5 =
          static_cast<index_type>(threadIdx.z) % m_policy.m_tile[5];

      for (index_type i = tile_id0; i < m_policy.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_policy.m_tile[0] + thr_id0 +
            static_cast<index_type>(m_policy.m_lower[0]);
        if (offset_0 < m_policy.m_upper[0] && thr_id0 < m_policy.m_tile[0]) {
          for (index_type j = tile_id1; j < m_policy.m_tile_end[1];
               j += numbl1) {
            const index_type offset_1 =
                j * m_policy.m_tile[1] + thr_id1 +
                static_cast<index_type>(m_policy.m_lower[1]);
            if (offset_1 < m_policy.m_upper[1] &&
                thr_id1 < m_policy.m_tile[1]) {
              for (index_type k = tile_id2; k < m_policy.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_policy.m_tile[2] + thr_id2 +
                    static_cast<index_type>(m_policy.m_lower[2]);
                if (offset_2 < m_policy.m_upper[2] &&
                    thr_id2 < m_policy.m_tile[2]) {
                  for (index_type l = tile_id3; l < m_policy.m_tile_end[3];
                       l += numbl3) {
                    const index_type offset_3 =
                        l * m_policy.m_tile[3] + thr_id3 +
                        static_cast<index_type>(m_policy.m_lower[3]);
                    if (offset_3 < m_policy.m_upper[3] &&
                        thr_id3 < m_policy.m_tile[3]) {
                      for (index_type m = tile_id4; m < m_policy.m_tile_end[4];
                           m += numbl4) {
                        const index_type offset_4 =
                            m * m_policy.m_tile[4] + thr_id4 +
                            static_cast<index_type>(m_policy.m_lower[4]);
                        if (offset_4 < m_policy.m_upper[4] &&
                            thr_id4 < m_policy.m_tile[4]) {
                          for (index_type n = tile_id5;
                               n < m_policy.m_tile_end[5]; n += numbl5) {
                            const index_type offset_5 =
                                n * m_policy.m_tile[5] + thr_id5 +
                                static_cast<index_type>(m_policy.m_lower[5]);
                            if (offset_5 < m_policy.m_upper[5] &&
                                thr_id5 < m_policy.m_tile[5]) {
                              Impl::_tag_invoke<Tag>(m_func, offset_0, offset_1,
                                                     offset_2, offset_3,
                                                     offset_4, offset_5);
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }  // end exec_range

 private:
  const PolicyType& m_policy;
  const Functor& m_func;
#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif
};

// ----------------------------------------------------------------------------------

namespace Reduce {

template <typename T>
using is_void = std::is_same<T, void>;

template <typename T>
struct is_array_type : std::false_type {
  using value_type = T;
};

template <typename T>
struct is_array_type<T*> : std::true_type {
  using value_type = T;
};

template <typename T>
struct is_array_type<T[]> : std::true_type {
  using value_type = T;
};

// ------------------------------------------------------------------ //

template <typename T>
using value_type_storage_t =
    typename std::conditional_t<is_array_type<T>::value, std::decay<T>,
                                std::add_lvalue_reference<T> >::type;

// ParallelReduce iteration pattern
// Scalar reductions

// num_blocks = min( num_tiles, max_num_blocks ); //i.e. determined by number of
// tiles and reduction algorithm constraints extract n-dim tile offsets (i.e.
// tile's global starting mulit-index) from the tileid = blockid using tile
// dimensions local indices within a tile extracted from (index_type)threadIdx_x
// using tile dims, constrained by blocksize combine tile and local id info for
// multi-dim global ids

// Pattern:
// Each block+thread is responsible for a tile+local_id combo (additional when
// striding by num_blocks)
// 1. create offset arrays
// 2. loop over number of tiles, striding by griddim (equal to num tiles, or max
// num blocks)
// 3. temps set for tile_idx and thrd_idx, which will be modified
// 4. if LL vs LR:
//      determine tile starting point offsets (multidim)
//      determine local index offsets (multidim)
//      concatentate tile offset + local offset for global multi-dim index
//    if offset withinin range bounds AND local offset within tile bounds, call
//    functor

template <int N, typename PolicyType, typename Functor, typename Tag,
          typename ValueType, typename Enable = void>
struct DeviceIterateTile {
  using index_type         = typename PolicyType::index_type;
  using value_type_storage = value_type_storage_t<ValueType>;

#ifdef KOKKOS_ENABLE_SYCL
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(
      const PolicyType& policy_, const Functor& f_, value_type_storage v_,
      const EmulateCUDADim3<index_type> gridDim_,
      const EmulateCUDADim3<index_type> blockIdx_,
      const EmulateCUDADim3<index_type> threadIdx_)
      : m_policy(policy_),
        m_func(f_),
        m_v(v_),
        gridDim(gridDim_),
        blockIdx(blockIdx_),
        threadIdx(threadIdx_) {}
#else
  KOKKOS_IMPL_DEVICE_FUNCTION DeviceIterateTile(const PolicyType& policy_,
                                                const Functor& f_,
                                                value_type_storage v_)
      : m_policy(policy_), m_func(f_), m_v(v_) {}
#endif

  KOKKOS_IMPL_DEVICE_FUNCTION
  void exec_range() const {
    if (static_cast<index_type>(blockIdx.x) < m_policy.m_num_tiles &&
        static_cast<index_type>(threadIdx.y) < m_policy.m_prod_tile_dims) {
      index_type m_offset[PolicyType::rank];  // tile starting global id offset
      index_type
          m_local_offset[PolicyType::rank];  // tile starting global id offset

      for (index_type tileidx = static_cast<index_type>(blockIdx.x);
           tileidx < m_policy.m_num_tiles; tileidx += gridDim.x) {
        index_type tile_idx =
            tileidx;  // temp because tile_idx will be modified while
                      // determining tile starting point offsets
        index_type thrd_idx = static_cast<index_type>(threadIdx.y);
        bool in_bounds      = true;

        // LL
        if (PolicyType::inner_direction == Iterate::Left) {
          for (int i = 0; i < PolicyType::rank; ++i) {
            m_offset[i] =
                (tile_idx % m_policy.m_tile_end[i]) * m_policy.m_tile[i] +
                m_policy.m_lower[i];
            tile_idx /= m_policy.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] = (thrd_idx % m_policy.m_tile[i]);
            thrd_idx /= m_policy.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if (!(m_offset[i] < m_policy.m_upper[i] &&
                  m_local_offset[i] < m_policy.m_tile[i])) {
              in_bounds = false;
            }
          }
          if (in_bounds) {
            Impl::_tag_invoke_array<Tag>(m_func, m_offset, m_v);
          }
        }
        // LR
        else {
          for (int i = PolicyType::rank - 1; i >= 0; --i) {
            m_offset[i] =
                (tile_idx % m_policy.m_tile_end[i]) * m_policy.m_tile[i] +
                m_policy.m_lower[i];
            tile_idx /= m_policy.m_tile_end[i];

            // tile-local indices identified with (index_type)threadIdx_y
            m_local_offset[i] =
                (thrd_idx %
                 m_policy.m_tile[i]);  // Move this to first computation,
                                       // add to m_offset right away
            thrd_idx /= m_policy.m_tile[i];

            m_offset[i] += m_local_offset[i];
            if (!(m_offset[i] < m_policy.m_upper[i] &&
                  m_local_offset[i] < m_policy.m_tile[i])) {
              in_bounds = false;
            }
          }
          if (in_bounds) {
            Impl::_tag_invoke_array<Tag>(m_func, m_offset, m_v);
          }
        }
      }
    }
  }  // end exec_range

 private:
  const PolicyType& m_policy;
  const Functor& m_func;
  value_type_storage m_v;
#ifdef KOKKOS_ENABLE_SYCL
  const EmulateCUDADim3<index_type> gridDim;
  const EmulateCUDADim3<index_type> blockIdx;
  const EmulateCUDADim3<index_type> threadIdx;
#endif
};

}  // namespace Reduce
}  // namespace Impl
}  // namespace Kokkos
#endif
