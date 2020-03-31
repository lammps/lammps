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

#ifndef KOKKOS_CUDA_EXP_ITERATE_TILE_HPP
#define KOKKOS_CUDA_EXP_ITERATE_TILE_HPP

#include <Kokkos_Macros.hpp>
#if defined(__CUDACC__) && defined(KOKKOS_ENABLE_CUDA)

#include <iostream>
#include <algorithm>
#include <cstdio>

#include <utility>

//#include<Cuda/Kokkos_CudaExec.hpp>
// Including the file above, leads to following type of errors:
// /home/ndellin/kokkos/core/src/Cuda/Kokkos_CudaExec.hpp(84): error: incomplete
// type is not allowed As a result, recreate cuda_parallel_launch and associated
// code

#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <typeinfo>
#endif

namespace Kokkos {
namespace Impl {

// ------------------------------------------------------------------ //

template <class DriverType>
__global__ static void cuda_parallel_launch(const DriverType driver) {
  driver();
}

template <class DriverType>
struct CudaLaunch {
  inline CudaLaunch(const DriverType& driver, const dim3& grid,
                    const dim3& block) {
    cuda_parallel_launch<DriverType><<<grid, block>>>(driver);
  }
};

// ------------------------------------------------------------------ //
template <int N, typename RP, typename Functor, typename Tag>
struct apply_impl;

// Rank 2
// Specializations for void tag type
template <typename RP, typename Functor>
struct apply_impl<2, RP, Functor, void> {
  using index_type = typename RP::index_type;

  __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  inline __device__ void exec_range() const {
    // LL
    if (RP::inner_direction == RP::Left) {
      for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
           tile_id1 += gridDim.y) {
        const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                    (index_type)threadIdx.y +
                                    (index_type)m_rp.m_lower[1];
        if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
          for (index_type tile_id0 = blockIdx.x; tile_id0 < m_rp.m_tile_end[0];
               tile_id0 += gridDim.x) {
            const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                        (index_type)threadIdx.x +
                                        (index_type)m_rp.m_lower[0];
            if (offset_0 < m_rp.m_upper[0] && threadIdx.x < m_rp.m_tile[0]) {
              m_func(offset_0, offset_1);
            }
          }
        }
      }
    }
    // LR
    else {
      for (index_type tile_id0 = blockIdx.x; tile_id0 < m_rp.m_tile_end[0];
           tile_id0 += gridDim.x) {
        const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                    (index_type)threadIdx.x +
                                    (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && threadIdx.x < m_rp.m_tile[0]) {
          for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
               tile_id1 += gridDim.y) {
            const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
              m_func(offset_0, offset_1);
            }
          }
        }
      }
    }

  }  // end exec_range

 private:
  const RP& m_rp;
  const Functor& m_func;
};

// Specializations for tag type
template <typename RP, typename Functor, typename Tag>
struct apply_impl<2, RP, Functor, Tag> {
  using index_type = typename RP::index_type;

  inline __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  inline __device__ void exec_range() const {
    if (RP::inner_direction == RP::Left) {
      // Loop over size maxnumblocks until full range covered
      for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
           tile_id1 += gridDim.y) {
        const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                    (index_type)threadIdx.y +
                                    (index_type)m_rp.m_lower[1];
        if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
          for (index_type tile_id0 = blockIdx.x; tile_id0 < m_rp.m_tile_end[0];
               tile_id0 += gridDim.x) {
            const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                        (index_type)threadIdx.x +
                                        (index_type)m_rp.m_lower[0];
            if (offset_0 < m_rp.m_upper[0] && threadIdx.x < m_rp.m_tile[0]) {
              m_func(Tag(), offset_0, offset_1);
            }
          }
        }
      }
    } else {
      for (index_type tile_id0 = blockIdx.x; tile_id0 < m_rp.m_tile_end[0];
           tile_id0 += gridDim.x) {
        const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                    (index_type)threadIdx.x +
                                    (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && threadIdx.x < m_rp.m_tile[0]) {
          for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
               tile_id1 += gridDim.y) {
            const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
              m_func(Tag(), offset_0, offset_1);
            }
          }
        }
      }
    }

  }  // end exec_range

 private:
  const RP& m_rp;
  const Functor& m_func;
};

// Rank 3
// Specializations for void tag type
template <typename RP, typename Functor>
struct apply_impl<3, RP, Functor, void> {
  using index_type = typename RP::index_type;

  __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  inline __device__ void exec_range() const {
    // LL
    if (RP::inner_direction == RP::Left) {
      for (index_type tile_id2 = blockIdx.z; tile_id2 < m_rp.m_tile_end[2];
           tile_id2 += gridDim.z) {
        const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                    (index_type)threadIdx.z +
                                    (index_type)m_rp.m_lower[2];
        if (offset_2 < m_rp.m_upper[2] && threadIdx.z < m_rp.m_tile[2]) {
          for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
               tile_id1 += gridDim.y) {
            const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
              for (index_type tile_id0 = blockIdx.x;
                   tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim.x) {
                const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                            (index_type)threadIdx.x +
                                            (index_type)m_rp.m_lower[0];
                if (offset_0 < m_rp.m_upper[0] &&
                    threadIdx.x < m_rp.m_tile[0]) {
                  m_func(offset_0, offset_1, offset_2);
                }
              }
            }
          }
        }
      }
    }
    // LR
    else {
      for (index_type tile_id0 = blockIdx.x; tile_id0 < m_rp.m_tile_end[0];
           tile_id0 += gridDim.x) {
        const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                    (index_type)threadIdx.x +
                                    (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && threadIdx.x < m_rp.m_tile[0]) {
          for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
               tile_id1 += gridDim.y) {
            const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
              for (index_type tile_id2 = blockIdx.z;
                   tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim.z) {
                const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                            (index_type)threadIdx.z +
                                            (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] &&
                    threadIdx.z < m_rp.m_tile[2]) {
                  m_func(offset_0, offset_1, offset_2);
                }
              }
            }
          }
        }
      }
    }

  }  // end exec_range

 private:
  const RP& m_rp;
  const Functor& m_func;
};

// Specializations for void tag type
template <typename RP, typename Functor, typename Tag>
struct apply_impl<3, RP, Functor, Tag> {
  using index_type = typename RP::index_type;

  inline __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  inline __device__ void exec_range() const {
    if (RP::inner_direction == RP::Left) {
      for (index_type tile_id2 = blockIdx.z; tile_id2 < m_rp.m_tile_end[2];
           tile_id2 += gridDim.z) {
        const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                    (index_type)threadIdx.z +
                                    (index_type)m_rp.m_lower[2];
        if (offset_2 < m_rp.m_upper[2] && threadIdx.z < m_rp.m_tile[2]) {
          for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
               tile_id1 += gridDim.y) {
            const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
              for (index_type tile_id0 = blockIdx.x;
                   tile_id0 < m_rp.m_tile_end[0]; tile_id0 += gridDim.x) {
                const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                            (index_type)threadIdx.x +
                                            (index_type)m_rp.m_lower[0];
                if (offset_0 < m_rp.m_upper[0] &&
                    threadIdx.x < m_rp.m_tile[0]) {
                  m_func(Tag(), offset_0, offset_1, offset_2);
                }
              }
            }
          }
        }
      }
    } else {
      for (index_type tile_id0 = blockIdx.x; tile_id0 < m_rp.m_tile_end[0];
           tile_id0 += gridDim.x) {
        const index_type offset_0 = tile_id0 * m_rp.m_tile[0] +
                                    (index_type)threadIdx.x +
                                    (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && threadIdx.x < m_rp.m_tile[0]) {
          for (index_type tile_id1 = blockIdx.y; tile_id1 < m_rp.m_tile_end[1];
               tile_id1 += gridDim.y) {
            const index_type offset_1 = tile_id1 * m_rp.m_tile[1] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && threadIdx.y < m_rp.m_tile[1]) {
              for (index_type tile_id2 = blockIdx.z;
                   tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim.z) {
                const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                            (index_type)threadIdx.z +
                                            (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] &&
                    threadIdx.z < m_rp.m_tile[2]) {
                  m_func(Tag(), offset_0, offset_1, offset_2);
                }
              }
            }
          }
        }
      }
    }

  }  // end exec_range

 private:
  const RP& m_rp;
  const Functor& m_func;
};

// Rank 4
// Specializations for void tag type
template <typename RP, typename Functor>
struct apply_impl<4, RP, Functor, void> {
  using index_type = typename RP::index_type;

  __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  static constexpr index_type max_blocks = 65535;

  inline __device__ void exec_range() const {
    // LL
    if (RP::inner_direction == RP::Left) {
      const index_type temp0  = m_rp.m_tile_end[0];
      const index_type temp1  = m_rp.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = blockIdx.x % numbl0;
      const index_type tile_id1 = blockIdx.x / numbl0;
      const index_type thr_id0  = threadIdx.x % m_rp.m_tile[0];
      const index_type thr_id1  = threadIdx.x / m_rp.m_tile[0];

      for (index_type tile_id3 = blockIdx.z; tile_id3 < m_rp.m_tile_end[3];
           tile_id3 += gridDim.z) {
        const index_type offset_3 = tile_id3 * m_rp.m_tile[3] +
                                    (index_type)threadIdx.z +
                                    (index_type)m_rp.m_lower[3];
        if (offset_3 < m_rp.m_upper[3] && threadIdx.z < m_rp.m_tile[3]) {
          for (index_type tile_id2 = blockIdx.y; tile_id2 < m_rp.m_tile_end[2];
               tile_id2 += gridDim.y) {
            const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[2];
            if (offset_2 < m_rp.m_upper[2] && threadIdx.y < m_rp.m_tile[2]) {
              for (index_type j = tile_id1; j < m_rp.m_tile_end[1];
                   j += numbl1) {
                const index_type offset_1 =
                    j * m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
                  for (index_type i = tile_id0; i < m_rp.m_tile_end[0];
                       i += numbl0) {
                    const index_type offset_0 = i * m_rp.m_tile[0] + thr_id0 +
                                                (index_type)m_rp.m_lower[0];
                    if (offset_0 < m_rp.m_upper[0] &&
                        thr_id0 < m_rp.m_tile[0]) {
                      m_func(offset_0, offset_1, offset_2, offset_3);
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
      const index_type temp0  = m_rp.m_tile_end[0];
      const index_type temp1  = m_rp.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = blockIdx.x / numbl1;
      const index_type tile_id1 = blockIdx.x % numbl1;
      const index_type thr_id0  = threadIdx.x / m_rp.m_tile[1];
      const index_type thr_id1  = threadIdx.x % m_rp.m_tile[1];

      for (index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0]) {
          for (index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1) {
            const index_type offset_1 =
                j * m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
              for (index_type tile_id2 = blockIdx.y;
                   tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim.y) {
                const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                            (index_type)threadIdx.y +
                                            (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] &&
                    threadIdx.y < m_rp.m_tile[2]) {
                  for (index_type tile_id3 = blockIdx.z;
                       tile_id3 < m_rp.m_tile_end[3]; tile_id3 += gridDim.z) {
                    const index_type offset_3 = tile_id3 * m_rp.m_tile[3] +
                                                (index_type)threadIdx.z +
                                                (index_type)m_rp.m_lower[3];
                    if (offset_3 < m_rp.m_upper[3] &&
                        threadIdx.z < m_rp.m_tile[3]) {
                      m_func(offset_0, offset_1, offset_2, offset_3);
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
  const RP& m_rp;
  const Functor& m_func;
};

// Specializations for void tag type
template <typename RP, typename Functor, typename Tag>
struct apply_impl<4, RP, Functor, Tag> {
  using index_type = typename RP::index_type;

  inline __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  static constexpr index_type max_blocks = 65535;

  inline __device__ void exec_range() const {
    if (RP::inner_direction == RP::Left) {
      const index_type temp0  = m_rp.m_tile_end[0];
      const index_type temp1  = m_rp.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = blockIdx.x % numbl0;
      const index_type tile_id1 = blockIdx.x / numbl0;
      const index_type thr_id0  = threadIdx.x % m_rp.m_tile[0];
      const index_type thr_id1  = threadIdx.x / m_rp.m_tile[0];

      for (index_type tile_id3 = blockIdx.z; tile_id3 < m_rp.m_tile_end[3];
           tile_id3 += gridDim.z) {
        const index_type offset_3 = tile_id3 * m_rp.m_tile[3] +
                                    (index_type)threadIdx.z +
                                    (index_type)m_rp.m_lower[3];
        if (offset_3 < m_rp.m_upper[3] && threadIdx.z < m_rp.m_tile[3]) {
          for (index_type tile_id2 = blockIdx.y; tile_id2 < m_rp.m_tile_end[2];
               tile_id2 += gridDim.y) {
            const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                        (index_type)threadIdx.y +
                                        (index_type)m_rp.m_lower[2];
            if (offset_2 < m_rp.m_upper[2] && threadIdx.y < m_rp.m_tile[2]) {
              for (index_type j = tile_id1; j < m_rp.m_tile_end[1];
                   j += numbl1) {
                const index_type offset_1 =
                    j * m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
                if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
                  for (index_type i = tile_id0; i < m_rp.m_tile_end[0];
                       i += numbl0) {
                    const index_type offset_0 = i * m_rp.m_tile[0] + thr_id0 +
                                                (index_type)m_rp.m_lower[0];
                    if (offset_0 < m_rp.m_upper[0] &&
                        thr_id0 < m_rp.m_tile[0]) {
                      m_func(Tag(), offset_0, offset_1, offset_2, offset_3);
                    }
                  }
                }
              }
            }
          }
        }
      }
    } else {
      const index_type temp0  = m_rp.m_tile_end[0];
      const index_type temp1  = m_rp.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = blockIdx.x / numbl1;
      const index_type tile_id1 = blockIdx.x % numbl1;
      const index_type thr_id0  = threadIdx.x / m_rp.m_tile[1];
      const index_type thr_id1  = threadIdx.x % m_rp.m_tile[1];

      for (index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0]) {
          for (index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1) {
            const index_type offset_1 = tile_id1 * m_rp.m_tile[1] + thr_id1 +
                                        (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
              for (index_type tile_id2 = blockIdx.y;
                   tile_id2 < m_rp.m_tile_end[2]; tile_id2 += gridDim.y) {
                const index_type offset_2 = tile_id2 * m_rp.m_tile[2] +
                                            (index_type)threadIdx.y +
                                            (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] &&
                    threadIdx.y < m_rp.m_tile[2]) {
                  for (index_type tile_id3 = blockIdx.z;
                       tile_id3 < m_rp.m_tile_end[3]; tile_id3 += gridDim.z) {
                    const index_type offset_3 = tile_id3 * m_rp.m_tile[3] +
                                                (index_type)threadIdx.z +
                                                (index_type)m_rp.m_lower[3];
                    if (offset_3 < m_rp.m_upper[3] &&
                        threadIdx.z < m_rp.m_tile[3]) {
                      m_func(Tag(), offset_0, offset_1, offset_2, offset_3);
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
  const RP& m_rp;
  const Functor& m_func;
};

// Rank 5
// Specializations for void tag type
template <typename RP, typename Functor>
struct apply_impl<5, RP, Functor, void> {
  using index_type = typename RP::index_type;

  __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  static constexpr index_type max_blocks = 65535;

  inline __device__ void exec_range() const {
    // LL
    if (RP::inner_direction == RP::Left) {
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = blockIdx.x % numbl0;
      const index_type tile_id1 = blockIdx.x / numbl0;
      const index_type thr_id0  = threadIdx.x % m_rp.m_tile[0];
      const index_type thr_id1  = threadIdx.x / m_rp.m_tile[0];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl2 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl3 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl2)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id2 = blockIdx.y % numbl2;
      const index_type tile_id3 = blockIdx.y / numbl2;
      const index_type thr_id2  = threadIdx.y % m_rp.m_tile[2];
      const index_type thr_id3  = threadIdx.y / m_rp.m_tile[2];

      for (index_type tile_id4 = blockIdx.z; tile_id4 < m_rp.m_tile_end[4];
           tile_id4 += gridDim.z) {
        const index_type offset_4 = tile_id4 * m_rp.m_tile[4] +
                                    (index_type)threadIdx.z +
                                    (index_type)m_rp.m_lower[4];
        if (offset_4 < m_rp.m_upper[4] && threadIdx.z < m_rp.m_tile[4]) {
          for (index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3) {
            const index_type offset_3 =
                l * m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
            if (offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3]) {
              for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2]) {
                  for (index_type j = tile_id1; j < m_rp.m_tile_end[1];
                       j += numbl1) {
                    const index_type offset_1 = j * m_rp.m_tile[1] + thr_id1 +
                                                (index_type)m_rp.m_lower[1];
                    if (offset_1 < m_rp.m_upper[1] &&
                        thr_id1 < m_rp.m_tile[1]) {
                      for (index_type i = tile_id0; i < m_rp.m_tile_end[0];
                           i += numbl0) {
                        const index_type offset_0 = i * m_rp.m_tile[0] +
                                                    thr_id0 +
                                                    (index_type)m_rp.m_lower[0];
                        if (offset_0 < m_rp.m_upper[0] &&
                            thr_id0 < m_rp.m_tile[0]) {
                          m_func(offset_0, offset_1, offset_2, offset_3,
                                 offset_4);
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
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = blockIdx.x / numbl1;
      const index_type tile_id1 = blockIdx.x % numbl1;
      const index_type thr_id0  = threadIdx.x / m_rp.m_tile[1];
      const index_type thr_id1  = threadIdx.x % m_rp.m_tile[1];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl3 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl2 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl3)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id2 = blockIdx.y / numbl3;
      const index_type tile_id3 = blockIdx.y % numbl3;
      const index_type thr_id2  = threadIdx.y / m_rp.m_tile[3];
      const index_type thr_id3  = threadIdx.y % m_rp.m_tile[3];

      for (index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0]) {
          for (index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1) {
            const index_type offset_1 =
                j * m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
              for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2]) {
                  for (index_type l = tile_id3; l < m_rp.m_tile_end[3];
                       l += numbl3) {
                    const index_type offset_3 = l * m_rp.m_tile[3] + thr_id3 +
                                                (index_type)m_rp.m_lower[3];
                    if (offset_3 < m_rp.m_upper[3] &&
                        thr_id3 < m_rp.m_tile[3]) {
                      for (index_type tile_id4 = blockIdx.z;
                           tile_id4 < m_rp.m_tile_end[4];
                           tile_id4 += gridDim.z) {
                        const index_type offset_4 = tile_id4 * m_rp.m_tile[4] +
                                                    (index_type)threadIdx.z +
                                                    (index_type)m_rp.m_lower[4];
                        if (offset_4 < m_rp.m_upper[4] &&
                            threadIdx.z < m_rp.m_tile[4]) {
                          m_func(offset_0, offset_1, offset_2, offset_3,
                                 offset_4);
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
  const RP& m_rp;
  const Functor& m_func;
};

// Specializations for tag type
template <typename RP, typename Functor, typename Tag>
struct apply_impl<5, RP, Functor, Tag> {
  using index_type = typename RP::index_type;

  __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  static constexpr index_type max_blocks = 65535;

  inline __device__ void exec_range() const {
    // LL
    if (RP::inner_direction == RP::Left) {
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = blockIdx.x % numbl0;
      const index_type tile_id1 = blockIdx.x / numbl0;
      const index_type thr_id0  = threadIdx.x % m_rp.m_tile[0];
      const index_type thr_id1  = threadIdx.x / m_rp.m_tile[0];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl2 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl3 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl2)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id2 = blockIdx.y % numbl2;
      const index_type tile_id3 = blockIdx.y / numbl2;
      const index_type thr_id2  = threadIdx.y % m_rp.m_tile[2];
      const index_type thr_id3  = threadIdx.y / m_rp.m_tile[2];

      for (index_type tile_id4 = blockIdx.z; tile_id4 < m_rp.m_tile_end[4];
           tile_id4 += gridDim.z) {
        const index_type offset_4 = tile_id4 * m_rp.m_tile[4] +
                                    (index_type)threadIdx.z +
                                    (index_type)m_rp.m_lower[4];
        if (offset_4 < m_rp.m_upper[4] && threadIdx.z < m_rp.m_tile[4]) {
          for (index_type l = tile_id3; l < m_rp.m_tile_end[3]; l += numbl3) {
            const index_type offset_3 =
                l * m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
            if (offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3]) {
              for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2]) {
                  for (index_type j = tile_id1; j < m_rp.m_tile_end[1];
                       j += numbl1) {
                    const index_type offset_1 = j * m_rp.m_tile[1] + thr_id1 +
                                                (index_type)m_rp.m_lower[1];
                    if (offset_1 < m_rp.m_upper[1] &&
                        thr_id1 < m_rp.m_tile[1]) {
                      for (index_type i = tile_id0; i < m_rp.m_tile_end[0];
                           i += numbl0) {
                        const index_type offset_0 = i * m_rp.m_tile[0] +
                                                    thr_id0 +
                                                    (index_type)m_rp.m_lower[0];
                        if (offset_0 < m_rp.m_upper[0] &&
                            thr_id0 < m_rp.m_tile[0]) {
                          m_func(Tag(), offset_0, offset_1, offset_2, offset_3,
                                 offset_4);
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
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = blockIdx.x / numbl1;
      const index_type tile_id1 = blockIdx.x % numbl1;
      const index_type thr_id0  = threadIdx.x / m_rp.m_tile[1];
      const index_type thr_id1  = threadIdx.x % m_rp.m_tile[1];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl3 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl2 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl3)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id2 = blockIdx.y / numbl3;
      const index_type tile_id3 = blockIdx.y % numbl3;
      const index_type thr_id2  = threadIdx.y / m_rp.m_tile[3];
      const index_type thr_id3  = threadIdx.y % m_rp.m_tile[3];

      for (index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0]) {
          for (index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1) {
            const index_type offset_1 =
                j * m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
              for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2]) {
                  for (index_type l = tile_id3; l < m_rp.m_tile_end[3];
                       l += numbl3) {
                    const index_type offset_3 = l * m_rp.m_tile[3] + thr_id3 +
                                                (index_type)m_rp.m_lower[3];
                    if (offset_3 < m_rp.m_upper[3] &&
                        thr_id3 < m_rp.m_tile[3]) {
                      for (index_type tile_id4 = blockIdx.z;
                           tile_id4 < m_rp.m_tile_end[4];
                           tile_id4 += gridDim.z) {
                        const index_type offset_4 = tile_id4 * m_rp.m_tile[4] +
                                                    (index_type)threadIdx.z +
                                                    (index_type)m_rp.m_lower[4];
                        if (offset_4 < m_rp.m_upper[4] &&
                            threadIdx.z < m_rp.m_tile[4]) {
                          m_func(Tag(), offset_0, offset_1, offset_2, offset_3,
                                 offset_4);
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
  const RP& m_rp;
  const Functor& m_func;
};

// Rank 6
// Specializations for void tag type
template <typename RP, typename Functor>
struct apply_impl<6, RP, Functor, void> {
  using index_type = typename RP::index_type;

  __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  static constexpr index_type max_blocks = 65535;

  inline __device__ void exec_range() const {
    // LL
    if (RP::inner_direction == RP::Left) {
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = blockIdx.x % numbl0;
      const index_type tile_id1 = blockIdx.x / numbl0;
      const index_type thr_id0  = threadIdx.x % m_rp.m_tile[0];
      const index_type thr_id1  = threadIdx.x / m_rp.m_tile[0];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl2 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl3 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl2)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id2 = blockIdx.y % numbl2;
      const index_type tile_id3 = blockIdx.y / numbl2;
      const index_type thr_id2  = threadIdx.y % m_rp.m_tile[2];
      const index_type thr_id3  = threadIdx.y / m_rp.m_tile[2];

      temp0                   = m_rp.m_tile_end[4];
      temp1                   = m_rp.m_tile_end[5];
      const index_type numbl4 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl5 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl4)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id4 = blockIdx.z % numbl4;
      const index_type tile_id5 = blockIdx.z / numbl4;
      const index_type thr_id4  = threadIdx.z % m_rp.m_tile[4];
      const index_type thr_id5  = threadIdx.z / m_rp.m_tile[4];

      for (index_type n = tile_id5; n < m_rp.m_tile_end[5]; n += numbl5) {
        const index_type offset_5 =
            n * m_rp.m_tile[5] + thr_id5 + (index_type)m_rp.m_lower[5];
        if (offset_5 < m_rp.m_upper[5] && thr_id5 < m_rp.m_tile[5]) {
          for (index_type m = tile_id4; m < m_rp.m_tile_end[4]; m += numbl4) {
            const index_type offset_4 =
                m * m_rp.m_tile[4] + thr_id4 + (index_type)m_rp.m_lower[4];
            if (offset_4 < m_rp.m_upper[4] && thr_id4 < m_rp.m_tile[4]) {
              for (index_type l = tile_id3; l < m_rp.m_tile_end[3];
                   l += numbl3) {
                const index_type offset_3 =
                    l * m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                if (offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3]) {
                  for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                       k += numbl2) {
                    const index_type offset_2 = k * m_rp.m_tile[2] + thr_id2 +
                                                (index_type)m_rp.m_lower[2];
                    if (offset_2 < m_rp.m_upper[2] &&
                        thr_id2 < m_rp.m_tile[2]) {
                      for (index_type j = tile_id1; j < m_rp.m_tile_end[1];
                           j += numbl1) {
                        const index_type offset_1 = j * m_rp.m_tile[1] +
                                                    thr_id1 +
                                                    (index_type)m_rp.m_lower[1];
                        if (offset_1 < m_rp.m_upper[1] &&
                            thr_id1 < m_rp.m_tile[1]) {
                          for (index_type i = tile_id0; i < m_rp.m_tile_end[0];
                               i += numbl0) {
                            const index_type offset_0 =
                                i * m_rp.m_tile[0] + thr_id0 +
                                (index_type)m_rp.m_lower[0];
                            if (offset_0 < m_rp.m_upper[0] &&
                                thr_id0 < m_rp.m_tile[0]) {
                              m_func(offset_0, offset_1, offset_2, offset_3,
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
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = blockIdx.x / numbl1;
      const index_type tile_id1 = blockIdx.x % numbl1;
      const index_type thr_id0  = threadIdx.x / m_rp.m_tile[1];
      const index_type thr_id1  = threadIdx.x % m_rp.m_tile[1];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl3 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl2 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl3)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id2 = blockIdx.y / numbl3;
      const index_type tile_id3 = blockIdx.y % numbl3;
      const index_type thr_id2  = threadIdx.y / m_rp.m_tile[3];
      const index_type thr_id3  = threadIdx.y % m_rp.m_tile[3];

      temp0                   = m_rp.m_tile_end[4];
      temp1                   = m_rp.m_tile_end[5];
      const index_type numbl5 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl4 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl5)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id4 = blockIdx.z / numbl5;
      const index_type tile_id5 = blockIdx.z % numbl5;
      const index_type thr_id4  = threadIdx.z / m_rp.m_tile[5];
      const index_type thr_id5  = threadIdx.z % m_rp.m_tile[5];

      for (index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0]) {
          for (index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1) {
            const index_type offset_1 =
                j * m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
              for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2]) {
                  for (index_type l = tile_id3; l < m_rp.m_tile_end[3];
                       l += numbl3) {
                    const index_type offset_3 = l * m_rp.m_tile[3] + thr_id3 +
                                                (index_type)m_rp.m_lower[3];
                    if (offset_3 < m_rp.m_upper[3] &&
                        thr_id3 < m_rp.m_tile[3]) {
                      for (index_type m = tile_id4; m < m_rp.m_tile_end[4];
                           m += numbl4) {
                        const index_type offset_4 = m * m_rp.m_tile[4] +
                                                    thr_id4 +
                                                    (index_type)m_rp.m_lower[4];
                        if (offset_4 < m_rp.m_upper[4] &&
                            thr_id4 < m_rp.m_tile[4]) {
                          for (index_type n = tile_id5; n < m_rp.m_tile_end[5];
                               n += numbl5) {
                            const index_type offset_5 =
                                n * m_rp.m_tile[5] + thr_id5 +
                                (index_type)m_rp.m_lower[5];
                            if (offset_5 < m_rp.m_upper[5] &&
                                thr_id5 < m_rp.m_tile[5]) {
                              m_func(offset_0, offset_1, offset_2, offset_3,
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
  const RP& m_rp;
  const Functor& m_func;
};

// Specializations for tag type
template <typename RP, typename Functor, typename Tag>
struct apply_impl<6, RP, Functor, Tag> {
  using index_type = typename RP::index_type;

  __device__ apply_impl(const RP& rp_, const Functor& f_)
      : m_rp(rp_), m_func(f_) {}

  static constexpr index_type max_blocks = 65535;

  inline __device__ void exec_range() const {
    // LL
    if (RP::inner_direction == RP::Left) {
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl0 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl1 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl0)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id0 = blockIdx.x % numbl0;
      const index_type tile_id1 = blockIdx.x / numbl0;
      const index_type thr_id0  = threadIdx.x % m_rp.m_tile[0];
      const index_type thr_id1  = threadIdx.x / m_rp.m_tile[0];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl2 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl3 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl2)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id2 = blockIdx.y % numbl2;
      const index_type tile_id3 = blockIdx.y / numbl2;
      const index_type thr_id2  = threadIdx.y % m_rp.m_tile[2];
      const index_type thr_id3  = threadIdx.y / m_rp.m_tile[2];

      temp0                   = m_rp.m_tile_end[4];
      temp1                   = m_rp.m_tile_end[5];
      const index_type numbl4 = (temp0 <= max_blocks ? temp0 : max_blocks);
      const index_type numbl5 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl4)
               : (temp1 <= max_blocks ? temp1 : max_blocks));

      const index_type tile_id4 = blockIdx.z % numbl4;
      const index_type tile_id5 = blockIdx.z / numbl4;
      const index_type thr_id4  = threadIdx.z % m_rp.m_tile[4];
      const index_type thr_id5  = threadIdx.z / m_rp.m_tile[4];

      for (index_type n = tile_id5; n < m_rp.m_tile_end[5]; n += numbl5) {
        const index_type offset_5 =
            n * m_rp.m_tile[5] + thr_id5 + (index_type)m_rp.m_lower[5];
        if (offset_5 < m_rp.m_upper[5] && thr_id5 < m_rp.m_tile[5]) {
          for (index_type m = tile_id4; m < m_rp.m_tile_end[4]; m += numbl4) {
            const index_type offset_4 =
                m * m_rp.m_tile[4] + thr_id4 + (index_type)m_rp.m_lower[4];
            if (offset_4 < m_rp.m_upper[4] && thr_id4 < m_rp.m_tile[4]) {
              for (index_type l = tile_id3; l < m_rp.m_tile_end[3];
                   l += numbl3) {
                const index_type offset_3 =
                    l * m_rp.m_tile[3] + thr_id3 + (index_type)m_rp.m_lower[3];
                if (offset_3 < m_rp.m_upper[3] && thr_id3 < m_rp.m_tile[3]) {
                  for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                       k += numbl2) {
                    const index_type offset_2 = k * m_rp.m_tile[2] + thr_id2 +
                                                (index_type)m_rp.m_lower[2];
                    if (offset_2 < m_rp.m_upper[2] &&
                        thr_id2 < m_rp.m_tile[2]) {
                      for (index_type j = tile_id1; j < m_rp.m_tile_end[1];
                           j += numbl1) {
                        const index_type offset_1 = j * m_rp.m_tile[1] +
                                                    thr_id1 +
                                                    (index_type)m_rp.m_lower[1];
                        if (offset_1 < m_rp.m_upper[1] &&
                            thr_id1 < m_rp.m_tile[1]) {
                          for (index_type i = tile_id0; i < m_rp.m_tile_end[0];
                               i += numbl0) {
                            const index_type offset_0 =
                                i * m_rp.m_tile[0] + thr_id0 +
                                (index_type)m_rp.m_lower[0];
                            if (offset_0 < m_rp.m_upper[0] &&
                                thr_id0 < m_rp.m_tile[0]) {
                              m_func(Tag(), offset_0, offset_1, offset_2,
                                     offset_3, offset_4, offset_5);
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
      index_type temp0        = m_rp.m_tile_end[0];
      index_type temp1        = m_rp.m_tile_end[1];
      const index_type numbl1 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl0 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl1)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id0 = blockIdx.x / numbl1;
      const index_type tile_id1 = blockIdx.x % numbl1;
      const index_type thr_id0  = threadIdx.x / m_rp.m_tile[1];
      const index_type thr_id1  = threadIdx.x % m_rp.m_tile[1];

      temp0                   = m_rp.m_tile_end[2];
      temp1                   = m_rp.m_tile_end[3];
      const index_type numbl3 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl2 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl3)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id2 = blockIdx.y / numbl3;
      const index_type tile_id3 = blockIdx.y % numbl3;
      const index_type thr_id2  = threadIdx.y / m_rp.m_tile[3];
      const index_type thr_id3  = threadIdx.y % m_rp.m_tile[3];

      temp0                   = m_rp.m_tile_end[4];
      temp1                   = m_rp.m_tile_end[5];
      const index_type numbl5 = (temp1 <= max_blocks ? temp1 : max_blocks);
      const index_type numbl4 =
          (temp0 * temp1 > max_blocks
               ? index_type(max_blocks / numbl5)
               : (temp0 <= max_blocks ? temp0 : max_blocks));

      const index_type tile_id4 = blockIdx.z / numbl5;
      const index_type tile_id5 = blockIdx.z % numbl5;
      const index_type thr_id4  = threadIdx.z / m_rp.m_tile[5];
      const index_type thr_id5  = threadIdx.z % m_rp.m_tile[5];

      for (index_type i = tile_id0; i < m_rp.m_tile_end[0]; i += numbl0) {
        const index_type offset_0 =
            i * m_rp.m_tile[0] + thr_id0 + (index_type)m_rp.m_lower[0];
        if (offset_0 < m_rp.m_upper[0] && thr_id0 < m_rp.m_tile[0]) {
          for (index_type j = tile_id1; j < m_rp.m_tile_end[1]; j += numbl1) {
            const index_type offset_1 =
                j * m_rp.m_tile[1] + thr_id1 + (index_type)m_rp.m_lower[1];
            if (offset_1 < m_rp.m_upper[1] && thr_id1 < m_rp.m_tile[1]) {
              for (index_type k = tile_id2; k < m_rp.m_tile_end[2];
                   k += numbl2) {
                const index_type offset_2 =
                    k * m_rp.m_tile[2] + thr_id2 + (index_type)m_rp.m_lower[2];
                if (offset_2 < m_rp.m_upper[2] && thr_id2 < m_rp.m_tile[2]) {
                  for (index_type l = tile_id3; l < m_rp.m_tile_end[3];
                       l += numbl3) {
                    const index_type offset_3 = l * m_rp.m_tile[3] + thr_id3 +
                                                (index_type)m_rp.m_lower[3];
                    if (offset_3 < m_rp.m_upper[3] &&
                        thr_id3 < m_rp.m_tile[3]) {
                      for (index_type m = tile_id4; m < m_rp.m_tile_end[4];
                           m += numbl4) {
                        const index_type offset_4 = m * m_rp.m_tile[4] +
                                                    thr_id4 +
                                                    (index_type)m_rp.m_lower[4];
                        if (offset_4 < m_rp.m_upper[4] &&
                            thr_id4 < m_rp.m_tile[4]) {
                          for (index_type n = tile_id5; n < m_rp.m_tile_end[5];
                               n += numbl5) {
                            const index_type offset_5 =
                                n * m_rp.m_tile[5] + thr_id5 +
                                (index_type)m_rp.m_lower[5];
                            if (offset_5 < m_rp.m_upper[5] &&
                                thr_id5 < m_rp.m_tile[5]) {
                              m_func(Tag(), offset_0, offset_1, offset_2,
                                     offset_3, offset_4, offset_5);
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
  const RP& m_rp;
  const Functor& m_func;
};

// ----------------------------------------------------------------------------------

template <typename RP, typename Functor, typename Tag>
struct DeviceIterateTile {
  using index_type       = typename RP::index_type;
  using array_index_type = typename RP::array_index_type;
  using point_type       = typename RP::point_type;

  struct VoidDummy {};
  typedef typename std::conditional<std::is_same<Tag, void>::value, VoidDummy,
                                    Tag>::type usable_tag;

  DeviceIterateTile(const RP& rp, const Functor& func)
      : m_rp{rp}, m_func{func} {}

 private:
  inline __device__ void apply() const {
    apply_impl<RP::rank, RP, Functor, Tag>(m_rp, m_func).exec_range();
  }  // end apply

 public:
  inline __device__ void operator()(void) const { this->apply(); }

  inline void execute() const {
    const array_index_type maxblocks =
        65535;  // not true for blockIdx.x for newer archs
    if (RP::rank == 2) {
      const dim3 block(m_rp.m_tile[0], m_rp.m_tile[1], 1);
      const dim3 grid(
          std::min((m_rp.m_upper[0] - m_rp.m_lower[0] + block.x - 1) / block.x,
                   maxblocks),
          std::min((m_rp.m_upper[1] - m_rp.m_lower[1] + block.y - 1) / block.y,
                   maxblocks),
          1);
      CudaLaunch<DeviceIterateTile>(*this, grid, block);
    } else if (RP::rank == 3) {
      const dim3 block(m_rp.m_tile[0], m_rp.m_tile[1], m_rp.m_tile[2]);
      const dim3 grid(
          std::min((m_rp.m_upper[0] - m_rp.m_lower[0] + block.x - 1) / block.x,
                   maxblocks),
          std::min((m_rp.m_upper[1] - m_rp.m_lower[1] + block.y - 1) / block.y,
                   maxblocks),
          std::min((m_rp.m_upper[2] - m_rp.m_lower[2] + block.z - 1) / block.z,
                   maxblocks));
      CudaLaunch<DeviceIterateTile>(*this, grid, block);
    } else if (RP::rank == 4) {
      // id0,id1 encoded within threadIdx.x; id2 to threadIdx.y; id3 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1], m_rp.m_tile[2],
                       m_rp.m_tile[3]);
      const dim3 grid(
          std::min(
              static_cast<index_type>(m_rp.m_tile_end[0] * m_rp.m_tile_end[1]),
              static_cast<index_type>(maxblocks)),
          std::min((m_rp.m_upper[2] - m_rp.m_lower[2] + block.y - 1) / block.y,
                   maxblocks),
          std::min((m_rp.m_upper[3] - m_rp.m_lower[3] + block.z - 1) / block.z,
                   maxblocks));
      CudaLaunch<DeviceIterateTile>(*this, grid, block);
    } else if (RP::rank == 5) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1],
                       m_rp.m_tile[2] * m_rp.m_tile[3], m_rp.m_tile[4]);
      const dim3 grid(
          std::min(
              static_cast<index_type>(m_rp.m_tile_end[0] * m_rp.m_tile_end[1]),
              static_cast<index_type>(maxblocks)),
          std::min(
              static_cast<index_type>(m_rp.m_tile_end[2] * m_rp.m_tile_end[3]),
              static_cast<index_type>(maxblocks)),
          std::min((m_rp.m_upper[4] - m_rp.m_lower[4] + block.z - 1) / block.z,
                   maxblocks));
      CudaLaunch<DeviceIterateTile>(*this, grid, block);
    } else if (RP::rank == 6) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4,id5 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1],
                       m_rp.m_tile[2] * m_rp.m_tile[3],
                       m_rp.m_tile[4] * m_rp.m_tile[5]);
      const dim3 grid(
          std::min(
              static_cast<index_type>(m_rp.m_tile_end[0] * m_rp.m_tile_end[1]),
              static_cast<index_type>(maxblocks)),
          std::min(
              static_cast<index_type>(m_rp.m_tile_end[2] * m_rp.m_tile_end[3]),
              static_cast<index_type>(maxblocks)),
          std::min(
              static_cast<index_type>(m_rp.m_tile_end[4] * m_rp.m_tile_end[5]),
              static_cast<index_type>(maxblocks)));
      CudaLaunch<DeviceIterateTile>(*this, grid, block);
    } else {
      printf("Kokkos::MDRange Error: Exceeded rank bounds with Cuda\n");
      Kokkos::abort("Aborting");
    }

  }  // end execute

 protected:
  const RP m_rp;
  const Functor m_func;
};

}  // namespace Impl
}  // namespace Kokkos

#endif
#endif
