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

#ifndef KOKKOS_OPENACC_PARALLEL_FOR_MDRANGE_HPP
#define KOKKOS_OPENACC_PARALLEL_FOR_MDRANGE_HPP

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>
#include <OpenACC/Kokkos_OpenACC_MDRangePolicy.hpp>
#include <Kokkos_Parallel.hpp>

namespace Kokkos::Experimental::Impl {

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<2> const& begin,
                                     OpenACCMDRangeEnd<2> const& end,
                                     int async_arg) {
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto i1 = m / dim0 + begin1;
    auto i0 = m % dim0 + begin0;
    functor(i0, i1);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(2) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i1 = begin1; i1 < end1; ++i1) {
    for (auto i0 = begin0; i0 < end0; ++i0) {
      functor(i0, i1);
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<2> const& begin,
                                     OpenACCMDRangeEnd<2> const& end,
                                     int async_arg) {
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto i0 = m / dim1 + begin0;
    auto i1 = m % dim1 + begin1;
    functor(i0, i1);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(2) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      functor(i0, i1);
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<2> const& begin,
                                     OpenACCMDRangeEnd<2> const& end,
                                     OpenACCMDRangeTile<2> const& tile,
                                     int async_arg) {
  auto tile0  = tile[0];
  auto tile1  = tile[1];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
// clang-format off
#pragma acc parallel loop gang vector tile(tile0,tile1) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i1 = begin1; i1 < end1; ++i1) {
    for (auto i0 = begin0; i0 < end0; ++i0) {
      functor(i0, i1);
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<2> const& begin,
                                     OpenACCMDRangeEnd<2> const& end,
                                     OpenACCMDRangeTile<2> const& tile,
                                     int async_arg) {
  auto tile1  = tile[1];
  auto tile0  = tile[0];
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
// clang-format off
#pragma acc parallel loop gang vector tile(tile1,tile0) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      functor(i0, i1);
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<3> const& begin,
                                     OpenACCMDRangeEnd<3> const& end,
                                     int async_arg) {
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim1 * dim0;
    auto i2   = m / tmp1 + begin2;
    auto tmp2 = m % tmp1;
    auto i1   = tmp2 / dim0 + begin1;
    auto i0   = tmp2 % dim0 + begin0;
    functor(i0, i1, i2);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(3) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i2 = begin2; i2 < end2; ++i2) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i0 = begin0; i0 < end0; ++i0) {
        functor(i0, i1, i2);
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<3> const& begin,
                                     OpenACCMDRangeEnd<3> const& end,
                                     int async_arg) {
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim2 * dim1;
    auto i0   = m / tmp1 + begin0;
    auto tmp2 = m % tmp1;
    auto i1   = tmp2 / dim2 + begin1;
    auto i2   = tmp2 % dim2 + begin2;
    functor(i0, i1, i2);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(3) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        functor(i0, i1, i2);
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<3> const& begin,
                                     OpenACCMDRangeEnd<3> const& end,
                                     OpenACCMDRangeTile<3> const& tile,
                                     int async_arg) {
  auto tile0  = tile[0];
  auto tile1  = tile[1];
  auto tile2  = tile[2];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
// clang-format off
#pragma acc parallel loop gang vector tile(tile0,tile1,tile2) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i2 = begin2; i2 < end2; ++i2) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i0 = begin0; i0 < end0; ++i0) {
        functor(i0, i1, i2);
      }
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<3> const& begin,
                                     OpenACCMDRangeEnd<3> const& end,
                                     OpenACCMDRangeTile<3> const& tile,
                                     int async_arg) {
  auto tile2  = tile[2];
  auto tile1  = tile[1];
  auto tile0  = tile[0];
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
// clang-format off
#pragma acc parallel loop gang vector tile(tile2,tile1,tile0) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        functor(i0, i1, i2);
      }
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<4> const& begin,
                                     OpenACCMDRangeEnd<4> const& end,
                                     int async_arg) {
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim3  = end3 - begin3;
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim3 * dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim2 * dim1 * dim0;
    auto i3   = m / tmp1 + begin3;
    auto tmp2 = m % tmp1;
    tmp1      = dim1 * dim0;
    auto i2   = tmp2 / tmp1 + begin2;
    tmp2      = tmp2 % tmp1;
    auto i1   = tmp2 / dim0 + begin1;
    auto i0   = tmp2 % dim0 + begin0;
    functor(i0, i1, i2, i3);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(4) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i3 = begin3; i3 < end3; ++i3) {
    for (auto i2 = begin2; i2 < end2; ++i2) {
      for (auto i1 = begin1; i1 < end1; ++i1) {
        for (auto i0 = begin0; i0 < end0; ++i0) {
          functor(i0, i1, i2, i3);
        }
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<4> const& begin,
                                     OpenACCMDRangeEnd<4> const& end,
                                     int async_arg) {
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin3 = begin[3];
  auto end3   = end[3];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim3  = end3 - begin3;
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim3 * dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim3 * dim2 * dim1;
    auto i0   = m / tmp1 + begin0;
    auto tmp2 = m % tmp1;
    tmp1      = dim3 * dim2;
    auto i1   = tmp2 / tmp1 + begin1;
    tmp2      = tmp2 % tmp1;
    auto i2   = tmp2 / dim3 + begin2;
    auto i3   = tmp2 % dim3 + begin3;
    functor(i0, i1, i2, i3);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(4) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i3 = begin3; i3 < end3; ++i3) {
          functor(i0, i1, i2, i3);
        }
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<4> const& begin,
                                     OpenACCMDRangeEnd<4> const& end,
                                     OpenACCMDRangeTile<4> const& tile,
                                     int async_arg) {
  auto tile0  = tile[0];
  auto tile1  = tile[1];
  auto tile2  = tile[2];
  auto tile3  = tile[3];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
// clang-format off
#pragma acc parallel loop gang vector tile(tile0,tile1,tile2,tile3) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i3 = begin3; i3 < end3; ++i3) {
    for (auto i2 = begin2; i2 < end2; ++i2) {
      for (auto i1 = begin1; i1 < end1; ++i1) {
        for (auto i0 = begin0; i0 < end0; ++i0) {
          functor(i0, i1, i2, i3);
        }
      }
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<4> const& begin,
                                     OpenACCMDRangeEnd<4> const& end,
                                     OpenACCMDRangeTile<4> const& tile,
                                     int async_arg) {
  auto tile3  = tile[3];
  auto tile2  = tile[2];
  auto tile1  = tile[1];
  auto tile0  = tile[0];
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin3 = begin[3];
  auto end3   = end[3];
// clang-format off
#pragma acc parallel loop gang vector tile(tile3,tile2,tile1,tile0) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i3 = begin3; i3 < end3; ++i3) {
          functor(i0, i1, i2, i3);
        }
      }
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<5> const& begin,
                                     OpenACCMDRangeEnd<5> const& end,
                                     int async_arg) {
  auto begin4 = begin[4];
  auto end4   = end[4];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim4  = end4 - begin4;
  auto dim3  = end3 - begin3;
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim4 * dim3 * dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim3 * dim2 * dim1 * dim0;
    auto i4   = m / tmp1 + begin4;
    auto tmp2 = m % tmp1;
    tmp1      = dim2 * dim1 * dim0;
    auto i3   = tmp2 / tmp1 + begin3;
    tmp2      = tmp2 % tmp1;
    tmp1      = dim1 * dim0;
    auto i2   = tmp2 / tmp1 + begin2;
    tmp2      = tmp2 % tmp1;
    auto i1   = tmp2 / dim0 + begin1;
    auto i0   = tmp2 % dim0 + begin0;
    functor(i0, i1, i2, i3, i4);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(5) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i4 = begin4; i4 < end4; ++i4) {
    for (auto i3 = begin3; i3 < end3; ++i3) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i1 = begin1; i1 < end1; ++i1) {
          for (auto i0 = begin0; i0 < end0; ++i0) {
            functor(i0, i1, i2, i3, i4);
          }
        }
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<5> const& begin,
                                     OpenACCMDRangeEnd<5> const& end,
                                     int async_arg) {
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin4 = begin[4];
  auto end4   = end[4];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim4  = end4 - begin4;
  auto dim3  = end3 - begin3;
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim4 * dim3 * dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim4 * dim3 * dim2 * dim1;
    auto i0   = m / tmp1 + begin0;
    auto tmp2 = m % tmp1;
    tmp1      = dim4 * dim3 * dim2;
    auto i1   = tmp2 / tmp1 + begin1;
    tmp2      = tmp2 % tmp1;
    tmp1      = dim4 * dim3;
    auto i2   = tmp2 / tmp1 + begin2;
    tmp2      = tmp2 % tmp1;
    auto i3   = tmp2 / dim4 + begin3;
    auto i4   = tmp2 % dim4 + begin4;
    functor(i0, i1, i2, i3, i4);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(5) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i3 = begin3; i3 < end3; ++i3) {
          for (auto i4 = begin4; i4 < end4; ++i4) {
            functor(i0, i1, i2, i3, i4);
          }
        }
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<5> const& begin,
                                     OpenACCMDRangeEnd<5> const& end,
                                     OpenACCMDRangeTile<5> const& tile,
                                     int async_arg) {
  auto tile0  = tile[0];
  auto tile1  = tile[1];
  auto tile2  = tile[2];
  auto tile3  = tile[3];
  auto tile4  = tile[4];
  auto begin4 = begin[4];
  auto end4   = end[4];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
// clang-format off
#pragma acc parallel loop gang vector tile(tile0,tile1,tile2,tile3,tile4) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i4 = begin4; i4 < end4; ++i4) {
    for (auto i3 = begin3; i3 < end3; ++i3) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i1 = begin1; i1 < end1; ++i1) {
          for (auto i0 = begin0; i0 < end0; ++i0) {
            functor(i0, i1, i2, i3, i4);
          }
        }
      }
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<5> const& begin,
                                     OpenACCMDRangeEnd<5> const& end,
                                     OpenACCMDRangeTile<5> const& tile,
                                     int async_arg) {
  auto tile4  = tile[4];
  auto tile3  = tile[3];
  auto tile2  = tile[2];
  auto tile1  = tile[1];
  auto tile0  = tile[0];
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin4 = begin[4];
  auto end4   = end[4];
// clang-format off
#pragma acc parallel loop gang vector tile(tile4,tile3,tile2,tile1,tile0) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i3 = begin3; i3 < end3; ++i3) {
          for (auto i4 = begin4; i4 < end4; ++i4) {
            functor(i0, i1, i2, i3, i4);
          }
        }
      }
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<6> const& begin,
                                     OpenACCMDRangeEnd<6> const& end,
                                     int async_arg) {
  auto begin5 = begin[5];
  auto end5   = end[5];
  auto begin4 = begin[4];
  auto end4   = end[4];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim5  = end5 - begin5;
  auto dim4  = end4 - begin4;
  auto dim3  = end3 - begin3;
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim5 * dim4 * dim3 * dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim4 * dim3 * dim2 * dim1 * dim0;
    auto i5   = m / tmp1 + begin5;
    auto tmp2 = m % tmp1;
    tmp1      = dim3 * dim2 * dim1 * dim0;
    auto i4   = tmp2 / tmp1 + begin4;
    tmp2      = tmp2 % tmp1;
    tmp1      = dim2 * dim1 * dim0;
    auto i3   = tmp2 / tmp1 + begin3;
    tmp2      = tmp2 % tmp1;
    tmp1      = dim1 * dim0;
    auto i2   = tmp2 / tmp1 + begin2;
    tmp2      = tmp2 % tmp1;
    auto i1   = tmp2 / dim0 + begin1;
    auto i0   = tmp2 % dim0 + begin0;
    functor(i0, i1, i2, i3, i4, i5);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(6) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i5 = begin5; i5 < end5; ++i5) {
    for (auto i4 = begin4; i4 < end4; ++i4) {
      for (auto i3 = begin3; i3 < end3; ++i3) {
        for (auto i2 = begin2; i2 < end2; ++i2) {
          for (auto i1 = begin1; i1 < end1; ++i1) {
            for (auto i0 = begin0; i0 < end0; ++i0) {
              functor(i0, i1, i2, i3, i4, i5);
            }
          }
        }
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<6> const& begin,
                                     OpenACCMDRangeEnd<6> const& end,
                                     int async_arg) {
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin4 = begin[4];
  auto end4   = end[4];
  auto begin5 = begin[5];
  auto end5   = end[5];
#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)
  auto dim5  = end5 - begin5;
  auto dim4  = end4 - begin4;
  auto dim3  = end3 - begin3;
  auto dim2  = end2 - begin2;
  auto dim1  = end1 - begin1;
  auto dim0  = end0 - begin0;
  auto nIter = dim5 * dim4 * dim3 * dim2 * dim1 * dim0;
// clang-format off
#pragma acc parallel loop gang vector copyin(functor) async(async_arg)
  // clang-format on
  for (decltype(nIter) m = 0; m < nIter; ++m) {
    auto tmp1 = dim5 * dim4 * dim3 * dim2 * dim1;
    auto i0   = m / tmp1 + begin0;
    auto tmp2 = m % tmp1;
    tmp1      = dim5 * dim4 * dim3 * dim2;
    auto i1   = tmp2 / tmp1 + begin1;
    tmp2      = tmp2 % tmp1;
    tmp1      = dim5 * dim4 * dim3;
    auto i2   = tmp2 / tmp1 + begin2;
    tmp2      = tmp2 % tmp1;
    tmp1      = dim5 * dim4;
    auto i3   = tmp2 / tmp1 + begin3;
    tmp2      = tmp2 % tmp1;
    auto i4   = tmp2 / dim5 + begin4;
    auto i5   = tmp2 % dim5 + begin5;
    functor(i0, i1, i2, i3, i4, i5);
  }
#else
// clang-format off
#pragma acc parallel loop gang vector collapse(6) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i3 = begin3; i3 < end3; ++i3) {
          for (auto i4 = begin4; i4 < end4; ++i4) {
            for (auto i5 = begin5; i5 < end5; ++i5) {
              functor(i0, i1, i2, i3, i4, i5);
            }
          }
        }
      }
    }
  }
#endif
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<6> const& begin,
                                     OpenACCMDRangeEnd<6> const& end,
                                     OpenACCMDRangeTile<6> const& tile,
                                     int async_arg) {
  auto tile0  = tile[0];
  auto tile1  = tile[1];
  auto tile2  = tile[2];
  auto tile3  = tile[3];
  auto tile4  = tile[4];
  auto tile5  = tile[5];
  auto begin5 = begin[5];
  auto end5   = end[5];
  auto begin4 = begin[4];
  auto end4   = end[4];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin0 = begin[0];
  auto end0   = end[0];
// clang-format off
#pragma acc parallel loop gang vector tile(tile0,tile1,tile2,tile3,tile4,tile5) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i5 = begin5; i5 < end5; ++i5) {
    for (auto i4 = begin4; i4 < end4; ++i4) {
      for (auto i3 = begin3; i3 < end3; ++i3) {
        for (auto i2 = begin2; i2 < end2; ++i2) {
          for (auto i1 = begin1; i1 < end1; ++i1) {
            for (auto i0 = begin0; i0 < end0; ++i0) {
              functor(i0, i1, i2, i3, i4, i5);
            }
          }
        }
      }
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<6> const& begin,
                                     OpenACCMDRangeEnd<6> const& end,
                                     OpenACCMDRangeTile<6> const& tile,
                                     int async_arg) {
  auto tile5  = tile[5];
  auto tile4  = tile[4];
  auto tile3  = tile[3];
  auto tile2  = tile[2];
  auto tile1  = tile[1];
  auto tile0  = tile[0];
  auto begin0 = begin[0];
  auto end0   = end[0];
  auto begin1 = begin[1];
  auto end1   = end[1];
  auto begin2 = begin[2];
  auto end2   = end[2];
  auto begin3 = begin[3];
  auto end3   = end[3];
  auto begin4 = begin[4];
  auto end4   = end[4];
  auto begin5 = begin[5];
  auto end5   = end[5];
// clang-format off
#pragma acc parallel loop gang vector tile(tile5,tile4,tile3,tile2,tile1,tile0) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      for (auto i2 = begin2; i2 < end2; ++i2) {
        for (auto i3 = begin3; i3 < end3; ++i3) {
          for (auto i4 = begin4; i4 < end4; ++i4) {
            for (auto i5 = begin5; i5 < end5; ++i5) {
              functor(i0, i1, i2, i3, i4, i5);
            }
          }
        }
      }
    }
  }
}

}  // namespace Kokkos::Experimental::Impl

template <class Functor, class... Traits>
class Kokkos::Impl::ParallelFor<Functor, Kokkos::MDRangePolicy<Traits...>,
                                Kokkos::Experimental::OpenACC> {
  using Policy = MDRangePolicy<Traits...>;
  Kokkos::Experimental::Impl::FunctorAdapter<
      Functor, Policy, Kokkos::Experimental::Impl::RoutineClause::seq>
      m_functor;
  Policy m_policy;

 public:
  ParallelFor(Functor const& functor, Policy const& policy)
      : m_functor(functor), m_policy(policy) {}

  void execute() const {
    static_assert(1 < Policy::rank && Policy::rank < 7);
    static_assert(Policy::inner_direction == Iterate::Left ||
                  Policy::inner_direction == Iterate::Right);
    constexpr int rank = Policy::rank;
    for (int i = 0; i < rank; ++i) {
      if (m_policy.m_lower[i] >= m_policy.m_upper[i]) {
        return;
      }
    }
    int const async_arg = m_policy.space().acc_async_queue();
#if 0  // FIXME_OPENACC: OpenACC requires tile size to be constant.
    for (int i = 0; i < rank; ++i) {
      if (m_policy.m_tile[i] < 1) {
        Kokkos::Experimental::Impl::OpenACCParallelForMDRangePolicy(
            Kokkos::Experimental::Impl::OpenACCCollapse(),
            std::integral_constant<Iterate, Policy::inner_direction>(),
            m_functor, m_policy.m_lower, m_policy.m_upper, async_arg);
        return;
      }
    }
    Kokkos::Experimental::Impl::OpenACCParallelForMDRangePolicy(
        Kokkos::Experimental::Impl::OpenACCTile(),
        std::integral_constant<Iterate, Policy::inner_direction>(), m_functor,
        m_policy.m_lower, m_policy.m_upper, m_policy.m_tile, async_arg);
#else
    Kokkos::Experimental::Impl::OpenACCParallelForMDRangePolicy(
        Kokkos::Experimental::Impl::OpenACCCollapse(),
        std::integral_constant<Iterate, Policy::inner_direction>(), m_functor,
        m_policy.m_lower, m_policy.m_upper, async_arg);
#endif
  }
};

#endif
