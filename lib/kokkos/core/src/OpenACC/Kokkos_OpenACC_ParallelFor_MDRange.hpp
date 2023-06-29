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
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
// clang-format off
#pragma acc parallel loop gang vector collapse(2) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i1 = begin1; i1 < end1; ++i1) {
    for (auto i0 = begin0; i0 < end0; ++i0) {
      functor(i0, i1);
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<2> const& begin,
                                     OpenACCMDRangeEnd<2> const& end,
                                     int async_arg) {
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
// clang-format off
#pragma acc parallel loop gang vector collapse(2) copyin(functor) async(async_arg)
  // clang-format on
  for (auto i0 = begin0; i0 < end0; ++i0) {
    for (auto i1 = begin1; i1 < end1; ++i1) {
      functor(i0, i1);
    }
  }
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<2> const& begin,
                                     OpenACCMDRangeEnd<2> const& end,
                                     OpenACCMDRangeTile<2> const& tile,
                                     int async_arg) {
  int tile0  = tile[0];
  int tile1  = tile[1];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
  int tile1  = tile[1];
  int tile0  = tile[0];
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
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
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<3> const& begin,
                                     OpenACCMDRangeEnd<3> const& end,
                                     int async_arg) {
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<3> const& begin,
                                     OpenACCMDRangeEnd<3> const& end,
                                     OpenACCMDRangeTile<3> const& tile,
                                     int async_arg) {
  int tile0  = tile[0];
  int tile1  = tile[1];
  int tile2  = tile[2];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
  int tile2  = tile[2];
  int tile1  = tile[1];
  int tile0  = tile[0];
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
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
  int begin3 = begin[3];
  int end3   = end[3];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<4> const& begin,
                                     OpenACCMDRangeEnd<4> const& end,
                                     int async_arg) {
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin3 = begin[3];
  int end3   = end[3];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<4> const& begin,
                                     OpenACCMDRangeEnd<4> const& end,
                                     OpenACCMDRangeTile<4> const& tile,
                                     int async_arg) {
  int tile0  = tile[0];
  int tile1  = tile[1];
  int tile2  = tile[2];
  int tile3  = tile[3];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
  int tile3  = tile[3];
  int tile2  = tile[2];
  int tile1  = tile[1];
  int tile0  = tile[0];
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin3 = begin[3];
  int end3   = end[3];
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
  int begin4 = begin[4];
  int end4   = end[4];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<5> const& begin,
                                     OpenACCMDRangeEnd<5> const& end,
                                     int async_arg) {
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin4 = begin[4];
  int end4   = end[4];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<5> const& begin,
                                     OpenACCMDRangeEnd<5> const& end,
                                     OpenACCMDRangeTile<5> const& tile,
                                     int async_arg) {
  int tile0  = tile[0];
  int tile1  = tile[1];
  int tile2  = tile[2];
  int tile3  = tile[3];
  int tile4  = tile[4];
  int begin4 = begin[4];
  int end4   = end[4];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
  int tile4  = tile[4];
  int tile3  = tile[3];
  int tile2  = tile[2];
  int tile1  = tile[1];
  int tile0  = tile[0];
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin4 = begin[4];
  int end4   = end[4];
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
  int begin5 = begin[5];
  int end5   = end[5];
  int begin4 = begin[4];
  int end4   = end[4];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCCollapse, OpenACCIterateRight,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<6> const& begin,
                                     OpenACCMDRangeEnd<6> const& end,
                                     int async_arg) {
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin4 = begin[4];
  int end4   = end[4];
  int begin5 = begin[5];
  int end5   = end[5];
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
}

template <class Functor>
void OpenACCParallelForMDRangePolicy(OpenACCTile, OpenACCIterateLeft,
                                     Functor const& functor,
                                     OpenACCMDRangeBegin<6> const& begin,
                                     OpenACCMDRangeEnd<6> const& end,
                                     OpenACCMDRangeTile<6> const& tile,
                                     int async_arg) {
  int tile0  = tile[0];
  int tile1  = tile[1];
  int tile2  = tile[2];
  int tile3  = tile[3];
  int tile4  = tile[4];
  int tile5  = tile[5];
  int begin5 = begin[5];
  int end5   = end[5];
  int begin4 = begin[4];
  int end4   = end[4];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin0 = begin[0];
  int end0   = end[0];
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
  int tile5  = tile[5];
  int tile4  = tile[4];
  int tile3  = tile[3];
  int tile2  = tile[2];
  int tile1  = tile[1];
  int tile0  = tile[0];
  int begin0 = begin[0];
  int end0   = end[0];
  int begin1 = begin[1];
  int end1   = end[1];
  int begin2 = begin[2];
  int end2   = end[2];
  int begin3 = begin[3];
  int end3   = end[3];
  int begin4 = begin[4];
  int end4   = end[4];
  int begin5 = begin[5];
  int end5   = end[5];
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
