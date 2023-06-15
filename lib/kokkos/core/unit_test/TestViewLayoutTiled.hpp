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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <cstdio>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_ViewLayoutTiled.hpp>

#include <type_traits>
#include <typeinfo>

namespace Test {

namespace {

template <typename ExecSpace>
struct TestViewLayoutTiled {
  using Scalar = double;

  static constexpr int T0 = 2;
  static constexpr int T1 = 4;
  static constexpr int T2 = 4;
  static constexpr int T3 = 2;
  static constexpr int T4 = 2;
  static constexpr int T5 = 2;
  static constexpr int T6 = 2;
  static constexpr int T7 = 2;

  // Rank 2
  using LayoutLL_2D_2x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Left,
                                        Kokkos::Iterate::Left, T0, T1>;
  using LayoutRL_2D_2x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Right,
                                        Kokkos::Iterate::Left, T0, T1>;
  using LayoutLR_2D_2x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Left,
                                        Kokkos::Iterate::Right, T0, T1>;
  using LayoutRR_2D_2x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Right,
                                        Kokkos::Iterate::Right, T0, T1>;

  // Rank 3
  using LayoutLL_3D_2x4x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Left,
                                        Kokkos::Iterate::Left, T0, T1, T2>;
  using LayoutRL_3D_2x4x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Right,
                                        Kokkos::Iterate::Left, T0, T1, T2>;
  using LayoutLR_3D_2x4x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Left,
                                        Kokkos::Iterate::Right, T0, T1, T2>;
  using LayoutRR_3D_2x4x4 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Right,
                                        Kokkos::Iterate::Right, T0, T1, T2>;

  // Rank 4
  using LayoutLL_4D_2x4x4x2 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Left,
                                        Kokkos::Iterate::Left, T0, T1, T2, T3>;
  using LayoutRL_4D_2x4x4x2 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Right,
                                        Kokkos::Iterate::Left, T0, T1, T2, T3>;
  using LayoutLR_4D_2x4x4x2 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Left,
                                        Kokkos::Iterate::Right, T0, T1, T2, T3>;
  using LayoutRR_4D_2x4x4x2 =
      Kokkos::Experimental::LayoutTiled<Kokkos::Iterate::Right,
                                        Kokkos::Iterate::Right, T0, T1, T2, T3>;

#if !defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  static void test_view_layout_tiled_2d(const int, const int) {
#else
  static void test_view_layout_tiled_2d(const int N0, const int N1) {
    const int FT = T0 * T1;

    const int NT0 = int(std::ceil(N0 / T0));
    const int NT1 = int(std::ceil(N1 / T1));
    // Test create_mirror_view, deep_copy
    // Create LL View
    {
      using ViewType =
          typename Kokkos::View<Scalar**, LayoutLL_2D_2x4, ExecSpace>;
      ViewType v("v", N0, N1);

      typename ViewType::HostMirror hv = Kokkos::create_mirror_view(v);

      // Initialize host-view
      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              hv(ti * T0 + i, tj * T1 + j) =
                  (ti + tj * NT0) * FT + (i + j * T0);
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(v, hv);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>,
          ExecSpace>
          mdrangepolicy({0, 0}, {NT0, NT1}, {T0, T1});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 2 LL", mdrangepolicy,
          KOKKOS_LAMBDA(const int ti, const int tj) {
            for (int j = 0; j < T1; ++j) {
              for (int i = 0; i < T0; ++i) {
                if ((ti * T0 + i < N0) && (tj * T1 + j < N1)) {
                  v(ti * T0 + i, tj * T1 + j) += 1;
                }
              }
            }
          });

      Kokkos::deep_copy(hv, v);

      long counter_subview = 0;
      long counter_inc     = 0;
      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          auto tile_subview = Kokkos::tile_subview(hv, ti, tj);
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              if (tile_subview(i, j) != hv(ti * T0 + i, tj * T1 + j)) {
                ++counter_subview;
              }
              if (tile_subview(i, j) !=
                  ((ti + tj * NT0) * FT + (i + j * T0) + 1)) {
                ++counter_inc;
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }

    // Create RL View
    {
      using ViewType =
          typename Kokkos::View<Scalar**, LayoutRL_2D_2x4, ExecSpace>;
      Kokkos::View<Scalar**, LayoutRL_2D_2x4, ExecSpace> v("v", N0, N1);

      typename ViewType::HostMirror hv = Kokkos::create_mirror_view(v);

      // Initialize host-view
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              hv(ti * T0 + i, tj * T1 + j) =
                  (ti * NT1 + tj) * FT + (i + j * T0);
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(v, hv);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Left>,
          ExecSpace>
          mdrangepolicy({0, 0}, {NT0, NT1}, {T0, T1});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 2 RL", mdrangepolicy,
          KOKKOS_LAMBDA(const int ti, const int tj) {
            for (int j = 0; j < T1; ++j) {
              for (int i = 0; i < T0; ++i) {
                if ((ti * T0 + i < N0) && (tj * T1 + j < N1)) {
                  v(ti * T0 + i, tj * T1 + j) += 1;
                }
              }
            }
          });

      Kokkos::deep_copy(hv, v);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          auto tile_subview = Kokkos::tile_subview(hv, ti, tj);
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              if (tile_subview(i, j) != hv(ti * T0 + i, tj * T1 + j)) {
                ++counter_subview;
              }
              if (tile_subview(i, j) !=
                  ((ti * NT1 + tj) * FT + (i + j * T0) + 1)) {
                ++counter_inc;
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create LR View
    {
      using ViewType =
          typename Kokkos::View<Scalar**, LayoutLR_2D_2x4, ExecSpace>;
      Kokkos::View<Scalar**, LayoutLR_2D_2x4, ExecSpace> v("v", N0, N1);

      typename ViewType::HostMirror hv = Kokkos::create_mirror_view(v);

      // Initialize host-view
      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              hv(ti * T0 + i, tj * T1 + j) =
                  (ti + tj * NT0) * FT + (i * T1 + j);
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(v, hv);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Right>,
          ExecSpace>
          mdrangepolicy({0, 0}, {NT0, NT1}, {T0, T1});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 2 LR", mdrangepolicy,
          KOKKOS_LAMBDA(const int ti, const int tj) {
            for (int j = 0; j < T1; ++j) {
              for (int i = 0; i < T0; ++i) {
                if ((ti * T0 + i < N0) && (tj * T1 + j < N1)) {
                  v(ti * T0 + i, tj * T1 + j) += 1;
                }
              }
            }
          });

      Kokkos::deep_copy(hv, v);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          auto tile_subview = Kokkos::tile_subview(hv, ti, tj);
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              if (tile_subview(i, j) != hv(ti * T0 + i, tj * T1 + j)) {
                ++counter_subview;
              }
              if (tile_subview(i, j) !=
                  ((ti + tj * NT0) * FT + (i * T1 + j) + 1)) {
                ++counter_inc;
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create RR View
    {
      using ViewType =
          typename Kokkos::View<Scalar**, LayoutRR_2D_2x4, ExecSpace>;
      Kokkos::View<Scalar**, LayoutRR_2D_2x4, ExecSpace> v("v", N0, N1);

      typename ViewType::HostMirror hv = Kokkos::create_mirror_view(v);

      // Initialize host-view
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              hv(ti * T0 + i, tj * T1 + j) =
                  (ti * NT1 + tj) * FT + (i * T1 + j);
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(v, hv);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Right>,
          ExecSpace>
          mdrangepolicy({0, 0}, {NT0, NT1}, {T0, T1});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 2 LR", mdrangepolicy,
          KOKKOS_LAMBDA(const int ti, const int tj) {
            for (int j = 0; j < T1; ++j) {
              for (int i = 0; i < T0; ++i) {
                if ((ti * T0 + i < N0) && (tj * T1 + j < N1)) {
                  v(ti * T0 + i, tj * T1 + j) += 1;
                }
              }
            }
          });

      Kokkos::deep_copy(hv, v);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          auto tile_subview = Kokkos::tile_subview(hv, ti, tj);
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              if (tile_subview(i, j) != hv(ti * T0 + i, tj * T1 + j)) {
                ++counter_subview;
              }
              if (tile_subview(i, j) !=
                  ((ti * NT1 + tj) * FT + (i * T1 + j) + 1)) {
                ++counter_inc;
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope
#endif
  }  // end test_view_layout_tiled_2d

#if !defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  static void test_view_layout_tiled_3d(const int, const int, const int) {
#else
  static void test_view_layout_tiled_3d(const int N0, const int N1,
                                        const int N2) {
    const int FT = T0 * T1 * T2;

    const int NT0 = int(std::ceil(N0 / T0));
    const int NT1 = int(std::ceil(N1 / T1));
    const int NT2 = int(std::ceil(N2 / T2));

    // Create LL View
    {
      using ViewType = Kokkos::View<Scalar***, LayoutLL_3D_2x4x4, ExecSpace>;
      Kokkos::View<Scalar***, LayoutLL_3D_2x4x4, ExecSpace> dv("dv", N0, N1,
                                                               N2);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti + tj * NT0 + tk * N0 * N1) * FT +
                      (i + j * T0 + k * T0 * T1);
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left>,
          ExecSpace>
          mdrangepolicy({0, 0, 0}, {N0, N1, N2}, {T0, T1, T2});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 3 LL", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k) {
            dv(i, j, k) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter_subview;
                  }
                  if (tile_subview(i, j, k) !=
                      ((ti + tj * NT0 + tk * N0 * N1) * FT +
                       (i + j * T0 + k * T0 * T1) + 1)) {
                    ++counter_inc;
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create RL View
    {
      using ViewType = Kokkos::View<Scalar***, LayoutRL_3D_2x4x4, ExecSpace>;
      Kokkos::View<Scalar***, LayoutRL_3D_2x4x4, ExecSpace> dv("dv", N0, N1,
                                                               N2);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                      (i + j * T0 + k * T0 * T1);
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Left>,
          ExecSpace>
          mdrangepolicy({0, 0, 0}, {N0, N1, N2}, {T0, T1, T2});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 3 RL", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k) {
            dv(i, j, k) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter_subview;
                  }
                  if (tile_subview(i, j, k) !=
                      ((ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                       (i + j * T0 + k * T0 * T1) + 1)) {
                    ++counter_inc;
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create LR View
    {
      using ViewType = Kokkos::View<Scalar***, LayoutLR_3D_2x4x4, ExecSpace>;
      Kokkos::View<Scalar***, LayoutLR_3D_2x4x4, ExecSpace> dv("dv", N0, N1,
                                                               N2);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti + tj * NT0 + tk * NT0 * NT1) * FT +
                      (i * T1 * T2 + j * T2 + k);
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Right>,
          ExecSpace>
          mdrangepolicy({0, 0, 0}, {N0, N1, N2}, {T0, T1, T2});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 3 LR", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k) {
            dv(i, j, k) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter_subview;
                  }
                  if (tile_subview(i, j, k) !=
                      ((ti + tj * NT0 + tk * NT0 * NT1) * FT +
                       (i * T1 * T2 + j * T2 + k) + 1)) {
                    ++counter_inc;
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create RR View
    {
      using ViewType = Kokkos::View<Scalar***, LayoutRR_3D_2x4x4, ExecSpace>;
      Kokkos::View<Scalar***, LayoutRR_3D_2x4x4, ExecSpace> dv("dv", N0, N1,
                                                               N2);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                      (i * T1 * T2 + j * T2 + k);
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>,
          ExecSpace>
          mdrangepolicy({0, 0, 0}, {N0, N1, N2}, {T0, T1, T2});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 3 RR", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k) {
            dv(i, j, k) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter_subview;
                  }
                  if (tile_subview(i, j, k) !=
                      ((ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                       (i * T1 * T2 + j * T2 + k) + 1)) {
                    ++counter_inc;
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope
#endif
  }  // end test_view_layout_tiled_3d

#if !defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  static void test_view_layout_tiled_4d(const int, const int, const int,
                                        const int){
#else
  static void test_view_layout_tiled_4d(const int N0, const int N1,
                                        const int N2, const int N3) {
    const int FT = T0 * T1 * T2 * T3;

    const int NT0 = int(std::ceil(N0 / T0));
    const int NT1 = int(std::ceil(N1 / T1));
    const int NT2 = int(std::ceil(N2 / T2));
    const int NT3 = int(std::ceil(N3 / T3));

    // Create LL View
    {
      using ViewType = Kokkos::View<Scalar****, LayoutLL_4D_2x4x4x2, ExecSpace>;
      Kokkos::View<Scalar****, LayoutLL_4D_2x4x4x2, ExecSpace> dv("dv", N0, N1,
                                                                  N2, N3);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti + tj * NT0 + tk * N0 * N1 + tl * N0 * N1 * N2) *
                              FT +
                          (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<4, Kokkos::Iterate::Left, Kokkos::Iterate::Left>,
          ExecSpace>
          mdrangepolicy({0, 0, 0, 0}, {N0, N1, N2, N3}, {T0, T1, T2, T3});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 4 LL", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
            dv(i, j, k, l) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter_subview;
                      }
                      if (tile_subview(i, j, k, l) !=
                          ((ti + tj * NT0 + tk * N0 * N1 + tl * N0 * N1 * N2) *
                               FT +
                           (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2) + 1)) {
                        ++counter_inc;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create RL View
    {
      using ViewType = Kokkos::View<Scalar****, LayoutRL_4D_2x4x4x2, ExecSpace>;
      Kokkos::View<Scalar****, LayoutRL_4D_2x4x4x2, ExecSpace> dv("dv", N0, N1,
                                                                  N2, N3);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti * NT1 * NT2 * N3 + tj * NT2 * N3 + tk * N3 + tl) *
                              FT +
                          (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<4, Kokkos::Iterate::Right, Kokkos::Iterate::Left>,
          ExecSpace>
          mdrangepolicy({0, 0, 0, 0}, {N0, N1, N2, N3}, {T0, T1, T2, T3});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 4 RL", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
            dv(i, j, k, l) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter_subview;
                      }
                      if (tile_subview(i, j, k, l) !=
                          ((ti * NT1 * NT2 * N3 + tj * NT2 * N3 + tk * N3 +
                            tl) *
                               FT +
                           (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2) + 1)) {
                        ++counter_inc;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create LR View
    {
      using ViewType = Kokkos::View<Scalar****, LayoutLR_4D_2x4x4x2, ExecSpace>;
      Kokkos::View<Scalar****, LayoutLR_4D_2x4x4x2, ExecSpace> dv("dv", N0, N1,
                                                                  N2, N3);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti + tj * NT0 + tk * NT0 * NT1 +
                           tl * NT0 * NT1 * NT2) *
                              FT +
                          (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<4, Kokkos::Iterate::Left, Kokkos::Iterate::Right>,
          ExecSpace>
          mdrangepolicy({0, 0, 0, 0}, {N0, N1, N2, N3}, {T0, T1, T2, T3});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 4 LR", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
            dv(i, j, k, l) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter_subview;
                      }
                      if (tile_subview(i, j, k, l) !=
                          ((ti + tj * NT0 + tk * NT0 * NT1 +
                            tl * NT0 * NT1 * NT2) *
                               FT +
                           (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l) + 1)) {
                        ++counter_inc;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope

    // Create RR View
    {
      using ViewType = Kokkos::View<Scalar****, LayoutRR_4D_2x4x4x2, ExecSpace>;
      Kokkos::View<Scalar****, LayoutRR_4D_2x4x4x2, ExecSpace> dv("dv", N0, N1,
                                                                  N2, N3);

      typename ViewType::HostMirror v = Kokkos::create_mirror_view(dv);

      // Initialize on host
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti * NT1 * NT2 * NT3 + tj * NT2 * NT3 + tk * NT3 +
                           tl) *
                              FT +
                          (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // copy to device
      Kokkos::deep_copy(dv, v);

      Kokkos::MDRangePolicy<
          Kokkos::Rank<4, Kokkos::Iterate::Right, Kokkos::Iterate::Right>,
          ExecSpace>
          mdrangepolicy({0, 0, 0, 0}, {N0, N1, N2, N3}, {T0, T1, T2, T3});

      // iterate by tile
      Kokkos::parallel_for(
          "ViewTile rank 4 RR", mdrangepolicy,
          KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
            dv(i, j, k, l) += 1;
          });

      Kokkos::deep_copy(v, dv);

      long counter_subview = 0;
      long counter_inc     = 0;

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter_subview;
                      }
                      if (tile_subview(i, j, k, l) !=
                          ((ti * NT1 * NT2 * NT3 + tj * NT2 * NT3 + tk * NT3 +
                            tl) *
                               FT +
                           (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l) + 1)) {
                        ++counter_inc;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      ASSERT_EQ(counter_subview, long(0));
      ASSERT_EQ(counter_inc, long(0));
    }  // end scope
#endif
  }  // end test_view_layout_tiled_4d

  static void test_view_layout_tiled_subtile_2d(const int N0, const int N1) {
    const int FT = T0 * T1;

    const int NT0 = int(std::ceil(N0 / T0));
    const int NT1 = int(std::ceil(N1 / T1));

    // Counter to check for errors at the end
    long counter[4] = {0};

    // Create LL View
    {
      Kokkos::View<Scalar**, LayoutLL_2D_2x4, Kokkos::HostSpace> v("v", N0, N1);
      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              v(ti * T0 + i, tj * T1 + j) = (ti + tj * NT0) * FT + (i + j * T0);
            }
          }
        }
      }

      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          auto tile_subview = Kokkos::tile_subview(v, ti, tj);
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              if (tile_subview(i, j) != v(ti * T0 + i, tj * T1 + j)) {
                ++counter[0];
              }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
              std::cout << "idx0,idx1 = " << ti * T0 + i << "," << tj * T1 + j
                        << std::endl;
              std::cout << "ti,tj,i,j: " << ti << "," << tj << "," << i << ","
                        << j << "  v = " << v(ti * T0 + i, tj * T1 + j)
                        << "  flat idx = "
                        << (ti + tj * NT0) * FT + (i + j * T0) << std::endl;
              std::cout << "subview_tile output = " << tile_subview(i, j)
                        << std::endl;
#endif
            }
          }
        }
      }
    }  // end scope

    // Create RL View
    {
      Kokkos::View<Scalar**, LayoutRL_2D_2x4, Kokkos::HostSpace> v("v", N0, N1);
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              v(ti * T0 + i, tj * T1 + j) = (ti * NT1 + tj) * FT + (i + j * T0);
            }
          }
        }
      }

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          auto tile_subview = Kokkos::tile_subview(v, ti, tj);
          for (int j = 0; j < T1; ++j) {
            for (int i = 0; i < T0; ++i) {
              if (tile_subview(i, j) != v(ti * T0 + i, tj * T1 + j)) {
                ++counter[1];
              }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
              std::cout << "idx0,idx1 = " << ti * T0 + i << "," << tj * T1 + j
                        << std::endl;
              std::cout << "ti,tj,i,j: " << ti << "," << tj << "," << i << ","
                        << j << "  v = " << v(ti * T0 + i, tj * T1 + j)
                        << "  flat idx = "
                        << (ti * NT1 + tj) * FT + (i + j * T0) << std::endl;
              std::cout << "subview_tile output = " << tile_subview(i, j)
                        << std::endl;
#endif
            }
          }
        }
      }
    }  // end scope

    // Create LR View
    {
      Kokkos::View<Scalar**, LayoutLR_2D_2x4, Kokkos::HostSpace> v("v", N0, N1);
      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              v(ti * T0 + i, tj * T1 + j) = (ti + tj * NT0) * FT + (i * T1 + j);
            }
          }
        }
      }

      for (int tj = 0; tj < NT1; ++tj) {
        for (int ti = 0; ti < NT0; ++ti) {
          auto tile_subview = Kokkos::tile_subview(v, ti, tj);
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              if (tile_subview(i, j) != v(ti * T0 + i, tj * T1 + j)) {
                ++counter[2];
              }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
              std::cout << "idx0,idx1 = " << ti * T0 + i << "," << tj * T1 + j
                        << std::endl;
              std::cout << "ti,tj,i,j: " << ti << "," << tj << "," << i << ","
                        << j << "  v = " << v(ti * T0 + i, tj * T1 + j)
                        << "  flat idx = "
                        << (ti + tj * NT0) * FT + (i * T1 + j) << std::endl;
              std::cout << "subview_tile output = " << tile_subview(i, j)
                        << std::endl;
#endif
            }
          }
        }
      }
    }  // end scope

    // Create RR View
    {
      Kokkos::View<Scalar**, LayoutRR_2D_2x4, Kokkos::HostSpace> v("v", N0, N1);
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              v(ti * T0 + i, tj * T1 + j) = (ti * NT1 + tj) * FT + (i * T1 + j);
            }
          }
        }
      }

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          auto tile_subview = Kokkos::tile_subview(v, ti, tj);
          for (int i = 0; i < T0; ++i) {
            for (int j = 0; j < T1; ++j) {
              if (tile_subview(i, j) != v(ti * T0 + i, tj * T1 + j)) {
                ++counter[3];
              }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
              std::cout << "idx0,idx1 = " << ti * T0 + i << "," << tj * T1 + j
                        << std::endl;
              std::cout << "ti,tj,i,j: " << ti << "," << tj << "," << i << ","
                        << j << "  v = " << v(ti * T0 + i, tj * T1 + j)
                        << "  flat idx = "
                        << (ti * NT1 + tj) * FT + (i * T1 + j) << std::endl;
              std::cout << "subview_tile output = " << tile_subview(i, j)
                        << std::endl;
              std::cout << "subview tile rank = " << Kokkos::rank(tile_subview)
                        << std::endl;
#endif
            }
          }
        }
      }
    }  // end scope

#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
    std::cout << "subview_tile vs view errors:\n"
              << " LL: " << counter[0] << " RL: " << counter[1]
              << " LR: " << counter[2] << " RR: " << counter[3] << std::endl;
#endif

    ASSERT_EQ(counter[0], long(0));
    ASSERT_EQ(counter[1], long(0));
    ASSERT_EQ(counter[2], long(0));
    ASSERT_EQ(counter[3], long(0));
  }  // end test_view_layout_tiled_subtile_2d

  static void test_view_layout_tiled_subtile_3d(const int N0, const int N1,
                                                const int N2) {
    const int FT = T0 * T1 * T2;

    const int NT0 = int(std::ceil(N0 / T0));
    const int NT1 = int(std::ceil(N1 / T1));
    const int NT2 = int(std::ceil(N2 / T2));

    // Counter to check for errors at the end
    long counter[4] = {0};
    // Create LL View
    {
      Kokkos::View<Scalar***, LayoutLL_3D_2x4x4, Kokkos::HostSpace> v("v", N0,
                                                                      N1, N2);
      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti + tj * NT0 + tk * N0 * N1) * FT +
                      (i + j * T0 + k * T0 * T1);
                }
              }
            }
          }
        }
      }

      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter[0];
                  }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                  std::cout << "idx0,idx1,idx2 = " << ti * T0 + i << ","
                            << tj * T1 + j << "," << tk * T2 + k << std::endl;
                  std::cout
                      << "ti,tj,tk,i,j,k: " << ti << "," << tj << "," << tk
                      << "," << i << "," << j << "," << k
                      << "  v = " << v(ti * T0 + i, tj * T1 + j, tk * T2 + k)
                      << "  flat idx = "
                      << (ti + tj * NT0 + tk * N0 * N1) * FT +
                             (i + j * T0 + k * T0 * T1)
                      << std::endl;
                  std::cout << "subview_tile output = " << tile_subview(i, j, k)
                            << std::endl;
                  std::cout
                      << "subview tile rank = " << Kokkos::rank(tile_subview)
                      << std::endl;
#endif
                }
              }
            }
          }
        }
      }
    }  // end scope

    // Create RL View
    {
      Kokkos::View<Scalar***, LayoutRL_3D_2x4x4, Kokkos::HostSpace> v("v", N0,
                                                                      N1, N2);
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                      (i + j * T0 + k * T0 * T1);
                }
              }
            }
          }
        }
      }

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int k = 0; k < T2; ++k) {
              for (int j = 0; j < T1; ++j) {
                for (int i = 0; i < T0; ++i) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter[1];
                  }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                  std::cout << "idx0,idx1,idx2 = " << ti * T0 + i << ","
                            << tj * T1 + j << "," << tk * T2 + k << std::endl;
                  std::cout
                      << "ti,tj,tk,i,j,k: " << ti << "," << tj << "," << tk
                      << "," << i << "," << j << "," << k
                      << "  v = " << v(ti * T0 + i, tj * T1 + j, tk * T2 + k)
                      << "  flat idx = "
                      << (ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                             (i + j * T0 + k * T0 * T1)
                      << std::endl;
                  std::cout << "subview_tile output = " << tile_subview(i, j, k)
                            << std::endl;
#endif
                }
              }
            }
          }
        }
      }
    }  // end scope

    // Create LR View
    {
      Kokkos::View<Scalar***, LayoutLR_3D_2x4x4, Kokkos::HostSpace> v("v", N0,
                                                                      N1, N2);
      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti + tj * NT0 + tk * NT0 * NT1) * FT +
                      (i * T1 * T2 + j * T2 + k);
                }
              }
            }
          }
        }
      }

      for (int tk = 0; tk < NT2; ++tk) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int ti = 0; ti < NT0; ++ti) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter[2];
                  }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                  std::cout << "idx0,idx1,idx2 = " << ti * T0 + i << ","
                            << tj * T1 + j << "," << tk * T2 + k << std::endl;
                  std::cout
                      << "ti,tj,tk,i,j,k: " << ti << "," << tj << "," << tk
                      << "," << i << "," << j << "," << k
                      << "  v = " << v(ti * T0 + i, tj * T1 + j, tk * T2 + k)
                      << "  flat idx = "
                      << (ti + tj * NT0 + tk * NT0 * NT1) * FT +
                             (i * T1 * T2 + j * T2 + k)
                      << std::endl;
                  std::cout << "subview_tile output = " << tile_subview(i, j, k)
                            << std::endl;
                  std::cout
                      << "subview tile rank = " << Kokkos::rank(tile_subview)
                      << std::endl;
#endif
                }
              }
            }
          }
        }
      }
    }  // end scope

    // Create RR View
    {
      Kokkos::View<Scalar***, LayoutRR_3D_2x4x4, Kokkos::HostSpace> v("v", N0,
                                                                      N1, N2);
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  v(ti * T0 + i, tj * T1 + j, tk * T2 + k) =
                      (ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                      (i * T1 * T2 + j * T2 + k);
                }
              }
            }
          }
        }
      }

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk);
            for (int i = 0; i < T0; ++i) {
              for (int j = 0; j < T1; ++j) {
                for (int k = 0; k < T2; ++k) {
                  if (tile_subview(i, j, k) !=
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k)) {
                    ++counter[3];
                  }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                  std::cout << "idx0,idx1,idx2 = " << ti * T0 + i << ","
                            << tj * T1 + j << "," << tk * T2 + k << std::endl;
                  std::cout
                      << "ti,tj,tk,i,j,k: " << ti << "," << tj << "," << tk
                      << "," << i << "," << j << "," << k
                      << "  v = " << v(ti * T0 + i, tj * T1 + j, tk * T2 + k)
                      << "  flat idx = "
                      << (ti * NT1 * NT2 + tj * NT2 + tk) * FT +
                             (i * T1 * T2 + j * T2 + k)
                      << std::endl;
                  std::cout << "subview_tile output = " << tile_subview(i, j, k)
                            << std::endl;
                  std::cout
                      << "subview tile rank = " << Kokkos::rank(tile_subview)
                      << std::endl;
#endif
                }
              }
            }
          }
        }
      }
    }  // end scope

#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
    std::cout << "subview_tile vs view errors:\n"
              << " LL: " << counter[0] << " RL: " << counter[1]
              << " LR: " << counter[2] << " RR: " << counter[3] << std::endl;
#endif

    ASSERT_EQ(counter[0], long(0));
    ASSERT_EQ(counter[1], long(0));
    ASSERT_EQ(counter[2], long(0));
    ASSERT_EQ(counter[3], long(0));

  }  // end test_view_layout_tiled_subtile_3d

  static void test_view_layout_tiled_subtile_4d(const int N0, const int N1,
                                                const int N2, const int N3) {
    const int FT = T0 * T1 * T2 * T3;

    const int NT0 = int(std::ceil(N0 / T0));
    const int NT1 = int(std::ceil(N1 / T1));
    const int NT2 = int(std::ceil(N2 / T2));
    const int NT3 = int(std::ceil(N3 / T3));

    // Counter to check for errors at the end
    long counter[4] = {0};
    // Create LL View
    {
      Kokkos::View<Scalar****, LayoutLL_4D_2x4x4x2, Kokkos::HostSpace> v(
          "v", N0, N1, N2, N3);
      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti + tj * NT0 + tk * N0 * N1 + tl * N0 * N1 * N2) *
                              FT +
                          (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2);
                    }
                  }
                }
              }
            }
          }
        }
      }

      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter[0];
                      }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                      std::cout << "idx0,idx1,idx2,idx3 = " << ti * T0 + i
                                << "," << tj * T1 + j << "," << tk * T2 + k
                                << "," << tl * T3 + l << std::endl;
                      std::cout
                          << "ti,tj,tk,tl: " << ti << "," << tj << "," << tk
                          << "," << tl << ","
                          << "  i,j,k,l: " << i << "," << j << "," << k << ","
                          << l << "  v = "
                          << v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                               tl * T3 + l)
                          << "  flat idx = "
                          << (ti + tj * NT0 + tk * N0 * N1 +
                              tl * N0 * N1 * N2) *
                                     FT +
                                 (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2)
                          << std::endl;
                      std::cout << "subview_tile output = "
                                << tile_subview(i, j, k, l) << std::endl;
                      std::cout << "subview tile rank = "
                                << Kokkos::rank(tile_subview) << std::endl;
#endif
                    }
                  }
                }
              }
            }
          }
        }
      }
    }  // end scope

    // Create RL View
    {
      Kokkos::View<Scalar****, LayoutRL_4D_2x4x4x2, Kokkos::HostSpace> v(
          "v", N0, N1, N2, N3);
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti * NT1 * NT2 * N3 + tj * NT2 * N3 + tk * N3 + tl) *
                              FT +
                          (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2);
                    }
                  }
                }
              }
            }
          }
        }
      }

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int l = 0; l < T3; ++l) {
                for (int k = 0; k < T2; ++k) {
                  for (int j = 0; j < T1; ++j) {
                    for (int i = 0; i < T0; ++i) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter[1];
                      }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                      std::cout << "idx0,idx1,idx2,idx3 = " << ti * T0 + i
                                << "," << tj * T1 + j << "," << tk * T2 + k
                                << "," << tl * T3 + l << std::endl;
                      std::cout
                          << "ti,tj,tk,tl: " << ti << "," << tj << "," << tk
                          << "," << tl << ","
                          << "  i,j,k,l: " << i << "," << j << "," << k << ","
                          << l << "  v = "
                          << v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                               tl * T3 + l)
                          << "  flat idx = "
                          << (ti * NT1 * NT2 * N3 + tj * NT2 * N3 + tk * N3 +
                              tl) * FT +
                                 (i + j * T0 + k * T0 * T1 + l * T0 * T1 * T2)
                          << std::endl;
                      std::cout << "subview_tile output = "
                                << tile_subview(i, j, k, l) << std::endl;
                      std::cout << "subview tile rank = "
                                << Kokkos::rank(tile_subview) << std::endl;
#endif
                    }
                  }
                }
              }
            }
          }
        }
      }
    }  // end scope

    // Create LR View
    {
      Kokkos::View<Scalar****, LayoutLR_4D_2x4x4x2, Kokkos::HostSpace> v(
          "v", N0, N1, N2, N3);
      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti + tj * NT0 + tk * NT0 * NT1 +
                           tl * NT0 * NT1 * NT2) *
                              FT +
                          (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l);
                    }
                  }
                }
              }
            }
          }
        }
      }

      for (int tl = 0; tl < NT3; ++tl) {
        for (int tk = 0; tk < NT2; ++tk) {
          for (int tj = 0; tj < NT1; ++tj) {
            for (int ti = 0; ti < NT0; ++ti) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter[2];
                      }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                      std::cout << "idx0,idx1,idx2,idx3 = " << ti * T0 + i
                                << "," << tj * T1 + j << "," << tk * T2 + k
                                << "," << tl * T3 + l << std::endl;
                      std::cout
                          << "ti,tj,tk,tl: " << ti << "," << tj << "," << tk
                          << "," << tl << ","
                          << "  i,j,k,l: " << i << "," << j << "," << k << ","
                          << l << "  v = "
                          << v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                               tl * T3 + l)
                          << "  flat idx = "
                          << (ti + tj * NT0 + tk * NT0 * NT1 +
                              tl * NT0 * NT1 * NT2) *
                                     FT +
                                 (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l)
                          << std::endl;
                      std::cout << "subview_tile output = "
                                << tile_subview(i, j, k, l) << std::endl;
                      std::cout << "subview tile rank = "
                                << Kokkos::rank(tile_subview) << std::endl;
#endif
                    }
                  }
                }
              }
            }
          }
        }
      }
    }  // end scope

    // Create RR View
    {
      Kokkos::View<Scalar****, LayoutRR_4D_2x4x4x2, Kokkos::HostSpace> v(
          "v", N0, N1, N2, N3);
      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      v(ti * T0 + i, tj * T1 + j, tk * T2 + k, tl * T3 + l) =
                          (ti * NT1 * NT2 * NT3 + tj * NT2 * NT3 + tk * NT3 +
                           tl) *
                              FT +
                          (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l);
                    }
                  }
                }
              }
            }
          }
        }
      }

      for (int ti = 0; ti < NT0; ++ti) {
        for (int tj = 0; tj < NT1; ++tj) {
          for (int tk = 0; tk < NT2; ++tk) {
            for (int tl = 0; tl < NT3; ++tl) {
              auto tile_subview = Kokkos::tile_subview(v, ti, tj, tk, tl);
              for (int i = 0; i < T0; ++i) {
                for (int j = 0; j < T1; ++j) {
                  for (int k = 0; k < T2; ++k) {
                    for (int l = 0; l < T3; ++l) {
                      if (tile_subview(i, j, k, l) !=
                          v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                            tl * T3 + l)) {
                        ++counter[3];
                      }
#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
                      std::cout << "idx0,idx1,idx2,idx3 = " << ti * T0 + i
                                << "," << tj * T1 + j << "," << tk * T2 + k
                                << "," << tl * T3 + l << std::endl;
                      std::cout
                          << "ti,tj,tk,tl: " << ti << "," << tj << "," << tk
                          << "," << tl << ","
                          << "  i,j,k,l: " << i << "," << j << "," << k << ","
                          << l << "  v = "
                          << v(ti * T0 + i, tj * T1 + j, tk * T2 + k,
                               tl * T3 + l)
                          << "  flat idx = "
                          << (ti * NT1 * NT2 * NT3 + tj * NT2 * NT3 + tk * NT3 +
                              tl) * FT +
                                 (i * T1 * T2 * T3 + j * T2 * T3 + k * T3 + l)
                          << std::endl;
                      std::cout << "subview_tile output = "
                                << tile_subview(i, j, k, l) << std::endl;
                      std::cout << "subview tile rank = "
                                << Kokkos::rank(tile_subview) << std::endl;
#endif
                    }
                  }
                }
              }
            }
          }
        }
      }
    }  // end scope

#ifdef KOKKOS_VERBOSE_LAYOUTTILED_OUTPUT
    std::cout << "subview_tile vs view errors:\n"
              << " LL: " << counter[0] << " RL: " << counter[1]
              << " LR: " << counter[2] << " RR: " << counter[3] << std::endl;
#endif

    ASSERT_EQ(counter[0], long(0));
    ASSERT_EQ(counter[1], long(0));
    ASSERT_EQ(counter[2], long(0));
    ASSERT_EQ(counter[3], long(0));

  }  // end test_view_layout_tiled_subtile_4d

};  // end TestViewLayoutTiled struct

}  // namespace

TEST(TEST_CATEGORY, view_layouttiled) {
  // These two examples are iterating by tile, then within a tile - not by
  // extents If N# is not a power of two, but want to iterate by tile then
  // within a tile, need to check that mapped index is within extent
  TestViewLayoutTiled<TEST_EXECSPACE>::test_view_layout_tiled_2d(4, 12);
  TestViewLayoutTiled<TEST_EXECSPACE>::test_view_layout_tiled_3d(4, 12, 16);
  TestViewLayoutTiled<TEST_EXECSPACE>::test_view_layout_tiled_4d(4, 12, 16, 12);
}
TEST(TEST_CATEGORY, view_layouttiled_subtile) {
  // These two examples are iterating by tile, then within a tile - not by
  // extents If N# is not a power of two, but want to iterate by tile then
  // within a tile, need to check that mapped index is within extent
  TestViewLayoutTiled<TEST_EXECSPACE>::test_view_layout_tiled_subtile_2d(4, 12);
  TestViewLayoutTiled<TEST_EXECSPACE>::test_view_layout_tiled_subtile_3d(4, 12,
                                                                         16);
  TestViewLayoutTiled<TEST_EXECSPACE>::test_view_layout_tiled_subtile_4d(
      4, 12, 16, 12);
}
}  // namespace Test

#undef KOKKOS_IMPL_PUBLIC_INCLUDE
