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

#include <Kokkos_Core.hpp>

namespace Test {

namespace TeamMDRange {

struct FillFlattenedIndex {
  explicit FillFlattenedIndex(int n0, int n1, int n2, int n3 = 1, int n4 = 1,
                              int n5 = 1, int n6 = 1, int n7 = 1)
      : initValue{n0, n1, n2, n3, n4, n5, n6, n7} {}

  KOKKOS_INLINE_FUNCTION
  int operator()(int n0, int n1, int n2, int n3 = 0, int n4 = 0, int n5 = 0,
                 int n6 = 0, int n7 = 0) const {
    return ((((((n7 * initValue[7] + n6) * initValue[6] + n5) * initValue[5] +
               n4) *
                  initValue[4] +
              n3) *
                 initValue[3] +
             n2) *
                initValue[2] +
            n1) *
               initValue[1] +
           n0;
  }

  int initValue[8];
};

struct TestTeamMDParallelFor {
  using DataType = int64_t;
  using DimsType = int[8];

  template <typename HostViewType, typename FillFunctor>
  static void check_result_3D(HostViewType h_view,
                              FillFunctor const& fillFunctor) {
    for (size_t i = 0; i < h_view.extent(0); ++i) {
      for (size_t j = 0; j < h_view.extent(1); ++j) {
        for (size_t k = 0; k < h_view.extent(2); ++k) {
          EXPECT_EQ(h_view(i, j, k), fillFunctor(i, j, k));
        }
      }
    }
  }

  template <typename HostViewType, typename FillFunctor>
  static void check_result_4D(HostViewType h_view, FillFunctor& fillFunctor) {
    for (size_t i = 0; i < h_view.extent(0); ++i) {
      for (size_t j = 0; j < h_view.extent(1); ++j) {
        for (size_t k = 0; k < h_view.extent(2); ++k) {
          for (size_t l = 0; l < h_view.extent(3); ++l) {
            EXPECT_EQ(h_view(i, j, k, l), fillFunctor(i, j, k, l));
          }
        }
      }
    }
  }

  template <typename HostViewType, typename FillFunctor>
  static void check_result_5D(HostViewType h_view, FillFunctor& fillFunctor) {
    for (size_t i = 0; i < h_view.extent(0); ++i) {
      for (size_t j = 0; j < h_view.extent(1); ++j) {
        for (size_t k = 0; k < h_view.extent(2); ++k) {
          for (size_t l = 0; l < h_view.extent(3); ++l) {
            for (size_t m = 0; m < h_view.extent(4); ++m) {
              EXPECT_EQ(h_view(i, j, k, l, m), fillFunctor(i, j, k, l, m));
            }
          }
        }
      }
    }
  }

  template <typename HostViewType, typename FillFunctor>
  static void check_result_6D(HostViewType h_view, FillFunctor& fillFunctor) {
    for (size_t i = 0; i < h_view.extent(0); ++i) {
      for (size_t j = 0; j < h_view.extent(1); ++j) {
        for (size_t k = 0; k < h_view.extent(2); ++k) {
          for (size_t l = 0; l < h_view.extent(3); ++l) {
            for (size_t m = 0; m < h_view.extent(4); ++m) {
              for (size_t n = 0; n < h_view.extent(5); ++n) {
                EXPECT_EQ(h_view(i, j, k, l, m, n),
                          fillFunctor(i, j, k, l, m, n));
              }
            }
          }
        }
      }
    }
  }

  template <typename HostViewType, typename FillFunctor>
  static void check_result_7D(HostViewType h_view, FillFunctor& fillFunctor) {
    for (size_t i = 0; i < h_view.extent(0); ++i) {
      for (size_t j = 0; j < h_view.extent(1); ++j) {
        for (size_t k = 0; k < h_view.extent(2); ++k) {
          for (size_t l = 0; l < h_view.extent(3); ++l) {
            for (size_t m = 0; m < h_view.extent(4); ++m) {
              for (size_t n = 0; n < h_view.extent(5); ++n) {
                for (size_t o = 0; o < h_view.extent(6); ++o) {
                  EXPECT_EQ(h_view(i, j, k, l, m, n, o),
                            fillFunctor(i, j, k, l, m, n, o));
                }
              }
            }
          }
        }
      }
    }
  }

  template <typename HostViewType, typename FillFunctor>
  static void check_result_8D(HostViewType h_view, FillFunctor& fillFunctor) {
    for (size_t i = 0; i < h_view.extent(0); ++i) {
      for (size_t j = 0; j < h_view.extent(1); ++j) {
        for (size_t k = 0; k < h_view.extent(2); ++k) {
          for (size_t l = 0; l < h_view.extent(3); ++l) {
            for (size_t m = 0; m < h_view.extent(4); ++m) {
              for (size_t n = 0; n < h_view.extent(5); ++n) {
                for (size_t o = 0; o < h_view.extent(6); ++o) {
                  for (size_t p = 0; p < h_view.extent(7); ++p) {
                    EXPECT_EQ(h_view(i, j, k, l, m, n, o, p),
                              fillFunctor(i, j, k, l, m, n, o, p));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
};

// If KOKKOS_ENABLE_CUDA_LAMBDA is off, extended lambdas used in parallel_for
// and parallel_reduce in these tests will not compile correctly
#if !defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)

template <typename ExecSpace>
struct TestTeamThreadMDRangeParallelFor : public TestTeamMDParallelFor {
  using TeamType = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_3D_TeamThreadMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType***, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];

    ViewType v("v", leagueSize, n0, n1);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<2, Direction>, TeamType>(
                  team, n0, n1);

          Kokkos::parallel_for(teamRange, [=](int i, int j) {
            v(leagueRank, i, j) += fillFlattenedIndex(leagueRank, i, j);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_3D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_4D_TeamThreadMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType****, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];

    ViewType v("v", leagueSize, n0, n1, n2);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3, Direction>, TeamType>(
                  team, n0, n1, n2);

          Kokkos::parallel_for(teamRange, [=](int i, int j, int k) {
            v(leagueRank, i, j, k) += fillFlattenedIndex(leagueRank, i, j, k);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_4D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_5D_TeamThreadMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType*****, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];

    ViewType v("v", leagueSize, n0, n1, n2, n3);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<4, Direction>, TeamType>(
                  team, n0, n1, n2, n3);

          Kokkos::parallel_for(teamRange, [=](int i, int j, int k, int l) {
            v(leagueRank, i, j, k, l) +=
                fillFlattenedIndex(leagueRank, i, j, k, l);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_5D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_6D_TeamThreadMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType******, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<5, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4);

          Kokkos::parallel_for(
              teamRange, [=](int i, int j, int k, int l, int m) {
                v(leagueRank, i, j, k, l, m) +=
                    fillFlattenedIndex(leagueRank, i, j, k, l, m);
              });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_6D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_7D_TeamThreadMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType*******, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<6, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5);

          Kokkos::parallel_for(
              teamRange, [=](int i, int j, int k, int l, int m, int n) {
                v(leagueRank, i, j, k, l, m, n) +=
                    fillFlattenedIndex(leagueRank, i, j, k, l, m, n);
              });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_7D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_8D_TeamThreadMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType********, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];
    int n6         = dims[7];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5, n6);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5,
                                          n6);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<7, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5, n6);

          Kokkos::parallel_for(
              teamRange, [=](int i, int j, int k, int l, int m, int n, int o) {
                v(leagueRank, i, j, k, l, m, n, o) +=
                    fillFlattenedIndex(leagueRank, i, j, k, l, m, n, o);
              });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_8D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_single_direction_test(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType***, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int n0 = dims[0];
    int n1 = dims[1];
    int n2 = dims[2];

    ViewType v("v", n0, n1, n2);
    FillFlattenedIndex fillFlattenedIndex(n0, n1, n2);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(1, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          auto teamRange =
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3, Direction>, TeamType>(
                  team, n0, n1, n2);

          Kokkos::parallel_for(teamRange, [=](int i, int j, int k) {
            v(i, j, k) += fillFlattenedIndex(i, j, k);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_3D(h_view, fillFlattenedIndex);
  }
};

template <typename ExecSpace>
struct TestThreadVectorMDRangeParallelFor : public TestTeamMDParallelFor {
  using TeamType = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_4D_ThreadVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType****, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];

    ViewType v("v", leagueSize, n0, n1, n2);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto teamRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<2, Direction>, TeamType>(
                  team, n1, n2);

          Kokkos::parallel_for(teamThreadRange, [=](int i) {
            Kokkos::parallel_for(teamRange, [=](int j, int k) {
              v(leagueRank, i, j, k) += fillFlattenedIndex(leagueRank, i, j, k);
            });
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_4D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_5D_ThreadVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType*****, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];

    ViewType v("v", leagueSize, n0, n1, n2, n3);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto teamRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<3, Direction>, TeamType>(
                  team, n1, n2, n3);

          Kokkos::parallel_for(teamThreadRange, [=](int i) {
            Kokkos::parallel_for(teamRange, [=](int j, int k, int l) {
              v(leagueRank, i, j, k, l) +=
                  fillFlattenedIndex(leagueRank, i, j, k, l);
            });
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_5D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_6D_ThreadVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType******, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto teamRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<4, Direction>, TeamType>(
                  team, n1, n2, n3, n4);

          Kokkos::parallel_for(teamThreadRange, [=](int i) {
            Kokkos::parallel_for(teamRange, [=](int j, int k, int l, int m) {
              v(leagueRank, i, j, k, l, m) +=
                  fillFlattenedIndex(leagueRank, i, j, k, l, m);
            });
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_6D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_7D_ThreadVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType*******, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto teamRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<5, Direction>, TeamType>(
                  team, n1, n2, n3, n4, n5);

          Kokkos::parallel_for(teamThreadRange, [=](int i) {
            Kokkos::parallel_for(
                teamRange, [=](int j, int k, int l, int m, int n) {
                  v(leagueRank, i, j, k, l, m, n) +=
                      fillFlattenedIndex(leagueRank, i, j, k, l, m, n);
                });
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_7D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_8D_ThreadVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType********, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];
    int n6         = dims[7];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5, n6);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5,
                                          n6);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto teamRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<6, Direction>, TeamType>(
                  team, n1, n2, n3, n4, n5, n6);

          Kokkos::parallel_for(teamThreadRange, [=](int i) {
            Kokkos::parallel_for(
                teamRange, [=](int j, int k, int l, int m, int n, int o) {
                  v(leagueRank, i, j, k, l, m, n, o) +=
                      fillFlattenedIndex(leagueRank, i, j, k, l, m, n, o);
                });
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_8D(h_view, fillFlattenedIndex);
  }
};

template <typename ExecSpace>
struct TestTeamVectorMDRangeParallelFor : public TestTeamMDParallelFor {
  using TeamType = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_3D_TeamVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType***, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];

    ViewType v("v", leagueSize, n0, n1);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<2, Direction>, TeamType>(
                  team, n0, n1);

          Kokkos::parallel_for(teamRange, [=](int i, int j) {
            v(leagueRank, i, j) += fillFlattenedIndex(leagueRank, i, j);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_3D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_4D_TeamVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType****, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];

    ViewType v("v", leagueSize, n0, n1, n2);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<3, Direction>, TeamType>(
                  team, n0, n1, n2);

          Kokkos::parallel_for(teamRange, [=](int i, int j, int k) {
            v(leagueRank, i, j, k) += fillFlattenedIndex(leagueRank, i, j, k);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_4D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_5D_TeamVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType*****, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];

    ViewType v("v", leagueSize, n0, n1, n2, n3);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<4, Direction>, TeamType>(
                  team, n0, n1, n2, n3);

          Kokkos::parallel_for(teamRange, [=](int i, int j, int k, int l) {
            v(leagueRank, i, j, k, l) +=
                fillFlattenedIndex(leagueRank, i, j, k, l);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_5D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_6D_TeamVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType******, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<5, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4);

          Kokkos::parallel_for(
              teamRange, [=](int i, int j, int k, int l, int m) {
                v(leagueRank, i, j, k, l, m) +=
                    fillFlattenedIndex(leagueRank, i, j, k, l, m);
              });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_6D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_7D_TeamVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType*******, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<6, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5);

          Kokkos::parallel_for(
              teamRange, [=](int i, int j, int k, int l, int m, int n) {
                v(leagueRank, i, j, k, l, m, n) +=
                    fillFlattenedIndex(leagueRank, i, j, k, l, m, n);
              });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_7D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_for_8D_TeamVectorMDRange(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType********, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];
    int n6         = dims[7];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5, n6);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5,
                                          n6);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          int leagueRank = team.league_rank();

          auto teamRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<7, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5, n6);

          Kokkos::parallel_for(
              teamRange, [=](int i, int j, int k, int l, int m, int n, int o) {
                v(leagueRank, i, j, k, l, m, n, o) +=
                    fillFlattenedIndex(leagueRank, i, j, k, l, m, n, o);
              });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_8D(h_view, fillFlattenedIndex);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_double_direction_test(DimsType const& dims) {
    using ViewType     = typename Kokkos::View<DataType****, ExecSpace>;
    using HostViewType = typename ViewType::HostMirror;

    int n0 = dims[0];
    int n1 = dims[1];
    int n2 = dims[2];
    int n3 = dims[3];

    ViewType v("v", n0, n1, n2, n3);
    FillFlattenedIndex fillFlattenedIndex(n0, n1, n2, n3);

    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(1, Kokkos::AUTO),
        KOKKOS_LAMBDA(const TeamType& team) {
          auto teamRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<4, Direction>, TeamType>(
                  team, n0, n1, n2, n3);

          Kokkos::parallel_for(teamRange, [=](int i, int j, int k, int l) {
            v(i, j, k, l) += fillFlattenedIndex(i, j, k, l);
          });
        });

    HostViewType h_view = Kokkos::create_mirror_view_and_copy(
        typename HostViewType::traits::memory_space(), v);

    check_result_4D(h_view, fillFlattenedIndex);
  }
};

struct TestTeamMDParallelReduce {
  using DataType = int64_t;
  using DimsType = int[8];

  template <typename F>
  constexpr static DataType get_expected_partial_sum(DimsType const& dims,
                                                     size_t maxRank, F const& f,
                                                     DimsType& indices,
                                                     size_t rank) {
    if (rank == maxRank) {
      return f(indices[0], indices[1], indices[2], indices[3], indices[4],
               indices[5], indices[6], indices[7]);
    }

    auto& index       = indices[rank];
    DataType accValue = 0;
    for (index = 0; index < dims[rank]; ++index) {
      accValue += get_expected_partial_sum(dims, maxRank, f, indices, rank + 1);
    }

    return accValue;
  }

  template <typename F>
  static DataType get_expected_sum(DimsType const& dims, size_t maxRank,
                                   F const& f) {
    DimsType indices = {};
    return get_expected_partial_sum(dims, maxRank, f, indices, 0);
  }
};

template <typename ExecSpace>
struct TestTeamThreadMDRangeParallelReduce : public TestTeamMDParallelReduce {
  using TeamType = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_3D_TeamThreadMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType***, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];

    ViewType v("v", leagueSize, n0, n1);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<3>>({0, 0, 0},
                                                          {leagueSize, n0, n1}),
        KOKKOS_LAMBDA(const int i, const int j, const int k) {
          v(i, j, k) = fillFlattenedIndex(i, j, k);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadMDRange<Kokkos::Rank<2, Direction>, TeamType>(
                  team, n0, n1),
              [=](const int& i, const int& j, DataType& threadSum) {
                threadSum += v(leagueRank, i, j);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 3, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_4D_TeamThreadMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType****, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];

    ViewType v("v", leagueSize, n0, n1, n2);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<4>>(
            {0, 0, 0, 0}, {leagueSize, n0, n1, n2}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
          v(i, j, k, l) = fillFlattenedIndex(i, j, k, l);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadMDRange<Kokkos::Rank<3, Direction>, TeamType>(
                  team, n0, n1, n2),
              [=](const int& i, const int& j, const int& k,
                  DataType& threadSum) { threadSum += v(leagueRank, i, j, k); },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 4, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_5D_TeamThreadMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType*****, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];

    ViewType v("v", leagueSize, n0, n1, n2, n3);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<5>>(
            {0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l,
                      const int m) {
          v(i, j, k, l, m) = fillFlattenedIndex(i, j, k, l, m);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadMDRange<Kokkos::Rank<4, Direction>, TeamType>(
                  team, n0, n1, n2, n3),
              [=](const int& i, const int& j, const int& k, const int& l,
                  DataType& threadSum) {
                threadSum += v(leagueRank, i, j, k, l);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 5, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_6D_TeamThreadMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType******, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
            {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l,
                      const int m, const int n) {
          v(i, j, k, l, m, n) = fillFlattenedIndex(i, j, k, l, m, n);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadMDRange<Kokkos::Rank<5, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4),
              [=](const int& i, const int& j, const int& k, const int& l,
                  const int& m, DataType& threadSum) {
                threadSum += v(leagueRank, i, j, k, l, m);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 6, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  // MDRangePolicy only allows up to rank of 6. Because of this, expectedSum
  // array had to be constructed from a nested parallel_for loop.
  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_7D_TeamThreadMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType*******, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5);
    auto mdRangePolicy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
        {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4});

    Kokkos::parallel_for(
        mdRangePolicy,
        KOKKOS_LAMBDA(const int leagueRank, const int i, const int j,
                      const int k, const int l, const int m) {
          for (int n = 0; n < n5; ++n) {
            v(leagueRank, i, j, k, l, m, n) =
                fillFlattenedIndex(leagueRank, i, j, k, l, m, n);
          }
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadMDRange<Kokkos::Rank<6, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5),
              [=](const int& i, const int& j, const int& k, const int& l,
                  const int& m, const int& n, DataType& threadSum) {
                threadSum += v(leagueRank, i, j, k, l, m, n);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 7, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_8D_TeamThreadMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType********, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];
    int n6         = dims[7];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5, n6);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5,
                                          n6);
    auto mdRangePolicy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
        {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4});

    Kokkos::parallel_for(
        mdRangePolicy,
        KOKKOS_LAMBDA(const int leagueRank, const int i, const int j,
                      const int k, const int l, const int m) {
          for (int n = 0; n < n5; ++n) {
            for (int o = 0; o < n6; ++o) {
              v(leagueRank, i, j, k, l, m, n, o) =
                  fillFlattenedIndex(leagueRank, i, j, k, l, m, n, o);
            }
          }
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadMDRange<Kokkos::Rank<7, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5, n6),
              [=](const int& i, const int& j, const int& k, const int& l,
                  const int& m, const int& n, const int& o,
                  DataType& threadSum) {
                threadSum += v(leagueRank, i, j, k, l, m, n, o);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 8, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }
};

template <typename ExecSpace>
struct TestThreadVectorMDRangeParallelReduce : public TestTeamMDParallelReduce {
  using TeamType = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_4D_ThreadVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType****, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];

    ViewType v("v", leagueSize, n0, n1, n2);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<4>>(
            {0, 0, 0, 0}, {leagueSize, n0, n1, n2}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
          v(i, j, k, l) = fillFlattenedIndex(i, j, k, l);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto threadVectorRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<2, Direction>, TeamType>(
                  team, n1, n2);

          Kokkos::parallel_for(teamThreadRange, [=, &teamSum](const int& i) {
            DataType threadSum = 0;
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](const int& j, const int& k, DataType& vectorSum) {
                  vectorSum += v(leagueRank, i, j, k);
                },
                threadSum);

            teamSum += threadSum;
          });

          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 4, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_5D_ThreadVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType*****, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];

    ViewType v("v", leagueSize, n0, n1, n2, n3);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<5>>(
            {0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l,
                      const int m) {
          v(i, j, k, l, m) = fillFlattenedIndex(i, j, k, l, m);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto threadVectorRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<3, Direction>, TeamType>(
                  team, n1, n2, n3);

          Kokkos::parallel_for(teamThreadRange, [=, &teamSum](const int& i) {
            DataType threadSum = 0;
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](const int& j, const int& k, const int& l,
                    DataType& vectorSum) {
                  vectorSum += v(leagueRank, i, j, k, l);
                },
                threadSum);

            teamSum += threadSum;
          });

          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 5, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_6D_ThreadVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType******, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
            {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l,
                      const int m, const int n) {
          v(i, j, k, l, m, n) = fillFlattenedIndex(i, j, k, l, m, n);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto threadVectorRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<4, Direction>, TeamType>(
                  team, n1, n2, n3, n4);

          Kokkos::parallel_for(teamThreadRange, [=, &teamSum](const int& i) {
            DataType threadSum = 0;
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](const int& j, const int& k, const int& l, const int& m,
                    DataType& vectorSum) {
                  vectorSum += v(leagueRank, i, j, k, l, m);
                },
                threadSum);

            teamSum += threadSum;
          });

          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 6, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_7D_ThreadVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType*******, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5);
    auto mdRangePolicy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
        {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4});

    Kokkos::parallel_for(
        mdRangePolicy,
        KOKKOS_LAMBDA(const int leagueRank, const int i, const int j,
                      const int k, const int l, const int m) {
          for (int n = 0; n < n5; ++n) {
            v(leagueRank, i, j, k, l, m, n) =
                fillFlattenedIndex(leagueRank, i, j, k, l, m, n);
          }
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto threadVectorRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<5, Direction>, TeamType>(
                  team, n1, n2, n3, n4, n5);

          Kokkos::parallel_for(teamThreadRange, [=, &teamSum](const int& i) {
            DataType threadSum = 0;
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](const int& j, const int& k, const int& l, const int& m,
                    const int& n, DataType& vectorSum) {
                  vectorSum += v(leagueRank, i, j, k, l, m, n);
                },
                threadSum);

            teamSum += threadSum;
          });

          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 7, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_8D_ThreadVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType********, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];
    int n6         = dims[7];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5, n6);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5,
                                          n6);
    auto mdRangePolicy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
        {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4});

    Kokkos::parallel_for(
        mdRangePolicy,
        KOKKOS_LAMBDA(const int leagueRank, const int i, const int j,
                      const int k, const int l, const int m) {
          for (int n = 0; n < n5; ++n) {
            for (int o = 0; o < n6; ++o) {
              v(leagueRank, i, j, k, l, m, n, o) =
                  fillFlattenedIndex(leagueRank, i, j, k, l, m, n, o);
            }
          }
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamThreadRange = Kokkos::TeamThreadRange(team, n0);
          auto threadVectorRange =
              Kokkos::ThreadVectorMDRange<Kokkos::Rank<6, Direction>, TeamType>(
                  team, n1, n2, n3, n4, n5, n6);

          Kokkos::parallel_for(teamThreadRange, [=, &teamSum](const int& i) {
            DataType threadSum = 0;
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](const int& j, const int& k, const int& l, const int& m,
                    const int& n, const int& o, DataType& vectorSum) {
                  vectorSum += v(leagueRank, i, j, k, l, m, n, o);
                },
                threadSum);

            teamSum += threadSum;
          });

          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 8, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }
};

template <typename ExecSpace>
struct TestTeamVectorMDRangeParallelReduce : public TestTeamMDParallelReduce {
  using TeamType = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_4D_TeamVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType****, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];

    ViewType v("v", leagueSize, n0, n1, n2);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<4>>(
            {0, 0, 0, 0}, {leagueSize, n0, n1, n2}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l) {
          v(i, j, k, l) = fillFlattenedIndex(i, j, k, l);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamVectorRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<3, Direction>, TeamType>(
                  team, n0, n1, n2);

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](const int& i, const int& j, const int& k,
                  DataType& vectorSum) { vectorSum += v(leagueRank, i, j, k); },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 4, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_5D_TeamVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType*****, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];

    ViewType v("v", leagueSize, n0, n1, n2, n3);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<5>>(
            {0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l,
                      const int m) {
          v(i, j, k, l, m) = fillFlattenedIndex(i, j, k, l, m);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamVectorRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<4, Direction>, TeamType>(
                  team, n0, n1, n2, n3);

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](const int& i, const int& j, const int& k, const int& l,
                  DataType& vectorSum) {
                vectorSum += v(leagueRank, i, j, k, l);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 5, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_6D_TeamVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType******, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4);

    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
            {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4}),
        KOKKOS_LAMBDA(const int i, const int j, const int k, const int l,
                      const int m, const int n) {
          v(i, j, k, l, m, n) = fillFlattenedIndex(i, j, k, l, m, n);
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamVectorRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<5, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4);

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](const int& i, const int& j, const int& k, const int& l,
                  const int& m, DataType& vectorSum) {
                vectorSum += v(leagueRank, i, j, k, l, m);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 6, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_7D_TeamVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType*******, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5);
    auto mdRangePolicy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
        {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4});

    Kokkos::parallel_for(
        mdRangePolicy,
        KOKKOS_LAMBDA(const int leagueRank, const int i, const int j,
                      const int k, const int l, const int m) {
          for (int n = 0; n < n5; ++n) {
            v(leagueRank, i, j, k, l, m, n) =
                fillFlattenedIndex(leagueRank, i, j, k, l, m, n);
          }
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamVectorRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<6, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5);

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](const int& i, const int& j, const int& k, const int& l,
                  const int& m, const int& n, DataType& vectorSum) {
                vectorSum += v(leagueRank, i, j, k, l, m, n);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 7, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }

  template <Kokkos::Iterate Direction = Kokkos::Iterate::Default>
  static void test_parallel_reduce_for_8D_TeamVectorMDRange(
      DimsType const& dims) {
    using ViewType = typename Kokkos::View<DataType********, ExecSpace>;

    int leagueSize = dims[0];
    int n0         = dims[1];
    int n1         = dims[2];
    int n2         = dims[3];
    int n3         = dims[4];
    int n4         = dims[5];
    int n5         = dims[6];
    int n6         = dims[7];

    ViewType v("v", leagueSize, n0, n1, n2, n3, n4, n5, n6);
    FillFlattenedIndex fillFlattenedIndex(leagueSize, n0, n1, n2, n3, n4, n5,
                                          n6);
    auto mdRangePolicy = Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<6>>(
        {0, 0, 0, 0, 0, 0}, {leagueSize, n0, n1, n2, n3, n4});

    Kokkos::parallel_for(
        mdRangePolicy,
        KOKKOS_LAMBDA(const int leagueRank, const int i, const int j,
                      const int k, const int l, const int m) {
          for (int n = 0; n < n5; ++n) {
            for (int o = 0; o < n6; ++o) {
              v(leagueRank, i, j, k, l, m, n, o) =
                  fillFlattenedIndex(leagueRank, i, j, k, l, m, n, o);
            }
          }
        });

    DataType finalSum = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamPolicy<ExecSpace>(leagueSize, Kokkos::AUTO),
        KOKKOS_LAMBDA(TeamType const& team, DataType& leagueSum) {
          auto leagueRank  = team.league_rank();
          DataType teamSum = 0;

          auto teamVectorRange =
              Kokkos::TeamVectorMDRange<Kokkos::Rank<7, Direction>, TeamType>(
                  team, n0, n1, n2, n3, n4, n5, n6);

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](const int& i, const int& j, const int& k, const int& l,
                  const int& m, const int& n, const int& o,
                  DataType& vectorSum) {
                vectorSum += v(leagueRank, i, j, k, l, m, n, o);
              },
              teamSum);
          leagueSum += teamSum;
        },
        finalSum);

    DataType expectedSum = get_expected_sum(dims, 8, fillFlattenedIndex);

    EXPECT_EQ(finalSum, expectedSum);
  }
};
/*--------------------------------------------------------------------------*/

constexpr auto Left  = Kokkos::Iterate::Left;
constexpr auto Right = Kokkos::Iterate::Right;

// Using prime numbers makes debugging easier
// small dimensions were needed for larger dimensions to reduce test run time
int dims[]      = {3, 5, 7, 11, 13, 17, 19, 23};
int smallDims[] = {2, 3, 2, 3, 5, 2, 3, 5};

TEST(TEST_CATEGORY, TeamThreadMDRangeParallelFor) {
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_3D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_3D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_4D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_4D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_5D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_5D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_6D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_6D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_7D_TeamThreadMDRange<Left>(smallDims);
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_7D_TeamThreadMDRange<Right>(smallDims);

  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_8D_TeamThreadMDRange<Left>(smallDims);
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_8D_TeamThreadMDRange<Right>(smallDims);

  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_single_direction_test<Left>(dims);
  TestTeamThreadMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_single_direction_test<Right>(dims);
}

TEST(TEST_CATEGORY, ThreadVectorMDRangeParallelFor) {
  TestThreadVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_4D_ThreadVectorMDRange<Left>(dims);
  TestThreadVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_4D_ThreadVectorMDRange<Right>(dims);

  TestThreadVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_5D_ThreadVectorMDRange<Left>(dims);
  TestThreadVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_5D_ThreadVectorMDRange<Right>(dims);

  TestThreadVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_6D_ThreadVectorMDRange<Left>(dims);
  TestThreadVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_6D_ThreadVectorMDRange<Right>(dims);

  TestThreadVectorMDRangeParallelFor<TEST_EXECSPACE>::
      test_parallel_for_7D_ThreadVectorMDRange<Left>(smallDims);
  TestThreadVectorMDRangeParallelFor<TEST_EXECSPACE>::
      test_parallel_for_7D_ThreadVectorMDRange<Right>(smallDims);

  TestThreadVectorMDRangeParallelFor<TEST_EXECSPACE>::
      test_parallel_for_8D_ThreadVectorMDRange<Left>(smallDims);
  TestThreadVectorMDRangeParallelFor<TEST_EXECSPACE>::
      test_parallel_for_8D_ThreadVectorMDRange<Right>(smallDims);
}

TEST(TEST_CATEGORY, TeamVectorMDRangeParallelFor) {
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_3D_TeamVectorMDRange<Left>(dims);
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_3D_TeamVectorMDRange<Right>(dims);

  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_4D_TeamVectorMDRange<Left>(dims);
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_4D_TeamVectorMDRange<Right>(dims);

  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_5D_TeamVectorMDRange<Left>(dims);
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_5D_TeamVectorMDRange<Right>(dims);

  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_6D_TeamVectorMDRange<Left>(dims);
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_6D_TeamVectorMDRange<Right>(dims);

  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_7D_TeamVectorMDRange<Left>(smallDims);
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_7D_TeamVectorMDRange<Right>(smallDims);

  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_8D_TeamVectorMDRange<Left>(smallDims);
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_for_8D_TeamVectorMDRange<Right>(smallDims);

  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_double_direction_test<Left>(dims);
  TestTeamVectorMDRangeParallelFor<
      TEST_EXECSPACE>::test_parallel_double_direction_test<Right>(dims);
}

TEST(TEST_CATEGORY, TeamThreadMDRangeParallelReduce) {
  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_3D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_3D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_4D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_4D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_5D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_5D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_6D_TeamThreadMDRange<Left>(dims);
  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_6D_TeamThreadMDRange<Right>(dims);

  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_7D_TeamThreadMDRange<Left>(smallDims);
  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_7D_TeamThreadMDRange<Right>(smallDims);

  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_8D_TeamThreadMDRange<Left>(smallDims);
  TestTeamThreadMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_8D_TeamThreadMDRange<Right>(smallDims);
}

TEST(TEST_CATEGORY, ThreadVectorMDRangeParallelReduce) {
// FIXME_SYCL sycl::group_barrier doesn't work correctly for non-Intel GPUs
#if defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ARCH_INTEL_GPU)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::SYCL>)
    GTEST_SKIP() << "skipping because of bug in group_barrier implementation";
#endif

// FIXME_OPENMPTARGET_CRAY: The unit tests fails correctness.
#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(KOKKOS_COMPILER_CRAYCLANG)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "Cray compiler fails correctness at runtime with the "
                    "OpenMPTarget backend.";
#endif

  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_4D_ThreadVectorMDRange<Left>(dims);
  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_4D_ThreadVectorMDRange<Right>(dims);

  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_5D_ThreadVectorMDRange<Left>(dims);
  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_5D_ThreadVectorMDRange<Right>(dims);

  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_6D_ThreadVectorMDRange<Left>(dims);
  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_6D_ThreadVectorMDRange<Right>(dims);

  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_7D_ThreadVectorMDRange<Left>(smallDims);
  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_7D_ThreadVectorMDRange<Right>(smallDims);

  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_8D_ThreadVectorMDRange<Left>(smallDims);
  TestThreadVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_8D_ThreadVectorMDRange<Right>(smallDims);
}

TEST(TEST_CATEGORY, TeamVectorMDRangeParallelReduce) {
// FIXME_SYCL sycl::group_barrier doesn't work correctly for non-Intel GPUs
#if defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ARCH_INTEL_GPU)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::SYCL>)
    GTEST_SKIP() << "skipping because of bug in group_barrier implementation";
#endif

// FIXME_OPENMPTARGET_CRAY: The unit tests fails correctness.
#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(KOKKOS_COMPILER_CRAYCLANG)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "Cray compiler fails correctness at runtime with the "
                    "OpenMPTarget backend.";
#endif

  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_4D_TeamVectorMDRange<Left>(dims);
  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_4D_TeamVectorMDRange<Right>(dims);

  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_5D_TeamVectorMDRange<Left>(dims);
  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_5D_TeamVectorMDRange<Right>(dims);

  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_6D_TeamVectorMDRange<Left>(dims);
  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_6D_TeamVectorMDRange<Right>(dims);

  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_7D_TeamVectorMDRange<Left>(smallDims);
  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_7D_TeamVectorMDRange<Right>(smallDims);

  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_8D_TeamVectorMDRange<Left>(smallDims);
  TestTeamVectorMDRangeParallelReduce<TEST_EXECSPACE>::
      test_parallel_reduce_for_8D_TeamVectorMDRange<Right>(smallDims);
}

#endif

}  // namespace TeamMDRange
}  // namespace Test
