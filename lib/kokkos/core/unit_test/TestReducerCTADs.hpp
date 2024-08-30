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

namespace {

struct TestReducerCTADS {
  using execspace   = TEST_EXECSPACE;
  using scalar_type = double;
  using index_type  = int;
  using memspace    = execspace::memory_space;

  struct CustomComparator {
    bool operator()(scalar_type, scalar_type) const;
  };
  static CustomComparator comparator;

  struct TestSum {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::Sum<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Sum(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Sum(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Sum(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Sum(unmanaged))>);
  };

  struct TestProd {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::Prod<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Prod(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Prod(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Prod(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Prod(unmanaged))>);
  };

  struct TestMin {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::Min<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Min(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Min(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Min(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Min(unmanaged))>);
  };

  struct TestMax {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::Max<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Max(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::Max(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Max(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::Max(unmanaged))>);
  };

  struct TestLAnd {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::LAnd<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::LAnd(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::LAnd(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::LAnd(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::LAnd(unmanaged))>);
  };

  struct TestLOr {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::LOr<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::LOr(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::LOr(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::LOr(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::LOr(unmanaged))>);
  };

  struct TestBAnd {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::BAnd<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::BAnd(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::BAnd(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::BAnd(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::BAnd(unmanaged))>);
  };

  struct TestBOr {
    static Kokkos::View<scalar_type, memspace> view;
    static Kokkos::View<scalar_type, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::BOr<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::BOr(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::BOr(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::BOr(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::BOr(unmanaged))>);
  };

  struct TestMinLoc {
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace>
        view;
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MinLoc<scalar_type, index_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::MinLoc(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::MinLoc(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinLoc(unmanaged))>);
  };

  struct TestMaxLoc {
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace>
        view;
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MaxLoc<scalar_type, index_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::MaxLoc(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::MaxLoc(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MaxLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MaxLoc(unmanaged))>);
  };

  struct TestMinMax {
    static Kokkos::View<Kokkos::MinMaxScalar<scalar_type>, memspace> view;
    static Kokkos::View<Kokkos::MinMaxScalar<scalar_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MinMax<scalar_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::MinMax(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::MinMax(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinMax(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinMax(unmanaged))>);
  };

  struct TestMinMaxLoc {
    static Kokkos::View<Kokkos::MinMaxLocScalar<scalar_type, index_type>,
                        memspace>
        view;
    static Kokkos::View<Kokkos::MinMaxLocScalar<scalar_type, index_type>,
                        memspace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MinMaxLoc<scalar_type, index_type, memspace> rt;

    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinMaxLoc(view))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinMaxLoc(rt))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MinMaxLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinMaxLoc(unmanaged))>);
  };

  struct TestMaxFirstLoc {
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace>
        view;
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MaxFirstLoc<scalar_type, index_type, memspace> rt;

    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MaxFirstLoc(view))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MaxFirstLoc(rt))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MaxFirstLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MaxFirstLoc(unmanaged))>);
  };

  struct TestMaxFirstLocCustomComparator {
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace>
        view;
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MaxFirstLocCustomComparator<scalar_type, index_type,
                                               CustomComparator, memspace>
        rt;

    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MaxFirstLocCustomComparator(
                                     view, comparator))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::MaxFirstLocCustomComparator(rt))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MaxFirstLocCustomComparator(
                                     std::move(rt)))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MaxFirstLocCustomComparator(
                                     unmanaged, comparator))>);
  };

  struct TestMinFirstLoc {
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace>
        view;
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MinFirstLoc<scalar_type, index_type, memspace> rt;

    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinFirstLoc(view))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinFirstLoc(rt))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MinFirstLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinFirstLoc(unmanaged))>);
  };

  struct TestMinFirstLocCustomComparator {
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace>
        view;
    static Kokkos::View<Kokkos::ValLocScalar<scalar_type, index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MinFirstLocCustomComparator<scalar_type, index_type,
                                               CustomComparator, memspace>
        rt;

    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MinFirstLocCustomComparator(
                                     view, comparator))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::MinFirstLocCustomComparator(rt))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MinFirstLocCustomComparator(
                                     std::move(rt)))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MinFirstLocCustomComparator(
                                     unmanaged, comparator))>);
  };

  struct TestMinMaxFirstLastLoc {
    static Kokkos::View<Kokkos::MinMaxLocScalar<scalar_type, index_type>,
                        memspace>
        view;
    static Kokkos::View<Kokkos::MinMaxLocScalar<scalar_type, index_type>,
                        memspace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MinMaxFirstLastLoc<scalar_type, index_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::MinMaxFirstLastLoc(view))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::MinMaxFirstLastLoc(rt))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::MinMaxFirstLastLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::MinMaxFirstLastLoc(unmanaged))>);
  };

  struct TestMinMaxFirstLastLocCustomComparator {
    static Kokkos::View<Kokkos::MinMaxLocScalar<scalar_type, index_type>,
                        memspace>
        view;
    static Kokkos::View<Kokkos::MinMaxLocScalar<scalar_type, index_type>,
                        memspace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::MinMaxFirstLastLocCustomComparator<
        scalar_type, index_type, CustomComparator, memspace>
        rt;

    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::MinMaxFirstLastLocCustomComparator(
                           view, comparator))>);
    static_assert(std::is_same_v<
                  decltype(rt),
                  decltype(Kokkos::MinMaxFirstLastLocCustomComparator(rt))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::MinMaxFirstLastLocCustomComparator(
                           std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::MinMaxFirstLastLocCustomComparator(
                           unmanaged, comparator))>);
  };

  struct TestFirstLoc {
    static Kokkos::View<Kokkos::FirstLocScalar<index_type>, memspace> view;
    static Kokkos::View<Kokkos::FirstLocScalar<index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::FirstLoc<index_type, memspace> rt;

    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::FirstLoc(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::FirstLoc(rt))>);
    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::FirstLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::FirstLoc(unmanaged))>);
  };

  struct TestLastLoc {
    static Kokkos::View<Kokkos::LastLocScalar<index_type>, memspace> view;
    static Kokkos::View<Kokkos::LastLocScalar<index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::LastLoc<index_type, memspace> rt;

    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::LastLoc(view))>);
    static_assert(std::is_same_v<decltype(rt), decltype(Kokkos::LastLoc(rt))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::LastLoc(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::LastLoc(unmanaged))>);
  };

  struct TestStdIsPartitioned {
    static Kokkos::View<Kokkos::StdIsPartScalar<index_type>, memspace> view;
    static Kokkos::View<Kokkos::StdIsPartScalar<index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::StdIsPartitioned<index_type, memspace> rt;

    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::StdIsPartitioned(view))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::StdIsPartitioned(rt))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::StdIsPartitioned(std::move(rt)))>);
    static_assert(std::is_same_v<
                  decltype(rt), decltype(Kokkos::StdIsPartitioned(unmanaged))>);
  };

  struct TestStdPartitionPoint {
    static Kokkos::View<Kokkos::StdPartPointScalar<index_type>, memspace> view;
    static Kokkos::View<Kokkos::StdPartPointScalar<index_type>, memspace,
                        Kokkos::MemoryTraits<Kokkos::Unmanaged>>
        unmanaged;
    static Kokkos::StdPartitionPoint<index_type, memspace> rt;

    static_assert(std::is_same_v<decltype(rt),
                                 decltype(Kokkos::StdPartitionPoint(view))>);
    static_assert(
        std::is_same_v<decltype(rt), decltype(Kokkos::StdPartitionPoint(rt))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::StdPartitionPoint(std::move(rt)))>);
    static_assert(
        std::is_same_v<decltype(rt),
                       decltype(Kokkos::StdPartitionPoint(unmanaged))>);
  };
};

}  // namespace
