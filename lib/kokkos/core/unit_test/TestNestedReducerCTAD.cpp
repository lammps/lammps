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

struct TestNestedReducerCTAD {
  using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;
  using ScalarType  = int;
  using IndexType   = int;
  using TeamPolicy  = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamHandle  = TeamPolicy::member_type;

  struct FakeComparator {
    template <class T>
    KOKKOS_FUNCTION bool operator()(T const&, T const&) const {
      return true;
    }
  };

  template <class ValueType>
  struct FakeFunctor {
    KOKKOS_FUNCTION void operator()(int, ValueType&) const {}
  };

  template <class ReducerTypeExpected, class ReducerTypeToCheck>
  KOKKOS_FUNCTION static void check_types([
      [maybe_unused]] ReducerTypeToCheck const& reducer) {
    static_assert(std::is_same_v<ReducerTypeExpected, ReducerTypeToCheck>);
  }

  KOKKOS_FUNCTION void operator()([
      [maybe_unused]] TeamHandle const& team_handle) const {
    {
      using ReducerTypeExpected = Kokkos::Sum<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::Sum reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::Prod<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::Prod reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::Min<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::Min reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::Max<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::Max reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::LAnd<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::LAnd reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::LOr<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::LOr reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::BAnd<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::BAnd reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::BOr<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::BOr reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MinLoc<ScalarType, IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::MinLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MaxLoc<ScalarType, IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::MaxLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::MinMax<ScalarType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::MinMax reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MinMaxLoc<ScalarType, IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::MinMaxLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MaxFirstLoc<ScalarType, IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::MaxFirstLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MaxFirstLocCustomComparator<ScalarType, IndexType,
                                              FakeComparator, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      FakeComparator comparator;
      Kokkos::MaxFirstLocCustomComparator reducer(view, comparator);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MinFirstLoc<ScalarType, IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::MinFirstLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MinFirstLocCustomComparator<ScalarType, IndexType,
                                              FakeComparator, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      FakeComparator comparator;
      Kokkos::MinFirstLocCustomComparator reducer(view, comparator);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::MinMaxFirstLastLoc<ScalarType, IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::MinMaxFirstLastLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::MinMaxFirstLastLocCustomComparator<
          ScalarType, IndexType, FakeComparator, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      FakeComparator comparator;
      Kokkos::MinMaxFirstLastLocCustomComparator reducer(view, comparator);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::FirstLoc<IndexType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::FirstLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected = Kokkos::LastLoc<IndexType, MemorySpace>;
      using ValueType           = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::LastLoc reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::StdIsPartitioned<IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::StdIsPartitioned reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }

    {
      using ReducerTypeExpected =
          Kokkos::StdPartitionPoint<IndexType, MemorySpace>;
      using ValueType = ReducerTypeExpected::value_type;
      Kokkos::View<ValueType, MemorySpace> view;
      Kokkos::StdPartitionPoint reducer(view);
      check_types<ReducerTypeExpected>(reducer);
    }
  }

  TestNestedReducerCTAD() {
    Kokkos::parallel_for(TeamPolicy(0, Kokkos::AUTO), *this);
  }
};

}  // namespace
