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

struct TestTeamThreadMDRangeCTAD {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamHandle = TeamPolicy::member_type;

  KOKKOS_FUNCTION void operator()(TeamHandle const& team_handle) const {
    {
      Kokkos::TeamThreadMDRange md_range(team_handle, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamThreadMDRange<Kokkos::Rank<2>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamThreadMDRange md_range(team_handle, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamThreadMDRange<Kokkos::Rank<3>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamThreadMDRange md_range(team_handle, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamThreadMDRange<Kokkos::Rank<4>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamThreadMDRange md_range(team_handle, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamThreadMDRange<Kokkos::Rank<5>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamThreadMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamThreadMDRange<Kokkos::Rank<6>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamThreadMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamThreadMDRange<Kokkos::Rank<7>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamThreadMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamThreadMDRange<Kokkos::Rank<8>, TeamHandle>,
                         decltype(md_range)>);
    }
  }

  TestTeamThreadMDRangeCTAD() {
    Kokkos::parallel_for(TeamPolicy(0, Kokkos::AUTO), *this);
  }
};

struct TestTeamVectorMDRangeCTAD {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamHandle = TeamPolicy::member_type;

  KOKKOS_FUNCTION void operator()(TeamHandle const& team_handle) const {
    {
      Kokkos::TeamVectorMDRange md_range(team_handle, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamVectorMDRange<Kokkos::Rank<2>, TeamHandle>,
                         decltype(md_range)>);
    }
    {
      Kokkos::TeamVectorMDRange md_range(team_handle, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamVectorMDRange<Kokkos::Rank<3>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamVectorMDRange md_range(team_handle, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamVectorMDRange<Kokkos::Rank<4>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamVectorMDRange<Kokkos::Rank<5>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamVectorMDRange<Kokkos::Rank<6>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamVectorMDRange<Kokkos::Rank<7>, TeamHandle>,
                         decltype(md_range)>);
    }

    {
      Kokkos::TeamVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0, 0, 0);
      static_assert(
          std::is_same_v<Kokkos::TeamVectorMDRange<Kokkos::Rank<8>, TeamHandle>,
                         decltype(md_range)>);
    }
  }

  TestTeamVectorMDRangeCTAD() {
    Kokkos::parallel_for(TeamPolicy(0, Kokkos::AUTO), *this);
  }
};

struct TestThreadVectorMDRangeCTAD {
  using TeamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>;
  using TeamHandle = TeamPolicy::member_type;

  template <class PolicyTypeExpected, class PolicyTypeToCheck>
  KOKKOS_FUNCTION static void check_types([
      [maybe_unused]] PolicyTypeToCheck const& team_handle) {
    static_assert(std::is_same_v<PolicyTypeExpected, PolicyTypeToCheck>);
  }

  KOKKOS_FUNCTION void operator()(TeamHandle const& team_handle) const {
    {
      Kokkos::ThreadVectorMDRange md_range(team_handle, 0, 0);
      check_types<Kokkos::ThreadVectorMDRange<Kokkos::Rank<2>, TeamHandle>>(
          md_range);
    }

    {
      Kokkos::ThreadVectorMDRange md_range(team_handle, 0, 0, 0);
      check_types<Kokkos::ThreadVectorMDRange<Kokkos::Rank<3>, TeamHandle>>(
          md_range);
    }

    {
      Kokkos::ThreadVectorMDRange md_range(team_handle, 0, 0, 0, 0);
      check_types<Kokkos::ThreadVectorMDRange<Kokkos::Rank<4>, TeamHandle>>(
          md_range);
    }

    {
      Kokkos::ThreadVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0);
      check_types<Kokkos::ThreadVectorMDRange<Kokkos::Rank<5>, TeamHandle>>(
          md_range);
    }

    {
      Kokkos::ThreadVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0);
      check_types<Kokkos::ThreadVectorMDRange<Kokkos::Rank<6>, TeamHandle>>(
          md_range);
    }

    {
      Kokkos::ThreadVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0, 0);
      check_types<Kokkos::ThreadVectorMDRange<Kokkos::Rank<7>, TeamHandle>>(
          md_range);
    }

    {
      Kokkos::ThreadVectorMDRange md_range(team_handle, 0, 0, 0, 0, 0, 0, 0, 0);
      check_types<Kokkos::ThreadVectorMDRange<Kokkos::Rank<8>, TeamHandle>>(
          md_range);
    }
  }

  TestThreadVectorMDRangeCTAD() {
    Kokkos::parallel_for(TeamPolicy(0, Kokkos::AUTO), *this);
  }
};

}  // namespace
