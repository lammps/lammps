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

#include <Kokkos_Core.hpp>

#include <Kokkos_Timer.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>

namespace TestTeamVectorRange {

struct my_complex {
  double re, im;
  int dummy;

  KOKKOS_INLINE_FUNCTION
  my_complex() {
    re    = 0.0;
    im    = 0.0;
    dummy = 0;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex(const my_complex& src) {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex& operator=(const my_complex& src) {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex(const double& val) {
    re    = val;
    im    = 0.0;
    dummy = 0;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex& operator+=(const my_complex& src) {
    re += src.re;
    im += src.im;
    dummy += src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex operator+(const my_complex& src) {
    my_complex tmp = *this;
    tmp.re += src.re;
    tmp.im += src.im;
    tmp.dummy += src.dummy;
    return tmp;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex& operator*=(const my_complex& src) {
    double re_tmp = re * src.re - im * src.im;
    double im_tmp = re * src.im + im * src.re;
    re            = re_tmp;
    im            = im_tmp;
    dummy *= src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(const my_complex& src) const {
    return (re == src.re) && (im == src.im) && (dummy == src.dummy);
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const my_complex& src) const {
    return (re != src.re) || (im != src.im) || (dummy != src.dummy);
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const double& val) const {
    return (re != val) || (im != 0) || (dummy != 0);
  }

  KOKKOS_INLINE_FUNCTION
  my_complex& operator=(const int& val) {
    re    = val;
    im    = 0.0;
    dummy = 0;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex& operator=(const double& val) {
    re    = val;
    im    = 0.0;
    dummy = 0;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  operator double() { return re; }
};
}  // namespace TestTeamVectorRange

namespace Kokkos {
template <>
struct reduction_identity<TestTeamVectorRange::my_complex> {
  using t_red_ident = reduction_identity<double>;
  KOKKOS_FORCEINLINE_FUNCTION static TestTeamVectorRange::my_complex sum() {
    return TestTeamVectorRange::my_complex(t_red_ident::sum());
  }
  KOKKOS_FORCEINLINE_FUNCTION static TestTeamVectorRange::my_complex prod() {
    return TestTeamVectorRange::my_complex(t_red_ident::prod());
  }
};
}  // namespace Kokkos

namespace TestTeamVectorRange {

template <typename Scalar, class ExecutionSpace>
struct functor_teamvector_for {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_teamvector_for(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_int =
      Kokkos::View<Scalar*, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int /*team_size*/) const {
    return shared_int::shmem_size(131);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    using size_type           = typename shmem_space::size_type;
    const size_type shmemSize = 131;
    shared_int values         = shared_int(team.team_shmem(), shmemSize);

    if (values.data() == nullptr || values.extent(0) < shmemSize) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "FAILED to allocate shared memory of size %u\n",
          static_cast<unsigned int>(shmemSize));
    } else {
      // Initialize shared memory.
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 131),
                           [&](int i) { values(i) = 0; });
      // Wait for all memory to be written.
      team.team_barrier();

      // Accumulate value into per thread shared memory.
      // This is non blocking.
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 131), [&](int i) {
        values(i) +=
            i - team.league_rank() + team.league_size() + team.team_size();
      });

      // Wait for all memory to be written.
      team.team_barrier();

      // One thread per team executes the comparison.
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Scalar test  = 0;
        Scalar value = 0;

        for (int i = 0; i < 131; ++i) {
          test +=
              i - team.league_rank() + team.league_size() + team.team_size();
        }

        for (int i = 0; i < 131; ++i) {
          value += values(i);
        }

        if (test != value) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED teamvector_parallel_for %i %i %lf %lf\n",
              team.league_rank(), team.team_rank(), static_cast<double>(test),
              static_cast<double>(value));
          flag() = 1;
        }
      });
    }
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_teamvector_reduce {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_teamvector_reduce(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_scalar_t =
      Kokkos::View<Scalar*, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_scalar_t::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Scalar value = Scalar();
    shared_scalar_t shared_value(team.team_scratch(0), 1);

    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, 131),
        [&](int i, Scalar& val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        shared_value(0));

    team.team_barrier();
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, 131),
        [&](int i, Scalar& val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        value);

    //    Kokkos::parallel_reduce( Kokkos::TeamVectorRange( team, 131 ), [&] (
    //    int i, Scalar & val )
    //    {
    //      val += i - team.league_rank() + team.league_size() +
    //      team.team_size();
    //    }, shared_value(0) );

    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      Scalar test = 0;

      for (int i = 0; i < 131; ++i) {
        test += i - team.league_rank() + team.league_size() + team.team_size();
      }

      if (test != value) {
        if (team.league_rank() == 0) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED teamvector_parallel_reduce %i %i %lf %lf %lu\n",
              (int)team.league_rank(), (int)team.team_rank(),
              static_cast<double>(test), static_cast<double>(value),
              static_cast<unsigned long>(sizeof(Scalar)));
        }

        flag() = 1;
      }
      if (test != shared_value(0)) {
        if (team.league_rank() == 0) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED teamvector_parallel_reduce with shared result %i %i %lf "
              "%lf %lu\n",
              static_cast<int>(team.league_rank()),
              static_cast<int>(team.team_rank()), static_cast<double>(test),
              static_cast<double>(shared_value(0)),
              static_cast<unsigned long>(sizeof(Scalar)));
        }

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_teamvector_reduce_reducer {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_teamvector_reduce_reducer(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_scalar_t =
      Kokkos::View<Scalar*, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_scalar_t::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Scalar value = 0;
    shared_scalar_t shared_value(team.team_scratch(0), 1);

    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, 131),
        [&](int i, Scalar& val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        Kokkos::Sum<Scalar>(value));

    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, 131),
        [&](int i, Scalar& val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        Kokkos::Sum<Scalar>(shared_value(0)));

    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      Scalar test = 0;

      for (int i = 0; i < 131; ++i) {
        test += i - team.league_rank() + team.league_size() + team.team_size();
      }

      if (test != value) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "FAILED teamvector_parallel_reduce_reducer %i %i %lf %lf\n",
            team.league_rank(), team.team_rank(), static_cast<double>(test),
            static_cast<double>(value));

        flag() = 1;
      }
      if (test != shared_value(0)) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "FAILED teamvector_parallel_reduce_reducer shared value %i %i %lf "
            "%lf\n",
            team.league_rank(), team.team_rank(), static_cast<double>(test),
            static_cast<double>(shared_value(0)));

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
bool test_scalar(int nteams, int team_size, int test) {
  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> d_flag("flag");
  typename Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace>::HostMirror
      h_flag("h_flag");
  h_flag() = 0;
  Kokkos::deep_copy(d_flag, h_flag);

  Kokkos::TeamPolicy<ExecutionSpace> policy(nteams, team_size, 8);

  // FIXME_OPENMPTARGET - Need to allocate scratch space via set_scratch_space
  // for the OPENMPTARGET backend.
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  using scratch_t = Kokkos::View<Scalar*, ExecutionSpace,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  int scratch_size = 0;
  if (test == 0) {
    scratch_size = scratch_t::shmem_size(131);
  } else {
    // FIXME_OPENMPTARGET - Currently allocating more than one team for nested
    // reduction leads to runtime errors of illegal memory access, caused mostly
    // due to the OpenMP memory allocation constraints.
    policy       = Kokkos::TeamPolicy<ExecutionSpace>(1, team_size, 8);
    scratch_size = scratch_t::shmem_size(1);
  }

  policy.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
#endif

  if (test == 0) {
    Kokkos::parallel_for(
        "Test::TeamVectorFor", policy,
        functor_teamvector_for<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 1) {
    Kokkos::parallel_for(
        "Test::TeamVectorReduce", policy,
        functor_teamvector_reduce<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 2) {
    Kokkos::parallel_for(
        "Test::TeamVectorReduceReducer",
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_teamvector_reduce_reducer<Scalar, ExecutionSpace>(d_flag));
  }

  Kokkos::deep_copy(h_flag, d_flag);

  return (h_flag() == 0);
}

template <class ExecutionSpace>
bool Test(int test) {
  bool passed = true;

// With SYCL 33*8 exceeds the maximum work group size
#ifdef KOKKOS_ENABLE_SYCL
  int team_size = 31;
#else
  int team_size = 33;
#endif
  if (team_size > int(ExecutionSpace::concurrency()))
    team_size = int(ExecutionSpace::concurrency());
  passed = passed && test_scalar<int, ExecutionSpace>(317, team_size, test);
  passed = passed &&
           test_scalar<long long int, ExecutionSpace>(317, team_size, test);
  passed = passed && test_scalar<float, ExecutionSpace>(317, team_size, test);
  passed = passed && test_scalar<double, ExecutionSpace>(317, team_size, test);
  // FIXME_OPENMPTARGET - Use of custom reducers currently results in runtime
  // memory errors.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
  passed =
      passed && test_scalar<my_complex, ExecutionSpace>(317, team_size, test);
#endif

  return passed;
}

}  // namespace TestTeamVectorRange

namespace Test {

TEST(TEST_CATEGORY, team_teamvector_range) {
  ASSERT_TRUE((TestTeamVectorRange::Test<TEST_EXECSPACE>(0)));
  ASSERT_TRUE((TestTeamVectorRange::Test<TEST_EXECSPACE>(1)));
  // FIXME_OPENMPTARGET - Use of kokkos reducers currently results in runtime
  // memory errors.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
  ASSERT_TRUE((TestTeamVectorRange::Test<TEST_EXECSPACE>(2)));
#endif
}
}  // namespace Test
