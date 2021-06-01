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

#include <impl/Kokkos_Timer.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>
#include <TestNonTrivialScalarTypes.hpp>

namespace TestTeamVector {

template <typename Scalar, class ExecutionSpace>
struct functor_team_for {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_team_for(Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_int =
      Kokkos::View<Scalar *, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_int::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    using size_type           = typename shmem_space::size_type;
    const size_type shmemSize = team.team_size() * 13;
    shared_int values         = shared_int(team.team_shmem(), shmemSize);

    if (values.data() == nullptr ||
        static_cast<size_type>(values.extent(0)) < shmemSize) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "FAILED to allocate shared memory of size %u\n",
          static_cast<unsigned int>(shmemSize));
    } else {
      // Initialize shared memory.
      values(team.team_rank()) = 0;

      // Accumulate value into per thread shared memory.
      // This is non blocking.
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 131), [&](int i) {
        values(team.team_rank()) +=
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

        for (int i = 0; i < team.team_size(); ++i) {
          value += values(i);
        }

        if (test != value) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED team_parallel_for %i %i %f %f\n", team.league_rank(),
              team.team_rank(), static_cast<double>(test),
              static_cast<double>(value));
          flag() = 1;
        }
      });
    }
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_team_reduce {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_team_reduce(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_scalar_t =
      Kokkos::View<Scalar *, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_scalar_t::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Scalar value = Scalar();
    shared_scalar_t shared_value(team.team_scratch(0), 1);

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 131),
        [&](int i, Scalar &val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        value);

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 131),
        [&](int i, Scalar &val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        shared_value(0));

    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      Scalar test = 0;

      for (int i = 0; i < 131; ++i) {
        test += i - team.league_rank() + team.league_size() + team.team_size();
      }

      if (test != value) {
        if (team.league_rank() == 0) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED team_parallel_reduce %i %i %lf %lf %lu\n",
              team.league_rank(), team.team_rank(), static_cast<double>(test),
              static_cast<double>(value),
              static_cast<unsigned long>(sizeof(Scalar)));
        }

        flag() = 1;
      }
      if (test != shared_value(0)) {
        if (team.league_rank() == 0) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED team_parallel_reduce with shared result %i %i %lf %lf "
              "%lu\n",
              team.league_rank(), team.team_rank(), static_cast<double>(test),
              static_cast<double>(shared_value(0)),
              static_cast<unsigned long>(sizeof(Scalar)));
        }

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_team_reduce_reducer {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_team_reduce_reducer(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_scalar_t =
      Kokkos::View<Scalar *, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_scalar_t::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Scalar value = 0;
    shared_scalar_t shared_value(team.team_scratch(0), 1);

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 131),
        [&](int i, Scalar &val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        Kokkos::Sum<Scalar>(value));

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 131),
        [&](int i, Scalar &val) {
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
            "FAILED team_vector_parallel_reduce_reducer %i %i %lf %lf\n",
            team.league_rank(), team.team_rank(), static_cast<double>(test),
            static_cast<double>(value));

        flag() = 1;
      }
      if (test != shared_value(0)) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "FAILED team_vector_parallel_reduce_reducer shared value %i %i %lf "
            "%lf\n",
            team.league_rank(), team.team_rank(), static_cast<double>(test),
            static_cast<double>(shared_value(0)));

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_team_vector_for {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_team_vector_for(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_int =
      Kokkos::View<Scalar *, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_int::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    using size_type = typename shared_int::size_type;

    const size_type shmemSize = team.team_size() * 13;
    shared_int values         = shared_int(team.team_shmem(), shmemSize);

    if (values.data() == nullptr ||
        static_cast<size_type>(values.extent(0)) < shmemSize) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "FAILED to allocate shared memory of size %u\n",
          static_cast<unsigned int>(shmemSize));
    } else {
      team.team_barrier();

      Kokkos::single(Kokkos::PerThread(team),
                     [&]() { values(team.team_rank()) = 0; });

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 131), [&](int i) {
        Kokkos::single(Kokkos::PerThread(team), [&]() {
          values(team.team_rank()) +=
              i - team.league_rank() + team.league_size() + team.team_size();
        });
      });

      team.team_barrier();

      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        Scalar test  = 0;
        Scalar value = 0;

        for (int i = 0; i < 131; ++i) {
          test +=
              i - team.league_rank() + team.league_size() + team.team_size();
        }

        for (int i = 0; i < team.team_size(); ++i) {
          value += values(i);
        }

        if (test != value) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED team_vector_parallel_for %i %i %f %f\n",
              team.league_rank(), team.team_rank(), static_cast<double>(test),
              static_cast<double>(value));

          flag() = 1;
        }
      });
    }
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_team_vector_reduce {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;
  functor_team_vector_reduce(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_int =
      Kokkos::View<Scalar *, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_int::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Scalar value = Scalar();

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 131),
        [&](int i, Scalar &val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        value);

    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      Scalar test = 0;

      for (int i = 0; i < 131; ++i) {
        test += i - team.league_rank() + team.league_size() + team.team_size();
      }

      if (test != value) {
        if (team.league_rank() == 0) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF(
              "FAILED team_vector_parallel_reduce %i %i %f %f %lu\n",
              team.league_rank(), team.team_rank(), static_cast<double>(test),
              static_cast<double>(value),
              static_cast<unsigned long>(sizeof(Scalar)));
        }

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_team_vector_reduce_reducer {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_team_vector_reduce_reducer(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_int =
      Kokkos::View<Scalar *, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_int::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Scalar value = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, 131),
        [&](int i, Scalar &val) {
          val += i - team.league_rank() + team.league_size() + team.team_size();
        },
        Kokkos::Sum<Scalar>(value));

    team.team_barrier();

    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      Scalar test = 0;

      for (int i = 0; i < 131; ++i) {
        test += i - team.league_rank() + team.league_size() + team.team_size();
      }

      if (test != value) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "FAILED team_vector_parallel_reduce_reducer %i %i %f %f\n",
            team.league_rank(), team.team_rank(), static_cast<double>(test),
            static_cast<double>(value));

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_vec_single {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;
  int nStart;
  int nEnd;

  functor_vec_single(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_,
      const int start_, const int end_)
      : flag(flag_), nStart(start_), nEnd(end_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    // Warning: this test case intentionally violates permissible semantics.
    // It is not valid to get references to members of the enclosing region
    // inside a parallel_for and write to it.
    Scalar value = 0;

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, nStart, nEnd),
                         [&](int i) {
                           value = i;  // This write is violating Kokkos
                                       // semantics for nested parallelism.
                         });

    Kokkos::single(
        Kokkos::PerThread(team), [&](Scalar &val) { val = 1; }, value);

    Scalar value2 = 0;
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(team, nStart, nEnd),
        [&](int /*i*/, Scalar &val) { val += value; }, value2);

    if (value2 != (value * Scalar(nEnd - nStart))) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "FAILED vector_single broadcast %i %i %f %f\n", team.league_rank(),
          team.team_rank(), (double)value2, (double)value);

      flag() = 1;
    }
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_vec_for {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_vec_for(Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  using shmem_space = typename ExecutionSpace::scratch_memory_space;
  using shared_int =
      Kokkos::View<Scalar *, shmem_space, Kokkos::MemoryUnmanaged>;
  unsigned team_shmem_size(int team_size) const {
    return shared_int::shmem_size(team_size * 13);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    shared_int values = shared_int(team.team_shmem(), team.team_size() * 13);

    if (values.data() == nullptr ||
        values.extent(0) < (unsigned)team.team_size() * 13) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("FAILED to allocate memory of size %i\n",
                                    static_cast<int>(team.team_size() * 13));
      flag() = 1;
    } else {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, 13), [&](int i) {
        values(13 * team.team_rank() + i) =
            i - team.team_rank() - team.league_rank() + team.league_size() +
            team.team_size();
      });

      Kokkos::single(Kokkos::PerThread(team), [&]() {
        Scalar test  = 0;
        Scalar value = 0;

        for (int i = 0; i < 13; ++i) {
          test += i - team.team_rank() - team.league_rank() +
                  team.league_size() + team.team_size();
          value += values(13 * team.team_rank() + i);
        }

        if (test != value) {
          KOKKOS_IMPL_DO_NOT_USE_PRINTF("FAILED vector_par_for %i %i %f %f\n",
                                        team.league_rank(), team.team_rank(),
                                        static_cast<double>(test),
                                        static_cast<double>(value));

          flag() = 1;
        }
      });
    }
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_vec_red {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_vec_red(Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Scalar value = 0;

    // When no reducer is given the default is summation.
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(team, 13),
        [&](int i, Scalar &val) { val += i; }, value);

    Kokkos::single(Kokkos::PerThread(team), [&]() {
      Scalar test = 0;

      for (int i = 0; i < 13; i++) test += i;

      if (test != value) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF("FAILED vector_par_reduce %i %i %f %f\n",
                                      team.league_rank(), team.team_rank(),
                                      (double)test, (double)value);

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_vec_red_reducer {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;

  functor_vec_red_reducer(
      Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    // Must initialize to the identity value for the reduce operation
    // for this test:
    //   ( identity, operation ) = ( 1 , *= )
    Scalar value = 1;

    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(team, 13),
        [&](int i, Scalar &val) { val *= (i % 5 + 1); },
        Kokkos::Prod<Scalar>(value));

    Kokkos::single(Kokkos::PerThread(team), [&]() {
      Scalar test = 1;

      for (int i = 0; i < 13; i++) test *= (i % 5 + 1);

      if (test != value) {
        KOKKOS_IMPL_DO_NOT_USE_PRINTF(
            "FAILED vector_par_reduce_reducer %i %i %f %f\n",
            team.league_rank(), team.team_rank(), (double)test, (double)value);

        flag() = 1;
      }
    });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_vec_scan {
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;
  functor_vec_scan(Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team) const {
    Kokkos::parallel_scan(Kokkos::ThreadVectorRange(team, 13),
                          [&](int i, Scalar &val, bool final) {
                            val += i;

                            if (final) {
                              Scalar test = 0;
                              for (int k = 0; k <= i; k++) test += k;

                              if (test != val) {
                                KOKKOS_IMPL_DO_NOT_USE_PRINTF(
                                    "FAILED vector_par_scan %i %i %f %f\n",
                                    team.league_rank(), team.team_rank(),
                                    (double)test, (double)val);

                                flag() = 1;
                              }
                            }
                          });
  }
};

template <typename Scalar, class ExecutionSpace>
struct functor_reduce {
  using value_type      = double;
  using policy_type     = Kokkos::TeamPolicy<ExecutionSpace>;
  using execution_space = ExecutionSpace;

  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag;
  functor_reduce(Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> flag_)
      : flag(flag_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(typename policy_type::member_type team, double &sum) const {
    sum += team.league_rank() * 100 + team.thread_rank();
  }
};

template <typename Scalar, class ExecutionSpace>
bool test_scalar(int nteams, int team_size, int test) {
  Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace> d_flag("flag");
  typename Kokkos::View<int, Kokkos::LayoutLeft, ExecutionSpace>::HostMirror
      h_flag("h_flag");
  h_flag() = 0;
  Kokkos::deep_copy(d_flag, h_flag);

  if (test == 0) {
    Kokkos::parallel_for(
        std::string("A"),
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_vec_red<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 1) {
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_vec_red_reducer<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 2) {
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_vec_scan<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 3) {
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_vec_for<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 4) {
    Kokkos::parallel_for(
        "B", Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_vec_single<Scalar, ExecutionSpace>(d_flag, 0, 13));
  } else if (test == 5) {
    Kokkos::parallel_for(Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size),
                         functor_team_for<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 6) {
    Kokkos::parallel_for(Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size),
                         functor_team_reduce<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 7) {
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size),
        functor_team_reduce_reducer<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 8) {
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_team_vector_for<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 9) {
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_team_vector_reduce<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 10) {
    Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_team_vector_reduce_reducer<Scalar, ExecutionSpace>(d_flag));
  } else if (test == 11) {
    Kokkos::parallel_for(
        "B", Kokkos::TeamPolicy<ExecutionSpace>(nteams, team_size, 8),
        functor_vec_single<Scalar, ExecutionSpace>(d_flag, 4, 13));
  }

  Kokkos::deep_copy(h_flag, d_flag);

  return (h_flag() == 0);
}

template <class ExecutionSpace>
bool Test(int test) {
  bool passed = true;

  int team_size = 33;
  if (team_size > int(ExecutionSpace::concurrency()))
    team_size = int(ExecutionSpace::concurrency());
  passed = passed && test_scalar<int, ExecutionSpace>(317, team_size, test);
  passed = passed &&
           test_scalar<long long int, ExecutionSpace>(317, team_size, test);
  passed = passed && test_scalar<float, ExecutionSpace>(317, team_size, test);
  passed = passed && test_scalar<double, ExecutionSpace>(317, team_size, test);
  passed = passed &&
           test_scalar<Test::my_complex, ExecutionSpace>(317, team_size, test);
  passed = passed && test_scalar<Test::array_reduce<double, 1>, ExecutionSpace>(
                         317, team_size, test);
  passed = passed && test_scalar<Test::array_reduce<float, 1>, ExecutionSpace>(
                         317, team_size, test);
  passed = passed && test_scalar<Test::array_reduce<double, 3>, ExecutionSpace>(
                         317, team_size, test);

  return passed;
}

}  // namespace TestTeamVector

namespace Test {

// Computes y^T*A*x
// ( modified from kokkos-tutorials/GTC2016/Exercises/ThreeLevelPar )

#if (!defined(KOKKOS_ENABLE_CUDA)) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)
template <typename ScalarType, class DeviceType>
class TestTripleNestedReduce {
 public:
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  TestTripleNestedReduce(const size_type &nrows, const size_type &ncols,
                         const size_type &team_size,
                         const size_type &vector_length) {
    run_test(nrows, ncols, team_size, vector_length);
  }

  void run_test(const size_type &nrows, const size_type &ncols,
                size_type team_size, const size_type &vector_length) {
    if (team_size > size_type(DeviceType::execution_space::concurrency()))
      team_size = size_type(DeviceType::execution_space::concurrency());

#ifdef KOKKOS_ENABLE_HPX
    team_size = 1;
    if (!std::is_same<execution_space, Kokkos::Experimental::HPX>::value) {
      team_size = 1;
    }
#endif

    // using Layout = Kokkos::LayoutLeft;
    using Layout = Kokkos::LayoutRight;

    using ViewVector = Kokkos::View<ScalarType *, DeviceType>;
    using ViewMatrix = Kokkos::View<ScalarType **, Layout, DeviceType>;

    ViewVector y("y", nrows);
    ViewVector x("x", ncols);
    ViewMatrix A("A", nrows, ncols);

    using range_policy = Kokkos::RangePolicy<DeviceType>;

    // Initialize y vector.
    Kokkos::parallel_for(
        range_policy(0, nrows), KOKKOS_LAMBDA(const int i) { y(i) = 1; });

    // Initialize x vector.
    Kokkos::parallel_for(
        range_policy(0, ncols), KOKKOS_LAMBDA(const int i) { x(i) = 1; });
    Kokkos::fence();

    using team_policy = Kokkos::TeamPolicy<DeviceType>;
    using member_type = typename Kokkos::TeamPolicy<DeviceType>::member_type;

    // Initialize A matrix, note 2D indexing computation.
    Kokkos::parallel_for(
        team_policy(nrows, Kokkos::AUTO),
        KOKKOS_LAMBDA(const member_type &teamMember) {
          const int j = teamMember.league_rank();
          Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, ncols),
                               [&](const int i) { A(j, i) = 1; });
        });
    Kokkos::fence();

    // Three level parallelism kernel to force caching of vector x.
    ScalarType result = 0.0;
    int chunk_size    = 128;
    Kokkos::parallel_reduce(
        team_policy(nrows / chunk_size, team_size, vector_length),
        KOKKOS_LAMBDA(const member_type &teamMember, double &update) {
          const int row_start = teamMember.league_rank() * chunk_size;
          const int row_end   = row_start + chunk_size;
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(teamMember, row_start, row_end),
              [&](const int i) {
                ScalarType sum_i = 0.0;
                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange(teamMember, ncols),
                    [&](const int j, ScalarType &innerUpdate) {
                      innerUpdate += A(i, j) * x(j);
                    },
                    sum_i);
                Kokkos::single(Kokkos::PerThread(teamMember),
                               [&]() { update += y(i) * sum_i; });
              });
        },
        result);
    Kokkos::fence();

    const ScalarType solution = (ScalarType)nrows * (ScalarType)ncols;
    if (int64_t(solution) != int64_t(result)) {
      printf("  TestTripleNestedReduce failed solution(%" PRId64
             ") != result(%" PRId64
             "),"
             " nrows(%" PRId32 ") ncols(%" PRId32 ") league_size(%" PRId32
             ") team_size(%" PRId32 ")\n",
             int64_t(solution), int64_t(result), int32_t(nrows), int32_t(ncols),
             int32_t(nrows / chunk_size), int32_t(team_size));
    }

    ASSERT_EQ(solution, result);
  }
};

#else  // #if ( ! defined( KOKKOS_ENABLE_CUDA ) ) || defined(
       // KOKKOS_ENABLE_CUDA_LAMBDA )

template <typename ScalarType, class DeviceType>
class TestTripleNestedReduce {
 public:
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  TestTripleNestedReduce(const size_type &, const size_type, const size_type &,
                         const size_type) {}
};

#endif

namespace VectorScanReducer {
enum class ScanType : bool { Inclusive, Exclusive };

template <typename ExecutionSpace, ScanType scan_type, int n,
          int n_vector_range, class Reducer>
struct checkScan {
  const int n_team_thread_range = 1000;
  const int n_per_team          = n_team_thread_range * n_vector_range;

  using size_type  = typename ExecutionSpace::size_type;
  using value_type = typename Reducer::value_type;
  using view_type  = Kokkos::View<value_type[n], ExecutionSpace>;

  view_type inputs  = view_type{"inputs"};
  view_type outputs = view_type{"outputs"};

  value_type result;
  Reducer reducer = {result};

  struct ThreadVectorFunctor {
    KOKKOS_FUNCTION void operator()(const size_type j, value_type &update,
                                    const bool final) const {
      const size_type element = j + m_team_offset + m_thread_offset;
      const auto tmp          = m_inputs(element);
      if (scan_type == ScanType::Inclusive) {
        m_reducer.join(update, tmp);
        if (final) {
          m_outputs(element) = update;
        }
      } else {
        if (final) {
          m_outputs(element) = update;
        }
        m_reducer.join(update, tmp);
      }
    }

    const Reducer &m_reducer;
    const size_type &m_team_offset;
    const size_type &m_thread_offset;
    const view_type &m_outputs;
    const view_type &m_inputs;
  };

  struct TeamThreadRangeFunctor {
    KOKKOS_FUNCTION void operator()(const size_type i) const {
      const size_type thread_offset = i * n_vector_range;
      Kokkos::parallel_scan(
          Kokkos::ThreadVectorRange(m_team, n_vector_range),
          ThreadVectorFunctor{m_reducer, m_team_offset, thread_offset,
                              m_outputs, m_inputs},
          m_reducer);
    }

    const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type &m_team;
    const Reducer &m_reducer;
    const size_type &m_team_offset;
    const view_type &m_outputs;
    const view_type &m_inputs;
  };

  KOKKOS_FUNCTION void operator()(
      const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type &team)
      const {
    const size_type iTeam       = team.league_rank();
    const size_type iTeamOffset = iTeam * n_per_team;
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, n_team_thread_range),
        TeamThreadRangeFunctor{team, reducer, iTeamOffset, outputs, inputs});
  }

  KOKKOS_FUNCTION void operator()(size_type i) const { inputs(i) = i * 1. / n; }

  void run() {
    const int n_teams = n / n_per_team;

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0, n), *this);

    // run ThreadVectorRange parallel_scan
    Kokkos::TeamPolicy<ExecutionSpace> policy(n_teams, Kokkos::AUTO,
                                              Kokkos::AUTO);
    const std::string label =
        (scan_type == ScanType::Inclusive ? std::string("inclusive")
                                          : std::string("exclusive")) +
        "Scan" + typeid(Reducer).name();
    Kokkos::parallel_for(label, policy, *this);
    Kokkos::fence();

    auto host_outputs =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, outputs);
    auto host_inputs =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, inputs);

    Kokkos::View<value_type[n], Kokkos::HostSpace> expected("expected");
    {
      value_type identity;
      reducer.init(identity);
      for (int i = 0; i < expected.extent_int(0); ++i) {
        const int vector       = i % n_vector_range;
        const value_type accum = vector == 0 ? identity : expected(i - 1);
        const value_type val =
            scan_type == ScanType::Inclusive
                ? host_inputs(i)
                : (vector == 0 ? identity : host_inputs(i - 1));
        expected(i) = accum;
        reducer.join(expected(i), val);
      }
    }
    for (int i = 0; i < host_outputs.extent_int(0); ++i)
      ASSERT_EQ(host_outputs(i), expected(i));
  }
};
}  // namespace VectorScanReducer

#if !(defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND) || defined(KOKKOS_ENABLE_HIP))
TEST(TEST_CATEGORY, team_vector) {
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(0)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(1)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(2)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(3)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(4)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(5)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(6)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(7)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(8)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(9)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(10)));
  ASSERT_TRUE((TestTeamVector::Test<TEST_EXECSPACE>(11)));
}
#endif

#if !defined(KOKKOS_IMPL_CUDA_CLANG_WORKAROUND)
TEST(TEST_CATEGORY, triple_nested_parallelism) {
// With KOKKOS_ENABLE_DEBUG enabled, the functor uses too many registers to run
// with a team size of 32 on GPUs, 16 is the max possible (at least on a K80
// GPU) See https://github.com/kokkos/kokkos/issues/1513
#if defined(KOKKOS_ENABLE_DEBUG) && defined(KOKKOS_ENABLE_CUDA)
  if (!std::is_same<TEST_EXECSPACE, Kokkos::Cuda>::value) {
#endif
    TestTripleNestedReduce<double, TEST_EXECSPACE>(8192, 2048, 32, 32);
    TestTripleNestedReduce<double, TEST_EXECSPACE>(8192, 2048, 32, 16);
#if defined(KOKKOS_ENABLE_DEBUG) && defined(KOKKOS_ENABLE_CUDA)
  }
#endif
  TestTripleNestedReduce<double, TEST_EXECSPACE>(8192, 2048, 16, 16);
  TestTripleNestedReduce<double, TEST_EXECSPACE>(8192, 2048, 16, 33);
  TestTripleNestedReduce<double, TEST_EXECSPACE>(8192, 2048, 16, 19);
  TestTripleNestedReduce<double, TEST_EXECSPACE>(8192, 2048, 7, 16);
}
#endif

TEST(TEST_CATEGORY, parallel_scan_with_reducers) {
  using T = double;
  using namespace VectorScanReducer;

  static constexpr int n              = 1000000;
  static constexpr int n_vector_range = 100;

  checkScan<TEST_EXECSPACE, ScanType::Exclusive, n, n_vector_range,
            Kokkos::Prod<T, TEST_EXECSPACE>>()
      .run();
  checkScan<TEST_EXECSPACE, ScanType::Inclusive, n, n_vector_range,
            Kokkos::Prod<T, TEST_EXECSPACE>>()
      .run();

  checkScan<TEST_EXECSPACE, ScanType::Exclusive, n, n_vector_range,
            Kokkos::Max<T, TEST_EXECSPACE>>()
      .run();
  checkScan<TEST_EXECSPACE, ScanType::Inclusive, n, n_vector_range,
            Kokkos::Max<T, TEST_EXECSPACE>>()
      .run();

  checkScan<TEST_EXECSPACE, ScanType::Exclusive, n, n_vector_range,
            Kokkos::Min<T, TEST_EXECSPACE>>()
      .run();
  checkScan<TEST_EXECSPACE, ScanType::Inclusive, n, n_vector_range,
            Kokkos::Min<T, TEST_EXECSPACE>>()
      .run();
}

}  // namespace Test
