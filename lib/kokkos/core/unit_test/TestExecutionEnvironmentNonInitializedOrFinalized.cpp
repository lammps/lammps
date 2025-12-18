// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <cstdlib>
#include <type_traits>

#include "KokkosExecutionEnvironmentNeverInitializedFixture.hpp"

namespace {

using ExecutionEnvironmentNonInitializedOrFinalized_DeathTest =
    KokkosExecutionEnvironmentNeverInitialized;

struct NonTrivial {
  KOKKOS_FUNCTION NonTrivial() {}
};
static_assert(!std::is_trivially_default_constructible_v<NonTrivial>);

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest,
       default_constructed_views) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  auto make_views = [] {
    Kokkos::View<int> v0;
    Kokkos::View<float*> v1;
    Kokkos::View<NonTrivial**> v2;
    return std::make_tuple(v0, v1, v2);
  };
  EXPECT_EXIT(
      {
        { auto views = make_views(); }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_EXIT(
      {
        {
          Kokkos::initialize();
          auto views =
              make_views();  // views outlive the Kokkos execution environment
          Kokkos::finalize();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_EXIT(
      {
        {
          Kokkos::initialize();
          Kokkos::finalize();
          auto views = make_views();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
}

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest, views) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_EXIT(
      {
        {
          Kokkos::View<int*> v;
          Kokkos::initialize();
          v = Kokkos::View<int*>("v", 10);
          v = Kokkos::View<int*>();
          Kokkos::finalize();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_EXIT(
      {
        {
          Kokkos::initialize();
          Kokkos::View<int*> v("v", 10);
          v = {};  // assign default constructed view
          Kokkos::finalize();
        }
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::View<int*> v("v", 0);
        Kokkos::finalize();
      },
      "Kokkos allocation \"v\" is being deallocated after Kokkos::finalize was "
      "called");

  EXPECT_DEATH(
      { Kokkos::View<int*> v("v", 0); },
      "View \\(label=\\\"v\\\"\\) is being constructed before initialize");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        Kokkos::View<int*> v("v", 0);
      },
      "View \\(label=\\\"v\\\"\\) is being constructed after finalize");
}

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest,
       c_style_memory_management) {
// FIXME_THREADS: Checking for calls to kokkos_malloc, kokkos_realloc,
// kokkos_free before initialize or after finalize is currently disabled
// for the Threads backend. Refer issue #7944.
#ifdef KOKKOS_ENABLE_THREADS
  GTEST_SKIP()
      << "skipping since initializing Threads backend calls kokkos_malloc";
#endif
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_DEATH(
      { [[maybe_unused]] void* ptr = Kokkos::kokkos_malloc(1); },
      "Kokkos ERROR: attempting to perform C-style memory management via "
      "kokkos_malloc\\(\\) \\*\\*before\\*\\* Kokkos::initialize\\(\\) was "
      "called");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        [[maybe_unused]] void* ptr = Kokkos::kokkos_malloc(1);
      },
      "Kokkos ERROR: attempting to perform C-style memory management via "
      "kokkos_malloc\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        void* ptr = Kokkos::kokkos_malloc(1);
        Kokkos::finalize();
        Kokkos::kokkos_free(ptr);
      },
      "Kokkos ERROR: attempting to perform C-style memory management via "
      "kokkos_free\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was called");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        void* prev = Kokkos::kokkos_malloc(1);
        Kokkos::finalize();
        [[maybe_unused]] void* next = Kokkos::kokkos_realloc(prev, 2);
      },
      "Kokkos ERROR: attempting to perform C-style memory management via "
      "kokkos_realloc\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called");
  EXPECT_DEATH(
      {
        // Take a fake pointer
        void* ptr = reinterpret_cast<void*>(0x8BADF00D);
        Kokkos::kokkos_free(ptr);
      },
      "Kokkos ERROR: attempting to perform C-style memory management via "
      "kokkos_free\\(\\) \\*\\*before\\*\\* Kokkos::initialize\\(\\) was "
      "called");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        // Take a fake pointer
        void* ptr = reinterpret_cast<void*>(0xB105F00D);
        Kokkos::kokkos_free(ptr);
      },
      "Kokkos ERROR: attempting to perform C-style memory management via "
      "kokkos_free\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was called");
}

template <typename Type>
class EmptyReduceFunctor {
 public:
  using size_type = typename Kokkos::DefaultExecutionSpace::size_type;
  KOKKOS_INLINE_FUNCTION
  void join(Type&, const Type&) const {}
  KOKKOS_INLINE_FUNCTION
  void operator()(size_type, Type&) const {}
  KOKKOS_INLINE_FUNCTION
  void final(Type&) const {}
};

struct ForFunctor {
  KOKKOS_FUNCTION void operator()(int) const {}
};

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest, parallel_for) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::RangePolicy<> policy(0, 0);
        Kokkos::finalize();

        Kokkos::parallel_for("for_policy", policy, ForFunctor());
      },
      "Kokkos ERROR: attempting to call "
      "parallel_for\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was called. "
      "Concerns for_policy with exec policy RangePolicy.");
  EXPECT_DEATH(
      { Kokkos::parallel_for("for_workcnt", 0, ForFunctor()); },
      "Kokkos ERROR: attempting to call "
      "parallel_for\\(\\) \\*\\*before\\*\\* Kokkos::initialize\\(\\) was "
      "called. Concerns for_workcnt with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        Kokkos::parallel_for("for_workcnt", 0, ForFunctor());
      },
      "Kokkos ERROR: attempting to call "
      "parallel_for\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was called. "
      "Concerns for_workcnt with exec policy RangePolicy.");
}

struct ReduceFunctor {
  KOKKOS_FUNCTION void operator()(int, float&) const {}
};

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest,
       parallel_reduce) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_DEATH(
      {
        using functor_type = EmptyReduceFunctor<float>;
        Kokkos::initialize();
        Kokkos::RangePolicy<> policy(0, 0);
        Kokkos::finalize();

        Kokkos::parallel_reduce("reduce_policy_func_noret", policy,
                                functor_type{});
      },
      "Kokkos ERROR: attempting to call parallel_reduce\\(\\) "
      "\\*\\*after\\*\\* Kokkos::finalize\\(\\) was called. "
      "Concerns reduce_policy_func_noret with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        using functor_type = EmptyReduceFunctor<float>;
        Kokkos::parallel_reduce("reduce_workcnt_func_noret", 0, functor_type{});
      },
      "Kokkos ERROR: attempting to call parallel_reduce\\(\\) "
      "\\*\\*before\\*\\* Kokkos::initialize\\(\\) was called. "
      "Concerns reduce_workcnt_func_noret with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        using functor_type = EmptyReduceFunctor<float>;
        Kokkos::parallel_reduce("reduce_workcnt_func_noret", 0, functor_type{});
      },
      "Kokkos ERROR: attempting to call parallel_reduce\\(\\) "
      "\\*\\*after\\*\\* Kokkos::finalize\\(\\) was called. "
      "Concerns reduce_workcnt_func_noret with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        float x;
        Kokkos::initialize();
        Kokkos::RangePolicy<> policy(0, 0);
        Kokkos::finalize();
        Kokkos::parallel_reduce("reduce_policy_ret_float", policy,
                                ReduceFunctor(), x);
      },
      "Kokkos ERROR: attempting to call "
      "parallel_reduce\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called. Concerns reduce_policy_ret_float with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        float x;
        Kokkos::parallel_reduce("reduce_workcnt_ret_float", 0, ReduceFunctor(),
                                x);
      },
      "Kokkos ERROR: attempting to call parallel_reduce\\(\\) "
      "\\*\\*before\\*\\* Kokkos::initialize\\(\\) was called. "
      "Concerns reduce_workcnt_ret_float with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        float x;
        Kokkos::parallel_reduce("reduce_workcnt_ret_float", 0, ReduceFunctor(),
                                x);
      },
      "Kokkos ERROR: attempting to call parallel_reduce\\(\\) "
      "\\*\\*after\\*\\* Kokkos::finalize\\(\\) was called. "
      "Concerns reduce_workcnt_ret_float with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::View<float> x{"x"};
        Kokkos::RangePolicy<> policy(0, 0);
        Kokkos::finalize();
        Kokkos::parallel_reduce("reduce_policy_ret_view", policy,
                                ReduceFunctor(), x);
      },
      "Kokkos ERROR: attempting to call "
      "parallel_reduce\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called. Concerns reduce_policy_ret_view with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::View<float> x{"x"};
        Kokkos::finalize();
        Kokkos::parallel_reduce("reduce_policy_ret_view", 0, ReduceFunctor(),
                                x);
      },
      "Kokkos ERROR: attempting to call "
      "parallel_reduce\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called. Concerns reduce_policy_ret_view with exec policy RangePolicy.");
}

struct ScanFunctor {
  KOKKOS_FUNCTION void operator()(int, float&, bool) const {}
};

TEST_F(ExecutionEnvironmentNonInitializedOrFinalized_DeathTest, parallel_scan) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::RangePolicy<> policy(0, 0);
        Kokkos::finalize();
        Kokkos::parallel_scan("scan_policy_noret", policy, ScanFunctor());
      },
      "Kokkos ERROR: attempting to call "
      "parallel_scan\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called. Concerns scan_policy_noret with exec policy RangePolicy.");
  EXPECT_DEATH(
      { Kokkos::parallel_scan("scan_workcnt_noret", 0, ScanFunctor()); },
      "Kokkos ERROR: attempting to call "
      "parallel_scan\\(\\) \\*\\*before\\*\\* Kokkos::initialize\\(\\) was "
      "called. Concerns scan_workcnt_noret with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        Kokkos::parallel_scan("scan_workcnt_noret", 0, ScanFunctor());
      },
      "Kokkos ERROR: attempting to call "
      "parallel_scan\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called. Concerns scan_workcnt_noret with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        float x;
        Kokkos::initialize();
        Kokkos::RangePolicy<> policy(0, 0);
        Kokkos::finalize();
        Kokkos::parallel_scan("scan_policy_ret_float", policy, ScanFunctor(),
                              x);
      },
      "Kokkos ERROR: attempting to call "
      "parallel_scan\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called. Concerns scan_policy_ret_float with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        float x;
        Kokkos::parallel_scan("scan_workcnt_ret_float", 0, ScanFunctor(), x);
      },
      "Kokkos ERROR: attempting to call "
      "parallel_scan\\(\\) \\*\\*before\\*\\* Kokkos::initialize\\(\\) was "
      "called. Concerns scan_workcnt_ret_float with exec policy RangePolicy.");
  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        float x;
        Kokkos::parallel_scan("scan_workcnt_ret_float", 0, ScanFunctor(), x);
      },
      "Kokkos ERROR: attempting to call "
      "parallel_scan\\(\\) \\*\\*after\\*\\* Kokkos::finalize\\(\\) was "
      "called. Concerns scan_workcnt_ret_float with exec policy RangePolicy.");
}
}  // namespace
