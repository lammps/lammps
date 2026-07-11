// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef TESTVIEWHOOKS_HPP_
#define TESTVIEWHOOKS_HPP_

#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

namespace Test {
template <class DeviceType, bool ExplicitSubscriber>
struct TestSubscriber {
  using test_view_type = std::conditional_t<
      ExplicitSubscriber,
      Kokkos::View<double **,
                   Kokkos::Experimental::SubscribableViewHooks<TestSubscriber>,
                   DeviceType>,
      Kokkos::View<double **, DeviceType>>;
  inline static test_view_type *self_ptr        = nullptr;
  inline static const test_view_type *other_ptr = nullptr;

  static void copy_constructed(test_view_type &self,
                               const test_view_type &other) {
    self_ptr  = &self;
    other_ptr = &other;
  }

  static void move_constructed(test_view_type &self,
                               const test_view_type &other) {
    self_ptr  = &self;
    other_ptr = &other;
  }

  static void copy_assigned(test_view_type &self, const test_view_type &other) {
    self_ptr  = &self;
    other_ptr = &other;
  }

  static void move_assigned(test_view_type &self, const test_view_type &other) {
    self_ptr  = &self;
    other_ptr = &other;
  }

  static void reset() {
    self_ptr  = nullptr;
    other_ptr = nullptr;
  }
};
}  // namespace Test

namespace TestADLViewHookCustomization {
template <typename DeviceType>
class TestMemSpace : public DeviceType::memory_space {
 public:
  using memory_space = TestMemSpace;
};

template <typename DeviceType>
using TestDevice = Kokkos::Device<typename DeviceType::execution_space,
                                  TestMemSpace<DeviceType>>;

template <class DataType, class... Properties>
constexpr auto customize_view_hooks() {
  using traits_type = Kokkos::ViewTraits<DataType, Properties...>;
  return Kokkos::Experimental::SubscribableViewHooks<
      Test::TestSubscriber<typename traits_type::device_type, false>>{};
}
}  // namespace TestADLViewHookCustomization

namespace Test {
template <class DeviceType>
struct TestViewHooks {
  static_assert(Kokkos::is_memory_space_v<
                TestADLViewHookCustomization::TestMemSpace<DeviceType>>);
  using subscriber_type = TestSubscriber<DeviceType, true>;
  using test_view_type  = typename subscriber_type::test_view_type;

  using adl_subscriber_type =
      TestSubscriber<TestADLViewHookCustomization::TestMemSpace<DeviceType>,
                     false>;
  using adl_test_view_type = typename adl_subscriber_type::test_view_type;

  static_assert(
      Kokkos::Experimental::is_hooks_policy<
          Kokkos::Experimental::SubscribableViewHooks<subscriber_type>>::value,
      "Must be a hooks policy");
  static_assert(
      std::same_as<typename test_view_type::type,
                   Kokkos::View<typename test_view_type::data_type,
                                typename test_view_type::array_layout,
                                typename test_view_type::device_type,
                                Kokkos::Experimental::SubscribableViewHooks<
                                    subscriber_type>,
                                typename test_view_type::memory_traits>>);

  static void testViewHooksCopyConstruct() {
    subscriber_type::reset();
    test_view_type testa;

    test_view_type testb(testa);
    EXPECT_EQ(subscriber_type::self_ptr, &testb);
    EXPECT_EQ(subscriber_type::other_ptr, &testa);
  }

  static void testViewHooksMoveConstruct() {
    subscriber_type::reset();
    test_view_type testa;

    test_view_type testb(std::move(testa));
    EXPECT_EQ(subscriber_type::self_ptr, &testb);

    // This is valid, even if the view is moved-from
    EXPECT_EQ(subscriber_type::other_ptr, &testa);
  }

  static void testViewHooksCopyAssign() {
    subscriber_type::reset();
    test_view_type testa;

    test_view_type testb;
    testb = testa;
    EXPECT_EQ(subscriber_type::self_ptr, &testb);
    EXPECT_EQ(subscriber_type::other_ptr, &testa);
  }

  static void testViewHooksMoveAssign() {
    subscriber_type::reset();
    test_view_type testa;

    test_view_type testb;
    testb = std::move(testa);
    EXPECT_EQ(subscriber_type::self_ptr, &testb);

    // This is valid, even if the view is moved-from
    EXPECT_EQ(subscriber_type::other_ptr, &testa);
  }

  static void testViewHooksCopyConstructWithADL() {
    adl_subscriber_type::reset();
    adl_test_view_type testa;

    adl_test_view_type testb(testa);
    EXPECT_EQ(adl_subscriber_type::self_ptr, &testb);
    EXPECT_EQ(adl_subscriber_type::other_ptr, &testa);
  }
};

TEST(TEST_CATEGORY, view_hooks) {
  using ExecSpace = TEST_EXECSPACE;
  TestViewHooks<ExecSpace>::testViewHooksCopyConstruct();
  TestViewHooks<ExecSpace>::testViewHooksMoveConstruct();
  TestViewHooks<ExecSpace>::testViewHooksCopyAssign();
  TestViewHooks<ExecSpace>::testViewHooksMoveAssign();
}

}  // namespace Test
#endif  // TESTVIEWHOOKS_HPP_
