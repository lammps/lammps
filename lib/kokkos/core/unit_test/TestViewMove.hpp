// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

namespace {

template <class ViewType>
void test_moving_view_use_count_and_label(ViewType v) {
  auto const ptr = v.data();
  auto const cnt = v.use_count();
  auto const lbl = v.label();

  // NOLINTBEGIN(bugprone-use-after-move)

  ViewType w(std::move(v));  // move construction
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  if (w.use_count() != 0)
    EXPECT_EQ(w.use_count(), cnt + 1);
  else
    EXPECT_EQ(w.use_count(), 0);
#else
  EXPECT_EQ(w.use_count(), cnt);
#endif
  EXPECT_EQ(w.data(), ptr);
  EXPECT_EQ(w.label(), lbl);
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  EXPECT_EQ(v.use_count(), w.use_count());
  EXPECT_EQ(v.label(), lbl);
#else
  EXPECT_EQ(v.use_count(), 0);
  EXPECT_EQ(v.label(), std::string(""));
#endif
  EXPECT_EQ(v.data(), ptr);

  v = std::move(w);  // move assignment
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  if (v.use_count() != 0)
    EXPECT_EQ(v.use_count(), cnt + 1);
  else
    EXPECT_EQ(w.use_count(), 0);
#else
  EXPECT_EQ(v.use_count(), cnt);
#endif
  EXPECT_EQ(v.data(), ptr);
  EXPECT_EQ(v.label(), lbl);
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  EXPECT_EQ(w.use_count(), v.use_count());
#else
  EXPECT_EQ(w.use_count(), 0);
#endif
  EXPECT_EQ(w.data(), v.data());
  EXPECT_EQ(v.label(), lbl);

  // NOLINTEND(bugprone-use-after-move)
}

TEST(TEST_CATEGORY, view_move_use_count_and_label) {
  using ExecutionSpace = TEST_EXECSPACE;

  test_moving_view_use_count_and_label(Kokkos::View<int, ExecutionSpace>("v0"));

  test_moving_view_use_count_and_label(
      Kokkos::View<float*, ExecutionSpace>("v1", 1));

  Kokkos::View<double**, ExecutionSpace> v2("v2", 1, 2);
  test_moving_view_use_count_and_label(Kokkos::View<double**, ExecutionSpace>(
      v2.data(), v2.extent(0), v2.extent(1)));
  test_moving_view_use_count_and_label(
      Kokkos::View<double**, ExecutionSpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
          v2.data(), v2.extent(0), v2.extent(1)));
}

template <class ViewType>
KOKKOS_FUNCTION int check_moved_from_view_state(ViewType v) {
  int err = 0;

  constexpr int rank = ViewType::rank();
  std::size_t
      exts[rank + 1];  // add a dummy trailing element to accomodate rank zero
  for (int i = 0; i < rank; ++i) {
    exts[i] = v.extent(i);
  }
  auto* const ptr = v.data();
  auto const span = v.span();

#define CHECK(CONTEXT, VIEW, PTR, SPAN, EXTS)                              \
  for (int i = 0; i < rank; ++i) {                                         \
    if (VIEW.extent(i) != EXTS[i]) {                                       \
      Kokkos::printf(CONTEXT "expected equality of " #VIEW                 \
                             ".extent(%d) which is %d and " #EXTS          \
                             "[%d] which is %d\n",                         \
                     i, VIEW.extent(i), i, EXTS[i]);                       \
      ++err;                                                               \
    }                                                                      \
  }                                                                        \
  if (VIEW.span() != SPAN) {                                               \
    Kokkos::printf(CONTEXT "expected equality of " #VIEW                   \
                           ".span() which is %d and " #SPAN "which is %d", \
                   VIEW.span(), span);                                     \
    ++err;                                                                 \
  }                                                                        \
  if (VIEW.data() != PTR) {                                                \
    Kokkos::printf(CONTEXT "expected equality of " #VIEW                   \
                           ".data() which is %p and " #PTR "which is %p",  \
                   VIEW.data(), ptr);                                      \
    ++err;                                                                 \
  }

  // NOLINTBEGIN(bugprone-use-after-move)

  ViewType w(std::move(v));  // move construction

  CHECK("failed moved-from view after calling move constructor\n", v, ptr, span,
        exts)
  CHECK("failed moved-from view after calling move constructor\n", w, ptr, span,
        exts)

  v = std::move(w);  // move assignment

  CHECK("failed moved-from view after calling move assignment operator\n", v,
        ptr, span, exts)
  CHECK("failed moved-from view after calling move assignment operator\n", w,
        ptr, span, exts)

  // NOLINTEND(bugprone-use-after-move)

  return err;
}

template <class ViewType>
void test_moved_from_view(ViewType v) {
  using ExecutionSpace = typename ViewType::execution_space;
  int errors;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int, int& err) { err += check_moved_from_view_state(v); },
      errors);
  EXPECT_EQ(errors, 0) << "within parallel region";
}

TEST(TEST_CATEGORY, view_moved_from) {
  using ExecutionSpace = TEST_EXECSPACE;

  test_moved_from_view(Kokkos::View<int, ExecutionSpace>("v0"));
  test_moved_from_view(Kokkos::View<float*, ExecutionSpace>("v1", 1));
  Kokkos::View<double**, ExecutionSpace> v2("v2", 1, 2);
  test_moved_from_view(Kokkos::View<double**, ExecutionSpace>(
      v2.data(), v2.extent(0), v2.extent(1)));
  test_moved_from_view(Kokkos::View<double**, ExecutionSpace,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
      v2.data(), v2.extent(0), v2.extent(1)));
}

}  // namespace
