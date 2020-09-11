#include <Kokkos_Core.hpp>

namespace Test {

#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
namespace Impl {
template <class MemorySpaceA, class MemorySpaceB>
struct TestDeepCopy {
  using a_base_t = Kokkos::View<double*, Kokkos::LayoutRight, MemorySpaceA>;
  using b_base_t = Kokkos::View<double*, Kokkos::LayoutRight, MemorySpaceB>;
  using a_char_t = Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA>;
  using b_char_t = Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceB>;

  using policyA_t = Kokkos::RangePolicy<typename MemorySpaceA::execution_space>;
  using policyB_t = Kokkos::RangePolicy<typename MemorySpaceB::execution_space>;

  static void reset_a_copy_and_b(
      Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char_copy,
      Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceB> b_char) {
    const int N = b_char.extent(0);
    Kokkos::parallel_for(
        "TestDeepCopy: FillA_copy", policyA_t(0, N),
        KOKKOS_LAMBDA(const int& i) { a_char_copy(i) = char(0); });
    Kokkos::parallel_for(
        "TestDeepCopy: FillB", policyB_t(0, N),
        KOKKOS_LAMBDA(const int& i) { b_char(i) = char(0); });
  }

  static int compare_equal(
      Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char_copy,
      Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char) {
    const int N = a_char.extent(0);
    int errors;
    Kokkos::parallel_reduce(
        "TestDeepCopy: FillA_copy", policyA_t(0, N),
        KOKKOS_LAMBDA(const int& i, int& lsum) {
          if (a_char_copy(i) != a_char(i)) lsum++;
        },
        errors);
    return errors;
  }

  static void run_test(int num_bytes) {
    a_base_t a_base("test_space_to_space", (num_bytes + 128) / 8);
    a_base_t a_base_copy("test_space_to_space", (num_bytes + 128) / 8);
    Kokkos::View<double*, Kokkos::LayoutRight, MemorySpaceB> b_base(
        "test_space_to_space", (num_bytes + 128) / 8);

    Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char(
        (char*)a_base.data(), a_base.extent(0) * 8);
    Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceA> a_char_copy(
        (char*)a_base_copy.data(), a_base.extent(0) * 8);
    Kokkos::View<char*, Kokkos::LayoutRight, MemorySpaceB> b_char(
        (char*)b_base.data(), b_base.extent(0) * 8);

    Kokkos::parallel_for(
        "TestDeepCopy: FillA", policyA_t(0, a_char.extent(0)),
        KOKKOS_LAMBDA(const int& i) {
          a_char(i) = static_cast<char>(i % 97) + 1;
        });

    reset_a_copy_and_b(a_char_copy, b_char);

    {
      int check = compare_equal(a_char_copy, a_char);
      ASSERT_EQ(check, a_char.extent(0));
    }

    // (a.data()%8, (a.data()+a.extent(0))%8, b.data()%8,
    // (b.data()+b.extent(0))%8 (0,0,0,0)
    {
      int a_begin = 0;
      int a_end   = 0;
      int b_begin = 0;
      int b_end   = 0;
      auto a      = Kokkos::subview(
          a_char, std::pair<int, int>(a_begin, a_char.extent(0) - a_end));
      auto b = Kokkos::subview(
          b_char, std::pair<int, int>(b_begin, b_char.extent(0) - b_end));
      auto a_copy = Kokkos::subview(
          a_char_copy,
          std::pair<int, int>(a_begin, a_char_copy.extent(0) - a_end));
      Kokkos::deep_copy(b, a);
      Kokkos::deep_copy(a_copy, b);
      int check = compare_equal(a_copy, a);
      ASSERT_EQ(check, 0);
    }

    {
      int a_begin = 0;
      int a_end   = 5;
      int b_begin = 0;
      int b_end   = 5;
      auto a      = Kokkos::subview(
          a_char, std::pair<int, int>(a_begin, a_char.extent(0) - a_end));
      auto b = Kokkos::subview(
          b_char, std::pair<int, int>(b_begin, b_char.extent(0) - b_end));
      auto a_copy = Kokkos::subview(
          a_char_copy,
          std::pair<int, int>(a_begin, a_char_copy.extent(0) - a_end));
      Kokkos::deep_copy(b, a);
      Kokkos::deep_copy(a_copy, b);
      int check = compare_equal(a_copy, a);
      ASSERT_EQ(check, 0);
    }

    {
      int a_begin = 3;
      int a_end   = 0;
      int b_begin = 3;
      int b_end   = 0;
      auto a      = Kokkos::subview(
          a_char, std::pair<int, int>(a_begin, a_char.extent(0) - a_end));
      auto b = Kokkos::subview(
          b_char, std::pair<int, int>(b_begin, b_char.extent(0) - b_end));
      auto a_copy = Kokkos::subview(
          a_char_copy,
          std::pair<int, int>(a_begin, a_char_copy.extent(0) - a_end));
      Kokkos::deep_copy(b, a);
      Kokkos::deep_copy(a_copy, b);
      int check = compare_equal(a_copy, a);
      ASSERT_EQ(check, 0);
    }

    {
      int a_begin = 3;
      int a_end   = 6;
      int b_begin = 3;
      int b_end   = 6;
      auto a      = Kokkos::subview(
          a_char, std::pair<int, int>(a_begin, a_char.extent(0) - a_end));
      auto b = Kokkos::subview(
          b_char, std::pair<int, int>(b_begin, b_char.extent(0) - b_end));
      auto a_copy = Kokkos::subview(
          a_char_copy,
          std::pair<int, int>(a_begin, a_char_copy.extent(0) - a_end));
      Kokkos::deep_copy(b, a);
      Kokkos::deep_copy(a_copy, b);
      int check = compare_equal(a_copy, a);
      ASSERT_EQ(check, 0);
    }

    {
      int a_begin = 5;
      int a_end   = 4;
      int b_begin = 3;
      int b_end   = 6;
      auto a      = Kokkos::subview(
          a_char, std::pair<int, int>(a_begin, a_char.extent(0) - a_end));
      auto b = Kokkos::subview(
          b_char, std::pair<int, int>(b_begin, b_char.extent(0) - b_end));
      auto a_copy = Kokkos::subview(
          a_char_copy,
          std::pair<int, int>(a_begin, a_char_copy.extent(0) - a_end));
      Kokkos::deep_copy(b, a);
      Kokkos::deep_copy(a_copy, b);
      int check = compare_equal(a_copy, a);
      ASSERT_EQ(check, 0);
    }

    {
      int a_begin = 0;
      int a_end   = 8;
      int b_begin = 2;
      int b_end   = 6;
      auto a      = Kokkos::subview(
          a_char, std::pair<int, int>(a_begin, a_char.extent(0) - a_end));
      auto b = Kokkos::subview(
          b_char, std::pair<int, int>(b_begin, b_char.extent(0) - b_end));
      auto a_copy = Kokkos::subview(
          a_char_copy,
          std::pair<int, int>(a_begin, a_char_copy.extent(0) - a_end));
      Kokkos::deep_copy(b, a);
      Kokkos::deep_copy(a_copy, b);
      int check = compare_equal(a_copy, a);
      ASSERT_EQ(check, 0);
    }

    {
      int a_begin = 2;
      int a_end   = 6;
      int b_begin = 0;
      int b_end   = 8;
      auto a      = Kokkos::subview(
          a_char, std::pair<int, int>(a_begin, a_char.extent(0) - a_end));
      auto b = Kokkos::subview(
          b_char, std::pair<int, int>(b_begin, b_char.extent(0) - b_end));
      auto a_copy = Kokkos::subview(
          a_char_copy,
          std::pair<int, int>(a_begin, a_char_copy.extent(0) - a_end));
      Kokkos::deep_copy(b, a);
      Kokkos::deep_copy(a_copy, b);
      int check = compare_equal(a_copy, a);
      ASSERT_EQ(check, 0);
    }
  }
};
}  // namespace Impl

TEST(TEST_CATEGORY, deep_copy_alignment) {
  {
    Impl::TestDeepCopy<TEST_EXECSPACE::memory_space,
                       TEST_EXECSPACE::memory_space>::run_test(100000);
  }
  {
    Impl::TestDeepCopy<Kokkos::HostSpace,
                       TEST_EXECSPACE::memory_space>::run_test(100000);
  }
  {
    Impl::TestDeepCopy<TEST_EXECSPACE::memory_space,
                       Kokkos::HostSpace>::run_test(100000);
  }
}
#endif

namespace Impl {
template <class Scalar1, class Scalar2, class Layout1, class Layout2>
struct TestDeepCopyScalarConversion {
  struct TagFill {};
  struct TagCompare {};

  using view_type_s1_1d = Kokkos::View<Scalar1*, Layout1, TEST_EXECSPACE>;
  using view_type_s2_1d = Kokkos::View<Scalar2*, Layout2, TEST_EXECSPACE>;
  using view_type_s1_2d = Kokkos::View<Scalar1**, Layout1, TEST_EXECSPACE>;
  using view_type_s2_2d = Kokkos::View<Scalar2**, Layout2, TEST_EXECSPACE>;

  using base_layout1 = typename std::conditional<
      std::is_same<Layout1, Kokkos::LayoutStride>::value, Kokkos::LayoutLeft,
      Layout1>::type;
  using base_layout2 = typename std::conditional<
      std::is_same<Layout2, Kokkos::LayoutStride>::value, Kokkos::LayoutLeft,
      Layout2>::type;

  using base_type_s1_1d = Kokkos::View<Scalar1*, base_layout1, TEST_EXECSPACE>;
  using base_type_s2_1d = Kokkos::View<Scalar2*, base_layout2, TEST_EXECSPACE>;
  using base_type_s1_2d = Kokkos::View<Scalar1**, base_layout1, TEST_EXECSPACE>;
  using base_type_s2_2d = Kokkos::View<Scalar2**, base_layout2, TEST_EXECSPACE>;

  view_type_s1_1d view_s1_1d;
  view_type_s2_1d view_s2_1d;
  view_type_s1_2d view_s1_2d;
  view_type_s2_2d view_s2_2d;

  Kokkos::View<int64_t, TEST_EXECSPACE> error_count;

  void create_views(int64_t N0, int64_t N1) {
    base_type_s1_1d b_s1_1d("TestDeepCopyConversion::b_s1_1d", N0);
    base_type_s2_1d b_s2_1d("TestDeepCopyConversion::b_s2_1d", N0);
    base_type_s1_2d b_s1_2d("TestDeepCopyConversion::b_s1_2d", N0, N1);
    base_type_s2_2d b_s2_2d("TestDeepCopyConversion::b_s2_2d", N0, N1);

    view_s1_1d = view_type_s1_1d(b_s1_1d, Kokkos::ALL);
    view_s2_1d = view_type_s2_1d(b_s2_1d, Kokkos::ALL);
    view_s1_2d = view_type_s1_2d(b_s1_2d, Kokkos::ALL, Kokkos::ALL);
    view_s2_2d = view_type_s2_2d(b_s2_2d, Kokkos::ALL, Kokkos::ALL);

    error_count = Kokkos::View<int64_t, TEST_EXECSPACE>(
        "TestDeepCopyConversion::error_count");
  }

  KOKKOS_FUNCTION
  void operator()(TagFill, const int64_t i) const {
    view_s2_1d(i) = static_cast<Scalar2>(i + 1);
    for (int64_t j = 0; j < static_cast<int64_t>(view_s2_2d.extent(1)); j++)
      view_s2_2d(i, j) = static_cast<Scalar2>((i + 1) * 1000 + j + 1);
  }

  KOKKOS_FUNCTION
  void operator()(TagCompare, const int64_t i) const {
    int64_t errors = 0;
    if (view_s1_1d(i) != static_cast<Scalar1>(static_cast<Scalar2>(i + 1)))
      errors++;
    for (int64_t j = 0; j < static_cast<int64_t>(view_s1_2d.extent(1)); j++) {
      if (view_s1_2d(i, j) !=
          static_cast<Scalar1>(static_cast<Scalar2>((i + 1) * 1000 + j + 1)))
        errors++;
    }
    if (errors > 0) Kokkos::atomic_add(&error_count(), errors);
  }

  void run_tests(int64_t N0, int64_t N1) {
    create_views(N0, N1);

    Kokkos::parallel_for("TestDeepCopyConversion::Fill",
                         Kokkos::RangePolicy<TEST_EXECSPACE, TagFill,
                                             Kokkos::IndexType<int64_t>>(0, N0),
                         *this);

    Kokkos::deep_copy(view_s1_1d, view_s2_1d);
    Kokkos::deep_copy(view_s1_2d, view_s2_2d);

    Kokkos::parallel_for("TestDeepCopyConversion::Compare",
                         Kokkos::RangePolicy<TEST_EXECSPACE, TagCompare,
                                             Kokkos::IndexType<int64_t>>(0, N0),
                         *this);

    int64_t errors = 0;
    Kokkos::deep_copy(errors, error_count);
    ASSERT_TRUE(errors == 0);

    Kokkos::deep_copy(view_s1_1d, static_cast<Scalar1>(0));
    Kokkos::deep_copy(view_s1_2d, static_cast<Scalar1>(0));

    Kokkos::parallel_for("TestDeepCopyConversion::Compare",
                         Kokkos::RangePolicy<TEST_EXECSPACE, TagCompare,
                                             Kokkos::IndexType<int64_t>>(0, N0),
                         *this);
    Kokkos::deep_copy(errors, error_count);
    ASSERT_TRUE(errors > 0);

    Kokkos::deep_copy(error_count, 0);
    Kokkos::deep_copy(TEST_EXECSPACE(), view_s1_1d, view_s2_1d);
    Kokkos::deep_copy(TEST_EXECSPACE(), view_s1_2d, view_s2_2d);

    Kokkos::parallel_for("TestDeepCopyConversion::Compare",
                         Kokkos::RangePolicy<TEST_EXECSPACE, TagCompare,
                                             Kokkos::IndexType<int64_t>>(0, N0),
                         *this);

    Kokkos::deep_copy(errors, error_count);
    ASSERT_TRUE(errors == 0);
  }
};
}  // namespace Impl

TEST(TEST_CATEGORY, deep_copy_conversion) {
  int64_t N0 = 19381;
  int64_t N1 = 17;

  using right  = Kokkos::LayoutRight;
  using left   = Kokkos::LayoutLeft;
  using stride = Kokkos::LayoutStride;

  Impl::TestDeepCopyScalarConversion<double, double, right, right>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, double, right, left>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, double, left, right>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, double, stride, right>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, double, right, stride>().run_tests(
      N0, N1);

  Impl::TestDeepCopyScalarConversion<double, float, right, right>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, float, right, left>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, float, left, right>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, float, stride, right>().run_tests(
      N0, N1);
  Impl::TestDeepCopyScalarConversion<double, float, right, stride>().run_tests(
      N0, N1);
}
}  // namespace Test
