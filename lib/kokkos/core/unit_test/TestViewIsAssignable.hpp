#include <Kokkos_Core.hpp>

namespace Test {
namespace Impl {
template <class ViewTypeDst, class ViewTypeSrc>
struct TestAssignability {
  using mapping_type =
      Kokkos::Impl::ViewMapping<typename ViewTypeDst::traits,
                                typename ViewTypeSrc::traits,
                                typename ViewTypeDst::specialize>;

  template <class MappingType>
  static void try_assign(
      ViewTypeDst& dst, ViewTypeSrc& src,
      typename std::enable_if<MappingType::is_assignable>::type* = nullptr) {
    dst = src;
  }

  template <class MappingType>
  static void try_assign(
      ViewTypeDst&, ViewTypeSrc&,
      typename std::enable_if<!MappingType::is_assignable>::type* = nullptr) {
    Kokkos::Impl::throw_runtime_exception(
        "TestAssignability::try_assign: Unexpected call path");
  }

  template <class... Dimensions>
  static void test(bool always, bool sometimes, Dimensions... dims) {
    ViewTypeDst dst;
    ViewTypeSrc src("SRC", dims...);

    bool is_always_assignable =
        Kokkos::is_always_assignable<ViewTypeDst, ViewTypeSrc>::value;
    bool is_assignable = Kokkos::is_assignable(dst, src);

    // Print out if there is an error with typeid so you can just filter the
    // output with c++filt -t to see which assignment causes the error.
    if (is_always_assignable != always || is_assignable != sometimes)
      printf(
          "is_always_assignable: %i (%i), is_assignable: %i (%i) [ %s ] to [ "
          "%s ]\n",
          is_always_assignable ? 1 : 0, always ? 1 : 0, is_assignable ? 1 : 0,
          sometimes ? 1 : 0, typeid(ViewTypeSrc).name(),
          typeid(ViewTypeDst).name());
    if (sometimes) {
      ASSERT_NO_THROW(try_assign<mapping_type>(dst, src));
    }
    ASSERT_EQ(always, is_always_assignable);
    ASSERT_EQ(sometimes, is_assignable);
  }
};

}  // namespace Impl

TEST(TEST_CATEGORY, view_is_assignable) {
  using namespace Kokkos;
  using h_exec = typename DefaultHostExecutionSpace::memory_space;
  using d_exec = typename TEST_EXECSPACE::memory_space;
  using left   = LayoutLeft;
  using right  = LayoutRight;
  using stride = LayoutStride;
  // Static/Dynamic Extents
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int*, left, d_exec>>::test(true, true, 10);
  Impl::TestAssignability<View<int[10], left, d_exec>,
                          View<int*, left, d_exec>>::test(false, true, 10);
  Impl::TestAssignability<View<int[5], left, d_exec>,
                          View<int*, left, d_exec>>::test(false, false, 10);
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int[10], left, d_exec>>::test(true, true);
  Impl::TestAssignability<View<int[10], left, d_exec>,
                          View<int[10], left, d_exec>>::test(true, true);
  Impl::TestAssignability<View<int[5], left, d_exec>,
                          View<int[10], left, d_exec>>::test(false, false);
  Impl::TestAssignability<View<int**, left, d_exec>,
                          View<int**, left, d_exec>>::test(true, true, 10, 10);
  Impl::TestAssignability<View<int * [10], left, d_exec>,
                          View<int**, left, d_exec>>::test(false, true, 10, 10);
  Impl::TestAssignability<View<int * [5], left, d_exec>,
                          View<int**, left, d_exec>>::test(false, false, 10,
                                                           10);
  Impl::TestAssignability<View<int**, left, d_exec>,
                          View<int * [10], left, d_exec>>::test(true, true, 10);
  Impl::TestAssignability<View<int * [10], left, d_exec>,
                          View<int * [10], left, d_exec>>::test(true, true, 10);
  Impl::TestAssignability<View<int * [5], left, d_exec>,
                          View<int * [10], left, d_exec>>::test(false, false,
                                                                10);

  // Mismatch value_type
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<double*, left, d_exec>>::test(false, false, 10);

  // Layout assignment
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int*, right, d_exec>>::test(true, true, 10);

  // This could be made possible (due to the degenerate nature of the views) but
  // we do not allow this yet
  // TestAssignability<View<int**,left,d_exec>,View<int**,right,d_exec>>::test(false,true,10,1);
  Impl::TestAssignability<View<int**, left, d_exec>,
                          View<int**, right, d_exec>>::test(false, false, 10,
                                                            2);
  Impl::TestAssignability<View<int**, stride, d_exec>,
                          View<int**, right, d_exec>>::test(true, true, 10, 2);
  Impl::TestAssignability<View<int**, stride, d_exec>,
                          View<int**, left, d_exec>>::test(true, true, 10, 2);

  // Space Assignment
  bool expected = Kokkos::Impl::MemorySpaceAccess<d_exec, h_exec>::assignable;
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int*, left, h_exec>>::test(expected, expected,
                                                          10);
  expected = Kokkos::Impl::MemorySpaceAccess<h_exec, d_exec>::assignable;
  Impl::TestAssignability<View<int*, left, h_exec>,
                          View<int*, left, d_exec>>::test(expected, expected,
                                                          10);

  // reference type and const-qualified types
  using SomeViewType = View<int*, left, d_exec>;
#if defined(KOKKOS_ENABLE_CXX17)
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType>);
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType&>);
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType const>);
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType const&>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType&>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType const>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType const&>);
#else
  static_assert(is_always_assignable<SomeViewType, SomeViewType>::value, "");
  static_assert(is_always_assignable<SomeViewType, SomeViewType&>::value, "");
  static_assert(is_always_assignable<SomeViewType, SomeViewType const>::value,
                "");
  static_assert(is_always_assignable<SomeViewType, SomeViewType const&>::value,
                "");
  static_assert(is_always_assignable<SomeViewType&, SomeViewType>::value, "");
  static_assert(is_always_assignable<SomeViewType&, SomeViewType&>::value, "");
  static_assert(is_always_assignable<SomeViewType&, SomeViewType const>::value,
                "");
  static_assert(is_always_assignable<SomeViewType&, SomeViewType const&>::value,
                "");
#endif
}
}  // namespace Test
