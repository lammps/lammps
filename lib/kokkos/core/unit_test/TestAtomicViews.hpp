// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_TypeInfo.hpp>

namespace {

//-------------------------------------------------
//-----------atomic view api tests-----------------
//-------------------------------------------------

template <class T, class... P>
size_t allocation_count(const Kokkos::View<T, P...>& view) {
  const size_t card  = view.size();
  const size_t alloc = view.span();

  const int memory_span = Kokkos::View<int*>::required_allocation_size(100);

  return (card <= alloc && memory_span == 400) ? alloc : 0;
}

template <class DataType, class DeviceType,
          unsigned Rank = Kokkos::ViewTraits<DataType>::rank>
struct TestViewOperator_LeftAndRight;

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 1> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type& update, const value_type& input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space,
                                 Kokkos::MemoryTraits<Kokkos::Atomic> >;

  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space,
                   Kokkos::MemoryTraits<Kokkos::Atomic> >;

  using stride_view =
      Kokkos::View<DataType, Kokkos::LayoutStride, execution_space,
                   Kokkos::MemoryTraits<Kokkos::Atomic> >;

  left_view left;
  right_view right;
  stride_view left_stride;
  stride_view right_stride;
  int64_t left_alloc;
  int64_t right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_stride(left),
        right_stride(right),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  static void testit() {
    TestViewOperator_LeftAndRight driver;

    int error_flag = 0;

    Kokkos::parallel_reduce(1, driver, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type& update) const {
    for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
      // Below checks that values match, but unable to check the references.
      // Should this be able to be checked?
      if (left(i0) != left.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (right(i0) != right.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (left(i0) != left_stride(i0)) {
        update |= 4;
      }
      if (right(i0) != right_stride(i0)) {
        update |= 8;
      }
      /*
            if ( &left( i0 )  != &left( i0, 0, 0, 0, 0, 0, 0, 0 ) )  { update |=
         3; } if ( &right( i0 ) != &right( i0, 0, 0, 0, 0, 0, 0, 0 ) ) { update
         |= 3; } if ( &left( i0 )  != &left_stride( i0 ) ) { update |= 4; } if (
         &right( i0 ) != &right_stride( i0 ) ) { update |= 8; }
      */
    }
  }
};

template <typename T, class DeviceType>
class TestAtomicViewAPI {
 public:
  using device = DeviceType;

  enum { N0 = 1000, N1 = 3, N2 = 5, N3 = 7 };

  using dView0           = Kokkos::View<T, device>;
  using dView1           = Kokkos::View<T*, device>;
  using dView2           = Kokkos::View<T* [N1], device>;
  using dView3           = Kokkos::View<T* [N1][N2], device>;
  using dView4           = Kokkos::View<T* [N1][N2][N3], device>;
  using const_dView4     = Kokkos::View<const T* [N1][N2][N3], device>;
  using dView4_unmanaged = Kokkos::View<T****, device, Kokkos::MemoryUnmanaged>;
  using host             = typename dView0::host_mirror_space;

  using aView0 = Kokkos::View<T, device, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using aView1 =
      Kokkos::View<T*, device, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using aView2 =
      Kokkos::View<T* [N1], device, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using aView3 =
      Kokkos::View<T* [N1][N2], device, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using aView4       = Kokkos::View<T* [N1][N2][N3], device,
                              Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using const_aView4 = Kokkos::View<const T* [N1][N2][N3], device,
                                    Kokkos::MemoryTraits<Kokkos::Atomic> >;

  using aView4_unmanaged =
      Kokkos::View<T****, device,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic> >;

  using host_atomic = typename aView0::host_mirror_space;

  TestAtomicViewAPI() {
    // FIXME_OPENMPTARGET
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    TestViewOperator_LeftAndRight<int[2], device>::testit();
#endif
    run_test_rank0();
    run_test_rank4();
    run_test_const();
  }

  static void run_test_rank0() {
    dView0 dx, dy;
    aView0 ax, ay, az;

    dx = dView0("dx");
    dy = dView0("dy");
    ASSERT_EQ(dx.use_count(), 1);
    ASSERT_EQ(dy.use_count(), 1);

    ax = dx;
    ay = dy;
    ASSERT_EQ(dx.use_count(), 2);
    ASSERT_EQ(dy.use_count(), 2);
    ASSERT_EQ(dx.use_count(), ax.use_count());

    az = ax;
    ASSERT_EQ(dx.use_count(), 3);
    ASSERT_EQ(ax.use_count(), 3);
    ASSERT_EQ(az.use_count(), 3);
    ASSERT_EQ(az.use_count(), ax.use_count());
  }

  static void run_test_rank4() {
    dView4 dx, dy;
    aView4 ax, ay, az;

    dx = dView4("dx", N0);
    dy = dView4("dy", N0);
    ASSERT_EQ(dx.use_count(), 1);
    ASSERT_EQ(dy.use_count(), 1);

    ax = dx;
    ay = dy;
    ASSERT_EQ(dx.use_count(), 2);
    ASSERT_EQ(dy.use_count(), 2);
    ASSERT_EQ(dx.use_count(), ax.use_count());

    dView4_unmanaged unmanaged_dx = dx;
    ASSERT_EQ(dx.use_count(), 2);
    // Legacy unmanaged view can still track the use count when being
    // constructed from a managed view New view behavior guarantees returning 0
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
    ASSERT_EQ(unmanaged_dx.use_count(), 2);
#else
    ASSERT_EQ(unmanaged_dx.use_count(), 0);
#endif

    az = ax;
    ASSERT_EQ(dx.use_count(), 3);
    ASSERT_EQ(ax.use_count(), 3);
    ASSERT_EQ(az.use_count(), 3);
    ASSERT_EQ(az.use_count(), ax.use_count());

    aView4_unmanaged unmanaged_ax = ax;
    ASSERT_EQ(ax.use_count(), 3);

    aView4_unmanaged unmanaged_ax_from_ptr_dx = aView4_unmanaged(
        dx.data(), dx.extent(0), dx.extent(1), dx.extent(2), dx.extent(3));
    ASSERT_EQ(ax.use_count(), 3);

    const_aView4 const_ax = ax;
    ASSERT_EQ(ax.use_count(), 4);
    ASSERT_EQ(const_ax.use_count(), ax.use_count());

    ASSERT_NE(ax.data(), nullptr);
    ASSERT_NE(const_ax.data(), nullptr);  // referenceable ptr
    ASSERT_NE(unmanaged_ax.data(), nullptr);
    ASSERT_NE(unmanaged_ax_from_ptr_dx.data(), nullptr);
    ASSERT_NE(ay.data(), nullptr);
    //    ASSERT_NE( ax, ay );
    //    Above test results in following runtime error from gtest:
    //    Expected: (ax) != (ay), actual: 32-byte object <30-01 D0-A0 D8-7F
    //    00-00 00-31 44-0C 01-00 00-00 E8-03 00-00 00-00 00-00 69-00 00-00
    //    00-00 00-00> vs 32-byte object <80-01 D0-A0 D8-7F 00-00 00-A1 4A-0C
    //    01-00 00-00 E8-03 00-00 00-00 00-00 69-00 00-00 00-00 00-00>

    ASSERT_EQ(ax.extent(0), unsigned(N0));
    ASSERT_EQ(ax.extent(1), unsigned(N1));
    ASSERT_EQ(ax.extent(2), unsigned(N2));
    ASSERT_EQ(ax.extent(3), unsigned(N3));

    ASSERT_EQ(ay.extent(0), unsigned(N0));
    ASSERT_EQ(ay.extent(1), unsigned(N1));
    ASSERT_EQ(ay.extent(2), unsigned(N2));
    ASSERT_EQ(ay.extent(3), unsigned(N3));

    ASSERT_EQ(unmanaged_ax_from_ptr_dx.span(),
              unsigned(N0) * unsigned(N1) * unsigned(N2) * unsigned(N3));
  }

  using DataType = T[2];

  static void check_auto_conversion_to_const(
      const Kokkos::View<const DataType, device,
                         Kokkos::MemoryTraits<Kokkos::Atomic> >& arg_const,
      const Kokkos::View<const DataType, device,
                         Kokkos::MemoryTraits<Kokkos::Atomic> >& arg) {
    ASSERT_EQ(arg_const, arg);
  }

  static void run_test_const() {
    using typeX =
        Kokkos::View<DataType, device, Kokkos::MemoryTraits<Kokkos::Atomic> >;
    using const_typeX = Kokkos::View<const DataType, device,
                                     Kokkos::MemoryTraits<Kokkos::Atomic> >;

    typeX x("X");
    const_typeX xc = x;

    // ASSERT_EQ( xc ,  x ); // const xc is referenceable, non-const x is not
    // ASSERT_EQ( x ,  xc );

    check_auto_conversion_to_const(x, xc);
  }
};

//---------------------------------------------------
//-----------initialization functors-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct InitFunctor_Seq {
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  const int64_t length;

  InitFunctor_Seq(view_type& input_, const int64_t length_)
      : input(input_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      input(i) = (T)i;
    }
  }
};

template <class T, class execution_space>
struct InitFunctor_ModTimes {
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  const int64_t length;
  const int64_t remainder;

  InitFunctor_ModTimes(view_type& input_, const int64_t length_,
                       const int64_t remainder_)
      : input(input_), length(length_), remainder(remainder_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % (remainder + 1) == remainder) {
        input(i) = (T)2;
      } else {
        input(i) = (T)1;
      }
    }
  }
};

template <class T, class execution_space>
struct InitFunctor_ModShift {
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  const int64_t length;
  const int64_t remainder;

  InitFunctor_ModShift(view_type& input_, const int64_t length_,
                       const int64_t remainder_)
      : input(input_), length(length_), remainder(remainder_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % (remainder + 1) == remainder) {
        input(i) = 1;
      }
    }
  }
};

//---------------------------------------------------
//-----------atomic view plus-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct PlusEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T*, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  atomic_view_type even_odd_result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator
  PlusEqualAtomicViewFunctor(const view_type& input_,
                             view_type& even_odd_result_, const int64_t length_)
      : input(input_), even_odd_result(even_odd_result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % 2 == 0) {
        even_odd_result(0) += input(i);
      } else {
        even_odd_result(1) += input(i);
      }
    }
  }
};

template <class T, class execution_space>
T PlusEqualAtomicView(const int64_t input_length) {
  using view_type      = Kokkos::View<T*, execution_space>;
  using host_view_type = typename view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  view_type result_view("result_view", 2);

  InitFunctor_Seq<T, execution_space> init_f(input, length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  PlusEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                         length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(0) + h_result_view(1));
}

template <class T>
T PlusEqualAtomicViewCheck(const int64_t input_length) {
  const int64_t N = input_length;
  T result[2];

  if (N % 2 == 0) {
    const T half_sum_end = static_cast<T>(N) / 2 - 1;
    const T full_sum_end = static_cast<T>(N) - 1;
    result[0]            = half_sum_end * (half_sum_end + 1) / 2;  // Even sum.
    result[1] =
        (full_sum_end * (full_sum_end + 1) / 2) - result[0];  // Odd sum.
  } else {
    const T half_sum_end = static_cast<T>(N) / 2;
    const T full_sum_end = static_cast<T>(N) - 2;
    result[0]            = half_sum_end * (half_sum_end - 1) / 2;  // Even sum.
    result[1] =
        (full_sum_end * (full_sum_end - 1) / 2) - result[0];  // Odd sum.
  }

  return (T)(result[0] + result[1]);
}

template <class T, class DeviceType>
void PlusEqualAtomicViewTest(int64_t input_length) {
  T res       = PlusEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = PlusEqualAtomicViewCheck<T>(input_length);

  ASSERT_EQ(res, resSerial)
      << "PlusEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ")";
}

//---------------------------------------------------
//-----------atomic view minus-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct MinusEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T*, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  atomic_view_type even_odd_result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator.
  MinusEqualAtomicViewFunctor(const view_type& input_,
                              view_type& even_odd_result_,
                              const int64_t length_)
      : input(input_), even_odd_result(even_odd_result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % 2 == 0) {
        even_odd_result(0) -= input(i);
      } else {
        even_odd_result(1) -= input(i);
      }
    }
  }
};

template <class T, class execution_space>
T MinusEqualAtomicView(const int64_t input_length) {
  using view_type      = Kokkos::View<T*, execution_space>;
  using host_view_type = typename view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  view_type result_view("result_view", 2);

  InitFunctor_Seq<T, execution_space> init_f(input, length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  MinusEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                          length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(0) + h_result_view(1));
}

template <class T>
T MinusEqualAtomicViewCheck(const int64_t input_length) {
  const int64_t N = input_length;
  T result[2];

  if (N % 2 == 0) {
    const T half_sum_end = static_cast<T>(N) / 2 - 1;
    const T full_sum_end = static_cast<T>(N) - 1;
    result[0] = -1 * (half_sum_end * (half_sum_end + 1) / 2);  // Even sum.
    result[1] =
        -1 * ((full_sum_end * (full_sum_end + 1) / 2) + result[0]);  // Odd sum.
  } else {
    const T half_sum_end = static_cast<T>(N) / 2;
    const T full_sum_end = static_cast<T>(N) - 2;
    result[0] = -1 * (half_sum_end * (half_sum_end - 1) / 2);  // Even sum.
    result[1] =
        -1 * ((full_sum_end * (full_sum_end - 1) / 2) + result[0]);  // Odd sum.
  }

  return (result[0] + result[1]);
}

template <class T, class DeviceType>
void MinusEqualAtomicViewTest(int64_t input_length) {
  T res       = MinusEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = MinusEqualAtomicViewCheck<T>(input_length);

  ASSERT_EQ(res, resSerial)
      << "MinusEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ")";
}

//---------------------------------------------------
//-----------atomic view times-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct TimesEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T*, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  atomic_view_type result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator
  TimesEqualAtomicViewFunctor(const view_type& input_, view_type& result_,
                              const int64_t length_)
      : input(input_), result(result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length && i > 0) {
      result(0) *= (double)input(i);
    }
  }
};

template <class T, class execution_space>
T TimesEqualAtomicView(const int64_t input_length, const int64_t remainder) {
  using view_type      = Kokkos::View<T*, execution_space>;
  using host_view_type = typename view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  view_type result_view("result_view", 1);
  deep_copy(result_view, 1.0);

  InitFunctor_ModTimes<T, execution_space> init_f(input, length, remainder);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  TimesEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                          length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(0));
}

template <class T>
T TimesEqualAtomicViewCheck(const int64_t input_length,
                            const int64_t remainder) {
  // Analytical result.
  const int64_t N = input_length;
  T result        = 1.0;

  for (int64_t i = 2; i < N; ++i) {
    if (i % (remainder + 1) == remainder) {
      result *= 2.0;
    } else {
      result *= 1.0;
    }
  }

  return (T)result;
}

template <class T, class DeviceType>
void TimesEqualAtomicViewTest(const int64_t input_length) {
  const int64_t remainder = 23;
  T res       = TimesEqualAtomicView<T, DeviceType>(input_length, remainder);
  T resSerial = TimesEqualAtomicViewCheck<T>(input_length, remainder);

  ASSERT_EQ(res, resSerial)
      << "TimesEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ")";
}

//---------------------------------------------------
//------------atomic view div-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct DivEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type        = Kokkos::View<T*, execution_space>;
  using scalar_view_type = Kokkos::View<T, execution_space>;

  view_type input;
  atomic_view_type result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator.
  DivEqualAtomicViewFunctor(const view_type& input_, scalar_view_type& result_,
                            const int64_t length_)
      : input(input_), result(result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length && i > 0) {
      result() /= (double)(input(i));
    }
  }
};

template <class T, class execution_space>
T DivEqualAtomicView(const int64_t input_length, const int64_t remainder) {
  using view_type             = Kokkos::View<T*, execution_space>;
  using scalar_view_type      = Kokkos::View<T, execution_space>;
  using host_scalar_view_type = typename scalar_view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  scalar_view_type result_view("result_view");
  Kokkos::deep_copy(result_view, 12121212121);

  InitFunctor_ModTimes<T, execution_space> init_f(input, length, remainder);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  DivEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                        length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_scalar_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view());
}

template <class T>
T DivEqualAtomicViewCheck(const int64_t input_length, const int64_t remainder) {
  const int64_t N = input_length;
  T result        = 12121212121.0;
  for (int64_t i = 2; i < N; ++i) {
    if (i % (remainder + 1) == remainder) {
      result /= 1.0;
    } else {
      result /= 2.0;
    }
  }

  return (T)result;
}

template <class T, class DeviceType>
void DivEqualAtomicViewTest(const int64_t input_length) {
  const int64_t remainder = 23;

  T res       = DivEqualAtomicView<T, DeviceType>(input_length, remainder);
  T resSerial = DivEqualAtomicViewCheck<T>(input_length, remainder);

  ASSERT_EQ(res, resSerial)
      << "DivEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ",remainder=" << remainder << ")";
}

//---------------------------------------------------
//------------atomic view mod-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct ModEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type        = Kokkos::View<T*, execution_space>;
  using scalar_view_type = Kokkos::View<T, execution_space>;

  view_type input;
  atomic_view_type result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator.
  ModEqualAtomicViewFunctor(const view_type& input_, scalar_view_type& result_,
                            const int64_t length_)
      : input(input_), result(result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length && i > 0) {
      result() %= (double)(input(i));
    }
  }
};

template <class T, class execution_space>
T ModEqualAtomicView(const int64_t input_length, const int64_t remainder) {
  using view_type             = Kokkos::View<T*, execution_space>;
  using scalar_view_type      = Kokkos::View<T, execution_space>;
  using host_scalar_view_type = typename scalar_view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  scalar_view_type result_view("result_view");
  Kokkos::deep_copy(result_view, 12121212121);

  InitFunctor_ModTimes<T, execution_space> init_f(input, length, remainder);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  ModEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                        length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_scalar_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view());
}

template <class T>
T ModEqualAtomicViewCheck(const int64_t input_length, const int64_t remainder) {
  const int64_t N = input_length;
  T result        = 12121212121;
  for (int64_t i = 2; i < N; ++i) {
    if (i % (remainder + 1) == remainder) {
      result %= 1;
    } else {
      result %= 2;
    }
  }

  return (T)result;
}

template <class T, class DeviceType>
void ModEqualAtomicViewTest(const int64_t input_length) {
  static_assert(std::is_integral_v<T>,
                "ModEqualAtomicView Error: Type must be integral type for this "
                "unit test");

  const int64_t remainder = 23;

  T res       = ModEqualAtomicView<T, DeviceType>(input_length, remainder);
  T resSerial = ModEqualAtomicViewCheck<T>(input_length, remainder);

  ASSERT_EQ(res, resSerial)
      << "ModEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ",remainder=" << remainder << ")";
}

//---------------------------------------------------
//------------atomic view rs-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct RSEqualAtomicViewFunctor {
  using atomic_view_type = Kokkos::View<T****, execution_space,
                                        Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type        = Kokkos::View<T*, execution_space>;
  using result_view_type = Kokkos::View<T****, execution_space>;

  const view_type input;
  atomic_view_type result;
  const int64_t length;
  const int64_t value;

  // Wrap the result view in an atomic view, use this for operator.
  RSEqualAtomicViewFunctor(const view_type& input_, result_view_type& result_,
                           const int64_t& length_, const int64_t& value_)
      : input(input_), result(result_), length(length_), value(value_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % 4 == 0) {
        result(1, 0, 0, 0) >>= input(i);
      } else if (i % 4 == 1) {
        result(0, 1, 0, 0) >>= input(i);
      } else if (i % 4 == 2) {
        result(0, 0, 1, 0) >>= input(i);
      } else if (i % 4 == 3) {
        result(0, 0, 0, 1) >>= input(i);
      }
    }
  }
};

template <class T, class execution_space>
T RSEqualAtomicView(const int64_t input_length, const int64_t value,
                    const int64_t remainder) {
  using view_type             = Kokkos::View<T*, execution_space>;
  using result_view_type      = Kokkos::View<T****, execution_space>;
  using host_scalar_view_type = typename result_view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  result_view_type result_view("result_view", 2, 2, 2, 2);
  host_scalar_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  h_result_view(1, 0, 0, 0)           = value;
  h_result_view(0, 1, 0, 0)           = value;
  h_result_view(0, 0, 1, 0)           = value;
  h_result_view(0, 0, 0, 1)           = value;
  Kokkos::deep_copy(result_view, h_result_view);

  InitFunctor_ModShift<T, execution_space> init_f(input, length, remainder);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  RSEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                       length, value);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(1, 0, 0, 0));
}

template <class T>
T RSEqualAtomicViewCheck(const int64_t input_length, const int64_t value,
                         const int64_t remainder) {
  T result[4];
  result[0] = value;
  result[1] = value;
  result[2] = value;
  result[3] = value;

  T* input = new T[input_length];
  for (int64_t i = 0; i < input_length; ++i) {
    if (i % (remainder + 1) == remainder) {
      input[i] = 1;
    } else {
      input[i] = 0;
    }
  }

  for (int64_t i = 0; i < input_length; ++i) {
    if (i % 4 == 0) {
      result[0] >>= input[i];
    } else if (i % 4 == 1) {
      result[1] >>= input[i];
    } else if (i % 4 == 2) {
      result[2] >>= input[i];
    } else if (i % 4 == 3) {
      result[3] >>= input[i];
    }
  }

  delete[] input;

  return (T)result[0];
}

template <class T, class DeviceType>
void RSEqualAtomicViewTest(const int64_t input_length) {
  static_assert(std::is_integral_v<T>,
                "RSEqualAtomicViewTest: Must be integral type for test");

  const int64_t remainder = 61042;       // prime - 1
  const int64_t value     = 1073741825;  //  2^30+1
  T res = RSEqualAtomicView<T, DeviceType>(input_length, value, remainder);
  T resSerial = RSEqualAtomicViewCheck<T>(input_length, value, remainder);

  ASSERT_EQ(res, resSerial)
      << "RSEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ",value=" << value
      << ",remainder=" << remainder << ")";
}

//---------------------------------------------------
//------------atomic view ls-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct LSEqualAtomicViewFunctor {
  using atomic_view_type = Kokkos::View<T****, execution_space,
                                        Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type        = Kokkos::View<T*, execution_space>;
  using result_view_type = Kokkos::View<T****, execution_space>;

  view_type input;
  atomic_view_type result;
  const int64_t length;
  const int64_t value;

  // Wrap the result view in an atomic view, use this for operator.
  LSEqualAtomicViewFunctor(const view_type& input_, result_view_type& result_,
                           const int64_t& length_, const int64_t& value_)
      : input(input_), result(result_), length(length_), value(value_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % 4 == 0) {
        result(1, 0, 0, 0) <<= input(i);
      } else if (i % 4 == 1) {
        result(0, 1, 0, 0) <<= input(i);
      } else if (i % 4 == 2) {
        result(0, 0, 1, 0) <<= input(i);
      } else if (i % 4 == 3) {
        result(0, 0, 0, 1) <<= input(i);
      }
    }
  }
};

template <class T, class execution_space>
T LSEqualAtomicView(const int64_t input_length, const int64_t value,
                    const int64_t remainder) {
  using view_type             = Kokkos::View<T*, execution_space>;
  using result_view_type      = Kokkos::View<T****, execution_space>;
  using host_scalar_view_type = typename result_view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  result_view_type result_view("result_view", 2, 2, 2, 2);
  host_scalar_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  h_result_view(1, 0, 0, 0)           = value;
  h_result_view(0, 1, 0, 0)           = value;
  h_result_view(0, 0, 1, 0)           = value;
  h_result_view(0, 0, 0, 1)           = value;
  Kokkos::deep_copy(result_view, h_result_view);

  InitFunctor_ModShift<T, execution_space> init_f(input, length, remainder);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  LSEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                       length, value);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(1, 0, 0, 0));
}

template <class T>
T LSEqualAtomicViewCheck(const int64_t input_length, const int64_t value,
                         const int64_t remainder) {
  T result[4];
  result[0] = value;
  result[1] = value;
  result[2] = value;
  result[3] = value;

  T* input = new T[input_length];
  for (int64_t i = 0; i < input_length; ++i) {
    if (i % (remainder + 1) == remainder) {
      input[i] = 1;
    } else {
      input[i] = 0;
    }
  }

  for (int64_t i = 0; i < input_length; ++i) {
    if (i % 4 == 0) {
      result[0] <<= input[i];
    } else if (i % 4 == 1) {
      result[1] <<= input[i];
    } else if (i % 4 == 2) {
      result[2] <<= input[i];
    } else if (i % 4 == 3) {
      result[3] <<= input[i];
    }
  }

  delete[] input;

  return (T)result[0];
}

template <class T, class DeviceType>
void LSEqualAtomicViewTest(const int64_t input_length) {
  static_assert(std::is_integral_v<T>,
                "LSEqualAtomicViewTest: Must be integral type for test");

  const int64_t remainder = 61042;  // prime - 1
  const int64_t value     = 1;      //  2^30+1
  T res = LSEqualAtomicView<T, DeviceType>(input_length, value, remainder);
  T resSerial = LSEqualAtomicViewCheck<T>(input_length, value, remainder);

  ASSERT_EQ(res, resSerial)
      << "LSEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ",value=" << value
      << ",remainder=" << remainder << ")";
}

//---------------------------------------------------
//-----------atomic view and-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct AndEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T*, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  atomic_view_type even_odd_result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator.
  AndEqualAtomicViewFunctor(const view_type& input_,
                            view_type& even_odd_result_, const int64_t length_)
      : input(input_), even_odd_result(even_odd_result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % 2 == 0) {
        even_odd_result(0) &= input(i);
      } else {
        even_odd_result(1) &= input(i);
      }
    }
  }
};

template <class T, class execution_space>
T AndEqualAtomicView(const int64_t input_length) {
  using view_type      = Kokkos::View<T*, execution_space>;
  using host_view_type = typename view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  view_type result_view("result_view", 2);
  Kokkos::deep_copy(result_view, 1);

  InitFunctor_Seq<T, execution_space> init_f(input, length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  AndEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                        length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(0));
}

template <class T>
T AndEqualAtomicViewCheck(const int64_t input_length) {
  const int64_t N = input_length;
  T result[2]     = {1};
  for (int64_t i = 0; i < N; ++i) {
    int64_t idx = N % 2;
    result[idx] &= (T)i;
  }
  return (result[0]);
}

template <class T, class DeviceType>
void AndEqualAtomicViewTest(int64_t input_length) {
  static_assert(std::is_integral_v<T>,
                "AndEqualAtomicViewTest: Must be integral type for test");

  T res       = AndEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = AndEqualAtomicViewCheck<T>(input_length);

  ASSERT_EQ(res, resSerial)
      << "AndEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ")";
}

//---------------------------------------------------
//-----------atomic view or-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct OrEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T*, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  atomic_view_type even_odd_result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator.
  OrEqualAtomicViewFunctor(const view_type& input_, view_type& even_odd_result_,
                           const int64_t length_)
      : input(input_), even_odd_result(even_odd_result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % 2 == 0) {
        even_odd_result(0) |= input(i);
      } else {
        even_odd_result(1) |= input(i);
      }
    }
  }
};

template <class T, class execution_space>
T OrEqualAtomicView(const int64_t input_length) {
  using view_type      = Kokkos::View<T*, execution_space>;
  using host_view_type = typename view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  view_type result_view("result_view", 2);

  InitFunctor_Seq<T, execution_space> init_f(input, length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  OrEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                       length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(0));
}

template <class T>
T OrEqualAtomicViewCheck(const int64_t input_length) {
  const int64_t N = input_length;
  T result[2]     = {0};
  for (int64_t i = 0; i < N; ++i) {
    if (i % 2 == 0) {
      result[0] |= (T)i;
    } else {
      result[1] |= (T)i;
    }
  }

  return (T)(result[0]);
}

template <class T, class DeviceType>
void OrEqualAtomicViewTest(int64_t input_length) {
  static_assert(std::is_integral_v<T>,
                "OrEqualAtomicViewTest: Must be integral type for test");

  T res       = OrEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = OrEqualAtomicViewCheck<T>(input_length);

  ASSERT_EQ(res, resSerial)
      << "OrEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ")";
}

//---------------------------------------------------
//-----------atomic view xor-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct XOrEqualAtomicViewFunctor {
  using atomic_view_type =
      Kokkos::View<T*, execution_space, Kokkos::MemoryTraits<Kokkos::Atomic> >;
  using view_type = Kokkos::View<T*, execution_space>;

  view_type input;
  atomic_view_type even_odd_result;
  const int64_t length;

  // Wrap the result view in an atomic view, use this for operator.
  XOrEqualAtomicViewFunctor(const view_type& input_,
                            view_type& even_odd_result_, const int64_t length_)
      : input(input_), even_odd_result(even_odd_result_), length(length_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int64_t i) const {
    if (i < length) {
      if (i % 2 == 0) {
        even_odd_result(0) ^= input(i);
      } else {
        even_odd_result(1) ^= input(i);
      }
    }
  }
};

template <class T, class execution_space>
T XOrEqualAtomicView(const int64_t input_length) {
  using view_type      = Kokkos::View<T*, execution_space>;
  using host_view_type = typename view_type::host_mirror_type;

  const int64_t length = input_length;

  view_type input("input_view", length);
  view_type result_view("result_view", 2);

  InitFunctor_Seq<T, execution_space> init_f(input, length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length), init_f);

  XOrEqualAtomicViewFunctor<T, execution_space> functor(input, result_view,
                                                        length);
  Kokkos::parallel_for(Kokkos::RangePolicy<execution_space>(0, length),
                       functor);
  Kokkos::fence();

  host_view_type h_result_view = Kokkos::create_mirror_view(result_view);
  Kokkos::deep_copy(h_result_view, result_view);

  return (T)(h_result_view(0));
}

template <class T>
T XOrEqualAtomicViewCheck(const int64_t input_length) {
  const int64_t N = input_length;
  T result[2]     = {0};
  for (int64_t i = 0; i < N; ++i) {
    if (i % 2 == 0) {
      result[0] ^= (T)i;
    } else {
      result[1] ^= (T)i;
    }
  }

  return (T)(result[0]);
}

template <class T, class DeviceType>
void XOrEqualAtomicViewTest(int64_t input_length) {
  static_assert(std::is_integral_v<T>,
                "XOrEqualAtomicViewTest: Must be integral type for test");

  T res       = XOrEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = XOrEqualAtomicViewCheck<T>(input_length);

  ASSERT_EQ(res, resSerial)
      << "XOrEqualAtomicViewTest<" << Kokkos::Impl::TypeInfo<T>::name()
      << ">(length=" << input_length << ")";
}

// inc/dec?

TEST(TEST_CATEGORY, atomic_views_integral) {
  const int64_t length = 1000000;
  PlusEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
  MinusEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
  RSEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
  LSEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
  ModEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
  AndEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
  OrEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
  XOrEqualAtomicViewTest<int64_t, TEST_EXECSPACE>(length);
}

TEST(TEST_CATEGORY, atomic_views_nonintegral) {
  const int64_t length = 1000000;
  PlusEqualAtomicViewTest<double, TEST_EXECSPACE>(length);
  MinusEqualAtomicViewTest<double, TEST_EXECSPACE>(length);
  TimesEqualAtomicViewTest<double, TEST_EXECSPACE>(length);
  DivEqualAtomicViewTest<double, TEST_EXECSPACE>(length);
}

TEST(TEST_CATEGORY, atomic_view_api) {
  TestAtomicViewAPI<int, TEST_EXECSPACE>();
}

}  // namespace
