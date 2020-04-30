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

namespace TestAtomicViews {

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
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& update,
                   const volatile value_type& input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      left_view;

  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      right_view;

  typedef Kokkos::View<DataType, Kokkos::LayoutStride, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      stride_view;

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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      if (left(i0) != left(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (right(i0) != right(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
#else
      if (left(i0) != left.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (right(i0) != right.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
#endif
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
  typedef DeviceType device;

  enum { N0 = 1000, N1 = 3, N2 = 5, N3 = 7 };

  typedef Kokkos::View<T, device> dView0;
  typedef Kokkos::View<T*, device> dView1;
  typedef Kokkos::View<T * [N1], device> dView2;
  typedef Kokkos::View<T * [N1][N2], device> dView3;
  typedef Kokkos::View<T * [N1][N2][N3], device> dView4;
  typedef Kokkos::View<const T * [N1][N2][N3], device> const_dView4;
  typedef Kokkos::View<T****, device, Kokkos::MemoryUnmanaged> dView4_unmanaged;
  typedef typename dView0::host_mirror_space host;

  typedef Kokkos::View<T, device, Kokkos::MemoryTraits<Kokkos::Atomic> > aView0;
  typedef Kokkos::View<T*, device, Kokkos::MemoryTraits<Kokkos::Atomic> >
      aView1;
  typedef Kokkos::View<T * [N1], device, Kokkos::MemoryTraits<Kokkos::Atomic> >
      aView2;
  typedef Kokkos::View<T * [N1][N2], device,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      aView3;
  typedef Kokkos::View<T * [N1][N2][N3], device,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      aView4;
  typedef Kokkos::View<const T * [N1][N2][N3], device,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      const_aView4;

  typedef Kokkos::View<
      T****, device, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic> >
      aView4_unmanaged;

  typedef typename aView0::host_mirror_space host_atomic;

  TestAtomicViewAPI() {
    TestViewOperator_LeftAndRight<int[2], device>::testit();
    run_test_rank0();
    run_test_rank4();
    run_test_const();
  }

  static void run_test_rank0() {
    dView0 dx, dy;
    aView0 ax, ay, az;

    dx = dView0("dx");
    dy = dView0("dy");
    ASSERT_EQ(dx.use_count(), size_t(1));
    ASSERT_EQ(dy.use_count(), size_t(1));

    ax = dx;
    ay = dy;
    ASSERT_EQ(dx.use_count(), size_t(2));
    ASSERT_EQ(dy.use_count(), size_t(2));
    ASSERT_EQ(dx.use_count(), ax.use_count());

    az = ax;
    ASSERT_EQ(dx.use_count(), size_t(3));
    ASSERT_EQ(ax.use_count(), size_t(3));
    ASSERT_EQ(az.use_count(), size_t(3));
    ASSERT_EQ(az.use_count(), ax.use_count());
  }

  static void run_test_rank4() {
    dView4 dx, dy;
    aView4 ax, ay, az;

    dx = dView4("dx", N0);
    dy = dView4("dy", N0);
    ASSERT_EQ(dx.use_count(), size_t(1));
    ASSERT_EQ(dy.use_count(), size_t(1));

    ax = dx;
    ay = dy;
    ASSERT_EQ(dx.use_count(), size_t(2));
    ASSERT_EQ(dy.use_count(), size_t(2));
    ASSERT_EQ(dx.use_count(), ax.use_count());

    dView4_unmanaged unmanaged_dx = dx;
    ASSERT_EQ(dx.use_count(), size_t(2));

    az = ax;
    ASSERT_EQ(dx.use_count(), size_t(3));
    ASSERT_EQ(ax.use_count(), size_t(3));
    ASSERT_EQ(az.use_count(), size_t(3));
    ASSERT_EQ(az.use_count(), ax.use_count());

    aView4_unmanaged unmanaged_ax = ax;
    ASSERT_EQ(ax.use_count(), size_t(3));

    aView4_unmanaged unmanaged_ax_from_ptr_dx = aView4_unmanaged(
        dx.data(), dx.extent(0), dx.extent(1), dx.extent(2), dx.extent(3));
    ASSERT_EQ(ax.use_count(), size_t(3));

    const_aView4 const_ax = ax;
    ASSERT_EQ(ax.use_count(), size_t(4));
    ASSERT_EQ(const_ax.use_count(), ax.use_count());

    ASSERT_FALSE(ax.data() == nullptr);
    ASSERT_FALSE(const_ax.data() == nullptr);  // referenceable ptr
    ASSERT_FALSE(unmanaged_ax.data() == nullptr);
    ASSERT_FALSE(unmanaged_ax_from_ptr_dx.data() == nullptr);
    ASSERT_FALSE(ay.data() == nullptr);
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

  typedef T DataType[2];

  static void check_auto_conversion_to_const(
      const Kokkos::View<const DataType, device,
                         Kokkos::MemoryTraits<Kokkos::Atomic> >& arg_const,
      const Kokkos::View<const DataType, device,
                         Kokkos::MemoryTraits<Kokkos::Atomic> >& arg) {
    ASSERT_TRUE(arg_const == arg);
  }

  static void run_test_const() {
    typedef Kokkos::View<DataType, device,
                         Kokkos::MemoryTraits<Kokkos::Atomic> >
        typeX;
    typedef Kokkos::View<const DataType, device,
                         Kokkos::MemoryTraits<Kokkos::Atomic> >
        const_typeX;

    typeX x("X");
    const_typeX xc = x;

    // ASSERT_TRUE( xc == x ); // const xc is referenceable, non-const x is not
    // ASSERT_TRUE( x == xc );

    check_auto_conversion_to_const(x, xc);
  }
};

//---------------------------------------------------
//-----------initialization functors-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct InitFunctor_Seq {
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef typename view_type::HostMirror host_view_type;

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
    const int64_t half_sum_end = (N / 2) - 1;
    const int64_t full_sum_end = N - 1;
    result[0] = half_sum_end * (half_sum_end + 1) / 2;  // Even sum.
    result[1] =
        (full_sum_end * (full_sum_end + 1) / 2) - result[0];  // Odd sum.
  } else {
    const int64_t half_sum_end = (T)(N / 2);
    const int64_t full_sum_end = N - 2;
    result[0] = half_sum_end * (half_sum_end - 1) / 2;  // Even sum.
    result[1] =
        (full_sum_end * (full_sum_end - 1) / 2) - result[0];  // Odd sum.
  }

  return (T)(result[0] + result[1]);
}

template <class T, class DeviceType>
bool PlusEqualAtomicViewTest(int64_t input_length) {
  T res       = PlusEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = PlusEqualAtomicViewCheck<T>(input_length);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = PlusEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//-----------atomic view minus-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct MinusEqualAtomicViewFunctor {
  typedef Kokkos::View<T*, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef typename view_type::HostMirror host_view_type;

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
    const int64_t half_sum_end = (N / 2) - 1;
    const int64_t full_sum_end = N - 1;
    result[0] = -1 * (half_sum_end * (half_sum_end + 1) / 2);  // Even sum.
    result[1] =
        -1 * ((full_sum_end * (full_sum_end + 1) / 2) + result[0]);  // Odd sum.
  } else {
    const int64_t half_sum_end = (int64_t)(N / 2);
    const int64_t full_sum_end = N - 2;
    result[0] = -1 * (half_sum_end * (half_sum_end - 1) / 2);  // Even sum.
    result[1] =
        -1 * ((full_sum_end * (full_sum_end - 1) / 2) + result[0]);  // Odd sum.
  }

  return (result[0] + result[1]);
}

template <class T, class DeviceType>
bool MinusEqualAtomicViewTest(int64_t input_length) {
  T res       = MinusEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = MinusEqualAtomicViewCheck<T>(input_length);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = MinusEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//-----------atomic view times-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct TimesEqualAtomicViewFunctor {
  typedef Kokkos::View<T*, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef typename view_type::HostMirror host_view_type;

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
bool TimesEqualAtomicViewTest(const int64_t input_length) {
  const int64_t remainder = 23;
  T res       = TimesEqualAtomicView<T, DeviceType>(input_length, remainder);
  T resSerial = TimesEqualAtomicViewCheck<T>(input_length, remainder);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = TimesEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//------------atomic view div-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct DivEqualAtomicViewFunctor {
  typedef Kokkos::View<T, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T, execution_space> scalar_view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T, execution_space> scalar_view_type;
  typedef typename scalar_view_type::HostMirror host_scalar_view_type;

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
bool DivEqualAtomicViewTest(const int64_t input_length) {
  const int64_t remainder = 23;

  T res       = DivEqualAtomicView<T, DeviceType>(input_length, remainder);
  T resSerial = DivEqualAtomicViewCheck<T>(input_length, remainder);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = DivEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//------------atomic view mod-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct ModEqualAtomicViewFunctor {
  typedef Kokkos::View<T, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T, execution_space> scalar_view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T, execution_space> scalar_view_type;
  typedef typename scalar_view_type::HostMirror host_scalar_view_type;

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
bool ModEqualAtomicViewTest(const int64_t input_length) {
  static_assert(std::is_integral<T>::value,
                "ModEqualAtomicView Error: Type must be integral type for this "
                "unit test");

  const int64_t remainder = 23;

  T res       = ModEqualAtomicView<T, DeviceType>(input_length, remainder);
  T resSerial = ModEqualAtomicViewCheck<T>(input_length, remainder);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = ModEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//------------atomic view rs-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct RSEqualAtomicViewFunctor {
  typedef Kokkos::View<T****, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T****, execution_space> result_view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T****, execution_space> result_view_type;
  typedef typename result_view_type::HostMirror host_scalar_view_type;

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
bool RSEqualAtomicViewTest(const int64_t input_length) {
  static_assert(std::is_integral<T>::value,
                "RSEqualAtomicViewTest: Must be integral type for test");

  const int64_t remainder = 61042;       // prime - 1
  const int64_t value     = 1073741825;  //  2^30+1
  T res = RSEqualAtomicView<T, DeviceType>(input_length, value, remainder);
  T resSerial = RSEqualAtomicViewCheck<T>(input_length, value, remainder);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = RSEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//------------atomic view ls-equal------------------
//---------------------------------------------------

template <class T, class execution_space>
struct LSEqualAtomicViewFunctor {
  typedef Kokkos::View<T****, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T****, execution_space> result_view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef Kokkos::View<T****, execution_space> result_view_type;
  typedef typename result_view_type::HostMirror host_scalar_view_type;

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
bool LSEqualAtomicViewTest(const int64_t input_length) {
  static_assert(std::is_integral<T>::value,
                "LSEqualAtomicViewTest: Must be integral type for test");

  const int64_t remainder = 61042;  // prime - 1
  const int64_t value     = 1;      //  2^30+1
  T res = LSEqualAtomicView<T, DeviceType>(input_length, value, remainder);
  T resSerial = LSEqualAtomicViewCheck<T>(input_length, value, remainder);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = RSEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//-----------atomic view and-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct AndEqualAtomicViewFunctor {
  typedef Kokkos::View<T*, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef typename view_type::HostMirror host_view_type;

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
    if (N % 2 == 0) {
      result[0] &= (T)i;
    } else {
      result[1] &= (T)i;
    }
  }

  return (result[0]);
}

template <class T, class DeviceType>
bool AndEqualAtomicViewTest(int64_t input_length) {
  static_assert(std::is_integral<T>::value,
                "AndEqualAtomicViewTest: Must be integral type for test");

  T res       = AndEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = AndEqualAtomicViewCheck<T>(input_length);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = AndEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//-----------atomic view or-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct OrEqualAtomicViewFunctor {
  typedef Kokkos::View<T*, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef typename view_type::HostMirror host_view_type;

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
bool OrEqualAtomicViewTest(int64_t input_length) {
  static_assert(std::is_integral<T>::value,
                "OrEqualAtomicViewTest: Must be integral type for test");

  T res       = OrEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = OrEqualAtomicViewCheck<T>(input_length);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = OrEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

//---------------------------------------------------
//-----------atomic view xor-equal-----------------
//---------------------------------------------------

template <class T, class execution_space>
struct XOrEqualAtomicViewFunctor {
  typedef Kokkos::View<T*, execution_space,
                       Kokkos::MemoryTraits<Kokkos::Atomic> >
      atomic_view_type;
  typedef Kokkos::View<T*, execution_space> view_type;

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
  typedef Kokkos::View<T*, execution_space> view_type;
  typedef typename view_type::HostMirror host_view_type;

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
bool XOrEqualAtomicViewTest(int64_t input_length) {
  static_assert(std::is_integral<T>::value,
                "XOrEqualAtomicViewTest: Must be integral type for test");

  T res       = XOrEqualAtomicView<T, DeviceType>(input_length);
  T resSerial = XOrEqualAtomicViewCheck<T>(input_length);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name()
              << ">( test = XOrEqualAtomicViewTest"
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

// inc/dec?

//---------------------------------------------------
//--------------atomic_test_control------------------
//---------------------------------------------------

template <class T, class DeviceType>
bool AtomicViewsTestIntegralType(const int length, int test) {
  static_assert(std::is_integral<T>::value,
                "TestAtomicViews Error: Non-integral type passed into "
                "IntegralType tests");

  switch (test) {
    case 1: return PlusEqualAtomicViewTest<T, DeviceType>(length);
    case 2: return MinusEqualAtomicViewTest<T, DeviceType>(length);
    case 3: return RSEqualAtomicViewTest<T, DeviceType>(length);
    case 4: return LSEqualAtomicViewTest<T, DeviceType>(length);
    case 5: return ModEqualAtomicViewTest<T, DeviceType>(length);
    case 6: return AndEqualAtomicViewTest<T, DeviceType>(length);
    case 7: return OrEqualAtomicViewTest<T, DeviceType>(length);
    case 8: return XOrEqualAtomicViewTest<T, DeviceType>(length);
  }

  return 0;
}

template <class T, class DeviceType>
bool AtomicViewsTestNonIntegralType(const int length, int test) {
  switch (test) {
    case 1: return PlusEqualAtomicViewTest<T, DeviceType>(length);
    case 2: return MinusEqualAtomicViewTest<T, DeviceType>(length);
    case 3: return TimesEqualAtomicViewTest<T, DeviceType>(length);
    case 4: return DivEqualAtomicViewTest<T, DeviceType>(length);
  }

  return 0;
}

}  // namespace TestAtomicViews

namespace Test {

TEST(TEST_CATEGORY, atomic_views_integral) {
  const int64_t length = 1000000;
  {
    // Integral Types.
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 1)));
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 2)));
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 3)));
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 4)));
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 5)));
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 6)));
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 7)));
    ASSERT_TRUE(
        (TestAtomicViews::AtomicViewsTestIntegralType<int64_t, TEST_EXECSPACE>(
            length, 8)));
  }
}

TEST(TEST_CATEGORY, atomic_views_nonintegral) {
  const int64_t length = 1000000;
  {
    // Non-Integral Types.
    ASSERT_TRUE((
        TestAtomicViews::AtomicViewsTestNonIntegralType<double, TEST_EXECSPACE>(
            length, 1)));
    ASSERT_TRUE((
        TestAtomicViews::AtomicViewsTestNonIntegralType<double, TEST_EXECSPACE>(
            length, 2)));
    ASSERT_TRUE((
        TestAtomicViews::AtomicViewsTestNonIntegralType<double, TEST_EXECSPACE>(
            length, 3)));
    ASSERT_TRUE((
        TestAtomicViews::AtomicViewsTestNonIntegralType<double, TEST_EXECSPACE>(
            length, 4)));
  }
}

TEST(TEST_CATEGORY, atomic_view_api) {
  TestAtomicViews::TestAtomicViewAPI<int, TEST_EXECSPACE>();
}
}  // namespace Test
