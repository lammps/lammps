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

namespace TestAtomic {

// Struct for testing arbitrary size atomics.

template <int N>
struct SuperScalar {
  double val[N];

  KOKKOS_INLINE_FUNCTION
  SuperScalar() {
    for (int i = 0; i < N; i++) {
      val[i] = 0.0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar(const SuperScalar& src) {
    for (int i = 0; i < N; i++) {
      val[i] = src.val[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar(const volatile SuperScalar& src) {
    for (int i = 0; i < N; i++) {
      val[i] = src.val[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar& operator=(const SuperScalar& src) {
    for (int i = 0; i < N; i++) {
      val[i] = src.val[i];
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar& operator=(const volatile SuperScalar& src) {
    for (int i = 0; i < N; i++) {
      val[i] = src.val[i];
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator=(const SuperScalar& src) volatile {
    for (int i = 0; i < N; i++) {
      val[i] = src.val[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar operator+(const SuperScalar& src) const {
    SuperScalar tmp = *this;
    for (int i = 0; i < N; i++) {
      tmp.val[i] += src.val[i];
    }
    return tmp;
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar& operator+=(const double& src) {
    for (int i = 0; i < N; i++) {
      val[i] += 1.0 * (i + 1) * src;
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar& operator+=(const SuperScalar& src) {
    for (int i = 0; i < N; i++) {
      val[i] += src.val[i];
    }
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(const SuperScalar& src) const {
    bool compare = true;
    for (int i = 0; i < N; i++) {
      compare = compare && (val[i] == src.val[i]);
    }
    return compare;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const SuperScalar& src) const {
    bool compare = true;
    for (int i = 0; i < N; i++) {
      compare = compare && (val[i] == src.val[i]);
    }
    return !compare;
  }

  KOKKOS_INLINE_FUNCTION
  SuperScalar(const double& src) {
    for (int i = 0; i < N; i++) {
      val[i] = 1.0 * (i + 1) * src;
    }
  }
};

template <int N>
std::ostream& operator<<(std::ostream& os, const SuperScalar<N>& dt) {
  os << "{ ";
  for (int i = 0; i < N - 1; i++) {
    os << dt.val[i] << ", ";
  }
  os << dt.val[N - 1] << "}";

  return os;
}

template <class T, class DEVICE_TYPE>
struct ZeroFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = typename Kokkos::View<T, execution_space>;
  using h_type          = typename Kokkos::View<T, execution_space>::HostMirror;

  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const { data() = 0; }
};

//---------------------------------------------------
//--------------atomic_fetch_add---------------------
//---------------------------------------------------

template <class T, class DEVICE_TYPE>
struct AddFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;

  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const { Kokkos::atomic_fetch_add(&data(), (T)1); }
};

template <class T, class DEVICE_TYPE>
struct AddFunctorReduce {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;

  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int, int&) const { Kokkos::atomic_fetch_add(&data(), (T)1); }
};

template <class T, class execution_space>
T AddLoop(int loop) {
  struct ZeroFunctor<T, execution_space> f_zero;
  typename ZeroFunctor<T, execution_space>::type data("Data");
  typename ZeroFunctor<T, execution_space>::h_type h_data("HData");

  f_zero.data = data;

  Kokkos::parallel_for(1, f_zero);
  execution_space().fence();

  struct AddFunctor<T, execution_space> f_add;

  f_add.data = data;
  Kokkos::parallel_for(loop, f_add);
  execution_space().fence();

  Kokkos::deep_copy(h_data, data);
  T val = h_data();

  struct AddFunctorReduce<T, execution_space> f_add_red;
  f_add_red.data = data;
  int dummy_result;
  Kokkos::parallel_reduce(loop, f_add_red, dummy_result);
  execution_space().fence();

  return val;
}

template <class T>
T AddLoopSerial(int loop) {
  T* data = new T[1];
  data[0] = 0;

  for (int i = 0; i < loop; i++) {
    *data += (T)1;
  }

  T val = *data;
  delete[] data;

  return val;
}

//------------------------------------------------------
//--------------atomic_compare_exchange-----------------
//------------------------------------------------------

template <class T, class DEVICE_TYPE>
struct CASFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;

  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const {
    T old = data();
    T newval, assumed;

    do {
      assumed = old;
      newval  = assumed + (T)1;
      old     = Kokkos::atomic_compare_exchange(&data(), assumed, newval);
    } while (old != assumed);
  }
};

template <class T, class DEVICE_TYPE>
struct CASFunctorReduce {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;

  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int, int&) const {
    T old = data();
    T newval, assumed;

    do {
      assumed = old;
      newval  = assumed + (T)1;
      old     = Kokkos::atomic_compare_exchange(&data(), assumed, newval);
    } while (old != assumed);
  }
};

template <class T, class execution_space>
T CASLoop(int loop) {
  struct ZeroFunctor<T, execution_space> f_zero;
  typename ZeroFunctor<T, execution_space>::type data("Data");
  typename ZeroFunctor<T, execution_space>::h_type h_data("HData");

  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  execution_space().fence();

  struct CASFunctor<T, execution_space> f_cas;
  f_cas.data = data;
  Kokkos::parallel_for(loop, f_cas);
  execution_space().fence();

  Kokkos::deep_copy(h_data, data);
  T val = h_data();

  struct CASFunctorReduce<T, execution_space> f_cas_red;
  f_cas_red.data = data;
  int dummy_result;
  Kokkos::parallel_reduce(loop, f_cas_red, dummy_result);
  execution_space().fence();

  return val;
}

template <class T>
T CASLoopSerial(int loop) {
  T* data = new T[1];
  data[0] = 0;

  for (int i = 0; i < loop; i++) {
    T assumed;
    T newval;
    T old;

    do {
      assumed = *data;
      newval  = assumed + (T)1;
      old     = *data;
      *data   = newval;
    } while (!(assumed == old));
  }

  T val = *data;
  delete[] data;

  return val;
}

//----------------------------------------------
//--------------atomic_exchange-----------------
//----------------------------------------------

template <class T, class DEVICE_TYPE>
struct ExchFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;

  type data, data2;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    T old = Kokkos::atomic_exchange(&data(), (T)i);
    Kokkos::atomic_fetch_add(&data2(), old);
  }
};

template <class T, class DEVICE_TYPE>
struct ExchFunctorReduce {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;

  type data, data2;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, int&) const {
    T old = Kokkos::atomic_exchange(&data(), (T)i);
    Kokkos::atomic_fetch_add(&data2(), old);
  }
};

template <class T, class execution_space>
T ExchLoop(int loop) {
  struct ZeroFunctor<T, execution_space> f_zero;
  typename ZeroFunctor<T, execution_space>::type data("Data");
  typename ZeroFunctor<T, execution_space>::h_type h_data("HData");

  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  execution_space().fence();

  typename ZeroFunctor<T, execution_space>::type data2("Data");
  typename ZeroFunctor<T, execution_space>::h_type h_data2("HData");

  f_zero.data = data2;
  Kokkos::parallel_for(1, f_zero);
  execution_space().fence();

  struct ExchFunctor<T, execution_space> f_exch;
  f_exch.data  = data;
  f_exch.data2 = data2;
  Kokkos::parallel_for(loop, f_exch);
  execution_space().fence();

  Kokkos::deep_copy(h_data, data);
  Kokkos::deep_copy(h_data2, data2);
  T val = h_data() + h_data2();

  struct ExchFunctorReduce<T, execution_space> f_exch_red;
  f_exch_red.data  = data;
  f_exch_red.data2 = data2;
  int dummy_result;
  Kokkos::parallel_reduce(loop, f_exch_red, dummy_result);
  execution_space().fence();

  return val;
}

template <class T>
T ExchLoopSerial(
    typename std::conditional<!std::is_same<T, Kokkos::complex<double> >::value,
                              int, void>::type loop) {
  T* data  = new T[1];
  T* data2 = new T[1];
  data[0]  = 0;
  data2[0] = 0;

  for (int i = 0; i < loop; i++) {
    T old = *data;
    *data = (T)i;
    *data2 += old;
  }

  T val = *data2 + *data;
  delete[] data;
  delete[] data2;

  return val;
}

template <class T>
T ExchLoopSerial(
    typename std::conditional<std::is_same<T, Kokkos::complex<double> >::value,
                              int, void>::type loop) {
  T* data  = new T[1];
  T* data2 = new T[1];
  data[0]  = 0;
  data2[0] = 0;

  for (int i = 0; i < loop; i++) {
    T old        = *data;
    data->real() = (static_cast<double>(i));
    data->imag() = 0;
    *data2 += old;
  }

  T val = *data2 + *data;
  delete[] data;
  delete[] data2;

  return val;
}

template <class T, class DeviceType>
T LoopVariant(int loop, int test) {
  switch (test) {
    case 1: return AddLoop<T, DeviceType>(loop);
    case 2: return CASLoop<T, DeviceType>(loop);
    case 3: return ExchLoop<T, DeviceType>(loop);
  }

  return 0;
}

template <class T>
T LoopVariantSerial(int loop, int test) {
  switch (test) {
    case 1: return AddLoopSerial<T>(loop);
    case 2: return CASLoopSerial<T>(loop);
    case 3: return ExchLoopSerial<T>(loop);
  }

  return 0;
}

template <class T, class DeviceType>
bool Loop(int loop, int test) {
  T res       = LoopVariant<T, DeviceType>(loop, test);
  T resSerial = LoopVariantSerial<T>(loop, test);

  bool passed = true;

  if (resSerial != res) {
    passed = false;

    std::cout << "Loop<" << typeid(T).name() << ">( test = " << test
              << " FAILED : " << resSerial << " != " << res << std::endl;
  }

  return passed;
}

}  // namespace TestAtomic

namespace Test {

TEST(TEST_CATEGORY, atomics) {
  const int loop_count = 1e4;

  ASSERT_TRUE((TestAtomic::Loop<int, TEST_EXECSPACE>(loop_count, 1)));
  ASSERT_TRUE((TestAtomic::Loop<int, TEST_EXECSPACE>(loop_count, 2)));
  ASSERT_TRUE((TestAtomic::Loop<int, TEST_EXECSPACE>(loop_count, 3)));

  ASSERT_TRUE((TestAtomic::Loop<unsigned int, TEST_EXECSPACE>(loop_count, 1)));
  ASSERT_TRUE((TestAtomic::Loop<unsigned int, TEST_EXECSPACE>(loop_count, 2)));
  ASSERT_TRUE((TestAtomic::Loop<unsigned int, TEST_EXECSPACE>(loop_count, 3)));

  ASSERT_TRUE((TestAtomic::Loop<long int, TEST_EXECSPACE>(loop_count, 1)));
  ASSERT_TRUE((TestAtomic::Loop<long int, TEST_EXECSPACE>(loop_count, 2)));
  ASSERT_TRUE((TestAtomic::Loop<long int, TEST_EXECSPACE>(loop_count, 3)));

  ASSERT_TRUE(
      (TestAtomic::Loop<unsigned long int, TEST_EXECSPACE>(loop_count, 1)));
  ASSERT_TRUE(
      (TestAtomic::Loop<unsigned long int, TEST_EXECSPACE>(loop_count, 2)));
  ASSERT_TRUE(
      (TestAtomic::Loop<unsigned long int, TEST_EXECSPACE>(loop_count, 3)));

  ASSERT_TRUE((TestAtomic::Loop<long long int, TEST_EXECSPACE>(loop_count, 1)));
  ASSERT_TRUE((TestAtomic::Loop<long long int, TEST_EXECSPACE>(loop_count, 2)));
  ASSERT_TRUE((TestAtomic::Loop<long long int, TEST_EXECSPACE>(loop_count, 3)));

  ASSERT_TRUE((TestAtomic::Loop<double, TEST_EXECSPACE>(loop_count, 1)));
  ASSERT_TRUE((TestAtomic::Loop<double, TEST_EXECSPACE>(loop_count, 2)));
  ASSERT_TRUE((TestAtomic::Loop<double, TEST_EXECSPACE>(loop_count, 3)));

  ASSERT_TRUE((TestAtomic::Loop<float, TEST_EXECSPACE>(100, 1)));
  ASSERT_TRUE((TestAtomic::Loop<float, TEST_EXECSPACE>(100, 2)));
  ASSERT_TRUE((TestAtomic::Loop<float, TEST_EXECSPACE>(100, 3)));

#ifndef KOKKOS_ENABLE_OPENMPTARGET
  ASSERT_TRUE((TestAtomic::Loop<Kokkos::complex<float>, TEST_EXECSPACE>(1, 1)));
  ASSERT_TRUE((TestAtomic::Loop<Kokkos::complex<float>, TEST_EXECSPACE>(1, 2)));
  ASSERT_TRUE((TestAtomic::Loop<Kokkos::complex<float>, TEST_EXECSPACE>(1, 3)));

  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<float>, TEST_EXECSPACE>(100, 1)));
  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<float>, TEST_EXECSPACE>(100, 2)));
  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<float>, TEST_EXECSPACE>(100, 3)));

// FIXME_SYCL atomics for large types to be implemented
#ifndef KOKKOS_ENABLE_SYCL
  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<double>, TEST_EXECSPACE>(1, 1)));
  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<double>, TEST_EXECSPACE>(1, 2)));
  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<double>, TEST_EXECSPACE>(1, 3)));

  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<double>, TEST_EXECSPACE>(100, 1)));
  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<double>, TEST_EXECSPACE>(100, 2)));
  ASSERT_TRUE(
      (TestAtomic::Loop<Kokkos::complex<double>, TEST_EXECSPACE>(100, 3)));

// WORKAROUND MSVC
#ifndef _WIN32
  ASSERT_TRUE(
      (TestAtomic::Loop<TestAtomic::SuperScalar<4>, TEST_EXECSPACE>(100, 1)));
  ASSERT_TRUE(
      (TestAtomic::Loop<TestAtomic::SuperScalar<4>, TEST_EXECSPACE>(100, 2)));
  ASSERT_TRUE(
      (TestAtomic::Loop<TestAtomic::SuperScalar<4>, TEST_EXECSPACE>(100, 3)));
#endif
#endif
#endif
}

}  // namespace Test
