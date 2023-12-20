//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <sstream>
#include <iostream>
#include <limits>

#include <Kokkos_Core.hpp>

namespace Test {

struct ReducerTag {};

template <typename ScalarType, class DeviceType>
class ReduceFunctor {
 public:
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  struct value_type {
    ScalarType value[3];
  };

  const size_type nwork;

  KOKKOS_INLINE_FUNCTION
  ReduceFunctor(const size_type& arg_nwork) : nwork(arg_nwork) {}

  KOKKOS_INLINE_FUNCTION
  ReduceFunctor(const ReduceFunctor& rhs) : nwork(rhs.nwork) {}

  /*
    KOKKOS_INLINE_FUNCTION
    void init( value_type & dst ) const
    {
      dst.value[0] = 0;
      dst.value[1] = 0;
      dst.value[2] = 0;
    }
  */

  KOKKOS_INLINE_FUNCTION
  void join(value_type& dst, const value_type& src) const {
    dst.value[0] += src.value[0];
    dst.value[1] += src.value[1];
    dst.value[2] += src.value[2];
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type iwork, value_type& dst) const {
    dst.value[0] += 1;
    dst.value[1] += iwork + 1;
    dst.value[2] += nwork - iwork;
  }
};

template <class DeviceType>
class ReduceFunctorFinal : public ReduceFunctor<int64_t, DeviceType> {
 public:
  using value_type = typename ReduceFunctor<int64_t, DeviceType>::value_type;

  KOKKOS_INLINE_FUNCTION
  ReduceFunctorFinal(const size_t n) : ReduceFunctor<int64_t, DeviceType>(n) {}

  KOKKOS_INLINE_FUNCTION
  void final(value_type& dst) const {
    dst.value[0] = -dst.value[0];
    dst.value[1] = -dst.value[1];
    dst.value[2] = -dst.value[2];
  }
};

template <class DeviceType>
class ReduceFunctorFinalTag {
 public:
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;
  using ScalarType      = int64_t;

  struct value_type {
    ScalarType value[3];
  };

  const size_type nwork;

  KOKKOS_INLINE_FUNCTION
  ReduceFunctorFinalTag(const size_type arg_nwork) : nwork(arg_nwork) {}

  KOKKOS_INLINE_FUNCTION
  void join(const ReducerTag, value_type& dst, const value_type& src) const {
    dst.value[0] += src.value[0];
    dst.value[1] += src.value[1];
    dst.value[2] += src.value[2];
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ReducerTag, size_type iwork, value_type& dst) const {
    dst.value[0] -= 1;
    dst.value[1] -= iwork + 1;
    dst.value[2] -= nwork - iwork;
  }

  KOKKOS_INLINE_FUNCTION
  void final(const ReducerTag, value_type& dst) const {
    ++dst.value[0];
    ++dst.value[1];
    ++dst.value[2];
  }
};

template <typename ScalarType, class DeviceType>
class RuntimeReduceFunctor {
 public:
  // Required for functor:
  using execution_space = DeviceType;
  using value_type      = ScalarType[];
  const unsigned value_count;

  // Unit test details:

  using size_type = typename execution_space::size_type;

  const size_type nwork;

  RuntimeReduceFunctor(const size_type arg_nwork, const size_type arg_count)
      : value_count(arg_count), nwork(arg_nwork) {}

  KOKKOS_INLINE_FUNCTION
  void init(ScalarType dst[]) const {
    for (unsigned i = 0; i < value_count; ++i) dst[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join(ScalarType dst[], const ScalarType src[]) const {
    for (unsigned i = 0; i < value_count; ++i) dst[i] += src[i];
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type iwork, ScalarType dst[]) const {
    const size_type tmp[3] = {1, iwork + 1, nwork - iwork};

    for (size_type i = 0; i < static_cast<size_type>(value_count); ++i) {
      dst[i] += tmp[i % 3];
    }
  }
};

template <typename ScalarType, class DeviceType>
class RuntimeReduceMinMax {
 public:
  // Required for functor:
  using execution_space = DeviceType;
  using value_type      = ScalarType[];
  const unsigned value_count;

  // Unit test details:

  using size_type = typename execution_space::size_type;

  const size_type nwork;
  const ScalarType amin;
  const ScalarType amax;

  RuntimeReduceMinMax(const size_type arg_nwork, const size_type arg_count)
      : value_count(arg_count),
        nwork(arg_nwork),
        amin(std::numeric_limits<ScalarType>::min()),
        amax(std::numeric_limits<ScalarType>::max()) {}

  KOKKOS_INLINE_FUNCTION
  void init(ScalarType dst[]) const {
    for (unsigned i = 0; i < value_count; ++i) {
      dst[i] = i % 2 ? amax : amin;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join(ScalarType dst[], const ScalarType src[]) const {
    for (unsigned i = 0; i < value_count; ++i) {
      dst[i] = i % 2 ? (dst[i] < src[i] ? dst[i] : src[i])   // min
                     : (dst[i] > src[i] ? dst[i] : src[i]);  // max
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type iwork, ScalarType dst[]) const {
    const ScalarType tmp[2] = {ScalarType(iwork + 1),
                               ScalarType(nwork - iwork)};

    for (size_type i = 0; i < static_cast<size_type>(value_count); ++i) {
      dst[i] = i % 2 ? (dst[i] < tmp[i % 2] ? dst[i] : tmp[i % 2])
                     : (dst[i] > tmp[i % 2] ? dst[i] : tmp[i % 2]);
    }
  }
};

template <class DeviceType>
class RuntimeReduceFunctorFinal
    : public RuntimeReduceFunctor<int64_t, DeviceType> {
 public:
  using base_type   = RuntimeReduceFunctor<int64_t, DeviceType>;
  using value_type  = typename base_type::value_type;
  using scalar_type = int64_t;

  RuntimeReduceFunctorFinal(const size_t theNwork, const size_t count)
      : base_type(theNwork, count) {}

  KOKKOS_INLINE_FUNCTION
  void final(value_type dst) const {
    for (unsigned i = 0; i < base_type::value_count; ++i) {
      dst[i] = -dst[i];
    }
  }
};

template <class ValueType, class DeviceType>
class CombinedReduceFunctorSameType {
 public:
  using execution_space = typename DeviceType::execution_space;
  using size_type       = typename execution_space::size_type;

  const size_type nwork;

  KOKKOS_INLINE_FUNCTION
  constexpr explicit CombinedReduceFunctorSameType(const size_type& arg_nwork)
      : nwork(arg_nwork) {}

  KOKKOS_DEFAULTED_FUNCTION
  constexpr CombinedReduceFunctorSameType(
      const CombinedReduceFunctorSameType& rhs) = default;

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type iwork, ValueType& dst1, ValueType& dst2,
                  ValueType& dst3) const {
    dst1 += 1;
    dst2 += iwork + 1;
    dst3 += nwork - iwork;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type iwork, size_type always_zero_1,
                  size_type always_zero_2, ValueType& dst1, ValueType& dst2,
                  ValueType& dst3) const {
    dst1 += 1 + always_zero_1;
    dst2 += iwork + 1 + always_zero_2;
    dst3 += nwork - iwork;
  }
};

namespace {

template <typename ScalarType, class DeviceType>
class TestReduce {
 public:
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  TestReduce(const size_type& nwork) {
    run_test(nwork);
    run_test_final(nwork);
    run_test_final_tag(nwork);
  }

  void run_test(const size_type& nwork) {
    using functor_type = Test::ReduceFunctor<ScalarType, execution_space>;
    using value_type   = typename functor_type::value_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    value_type result[Repeat];

    const uint64_t nw   = nwork;
    const uint64_t nsum = nw % 2 ? nw * ((nw + 1) / 2) : (nw / 2) * (nw + 1);

    for (unsigned i = 0; i < Repeat; ++i) {
      Kokkos::parallel_reduce(nwork, functor_type(nwork), result[i]);
    }

    for (unsigned i = 0; i < Repeat; ++i) {
      for (unsigned j = 0; j < Count; ++j) {
        const uint64_t correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ((ScalarType)correct, result[i].value[j]);
      }
    }
  }

  void run_test_final(const size_type& nwork) {
    using functor_type = Test::ReduceFunctorFinal<execution_space>;
    using value_type   = typename functor_type::value_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    value_type result[Repeat];

    const uint64_t nw   = nwork;
    const uint64_t nsum = nw % 2 ? nw * ((nw + 1) / 2) : (nw / 2) * (nw + 1);

    for (unsigned i = 0; i < Repeat; ++i) {
      if (i % 2 == 0) {
        Kokkos::parallel_reduce(nwork, functor_type(nwork), result[i]);
      } else {
        Kokkos::parallel_reduce("Reduce", nwork, functor_type(nwork),
                                result[i]);
      }
    }

    for (unsigned i = 0; i < Repeat; ++i) {
      for (unsigned j = 0; j < Count; ++j) {
        const uint64_t correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ((ScalarType)correct, -result[i].value[j]);
      }
    }
  }

  void run_test_final_tag(const size_type& nwork) {
    using functor_type = Test::ReduceFunctorFinalTag<execution_space>;
    using value_type   = typename functor_type::value_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    value_type result[Repeat];

    const uint64_t nw   = nwork;
    const uint64_t nsum = nw % 2 ? nw * ((nw + 1) / 2) : (nw / 2) * (nw + 1);

    for (unsigned i = 0; i < Repeat; ++i) {
      if (i % 2 == 0) {
        Kokkos::parallel_reduce(
            Kokkos::RangePolicy<execution_space, ReducerTag>(0, nwork),
            functor_type(nwork), result[i]);
      } else {
        Kokkos::parallel_reduce(
            "Reduce",
            Kokkos::RangePolicy<execution_space, ReducerTag>(0, nwork),
            functor_type(nwork), result[i]);
      }
    }

    for (unsigned i = 0; i < Repeat; ++i) {
      for (unsigned j = 0; j < Count; ++j) {
        const uint64_t correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ((ScalarType)correct, 1 - result[i].value[j]);
      }
    }
  }
};

template <typename ScalarType, class DeviceType>
class TestReduceDynamic {
 public:
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  TestReduceDynamic(const size_type nwork) {
    run_test_dynamic(nwork);
#ifndef KOKKOS_ENABLE_OPENACC
    // FIXME_OPENACC - OpenACC (V3.3) does not support custom reductions.
    run_test_dynamic_minmax(nwork);
#endif
    run_test_dynamic_final(nwork);
  }

  void run_test_dynamic(const size_type nwork) {
    using functor_type =
        Test::RuntimeReduceFunctor<ScalarType, execution_space>;

    enum { Count = 3 };
    enum { Repeat = 100 };

    ScalarType result[Repeat][Count];

    const uint64_t nw   = nwork;
    const uint64_t nsum = nw % 2 ? nw * ((nw + 1) / 2) : (nw / 2) * (nw + 1);

    for (unsigned i = 0; i < Repeat; ++i) {
      if (i % 2 == 0) {
        Kokkos::parallel_reduce(nwork, functor_type(nwork, Count), result[i]);
      } else {
        Kokkos::parallel_reduce("Reduce", nwork, functor_type(nwork, Count),
                                result[i]);
      }
    }

    for (unsigned i = 0; i < Repeat; ++i) {
      for (unsigned j = 0; j < Count; ++j) {
        const uint64_t correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ((ScalarType)correct, result[i][j]);
      }
    }
  }

  void run_test_dynamic_minmax(const size_type nwork) {
    using functor_type = Test::RuntimeReduceMinMax<ScalarType, execution_space>;

    enum { Count = 2 };
    enum { Repeat = 100 };

    ScalarType result[Repeat][Count];

    for (unsigned i = 0; i < Repeat; ++i) {
      if (i % 2 == 0) {
        Kokkos::parallel_reduce(nwork, functor_type(nwork, Count), result[i]);
      } else {
        Kokkos::parallel_reduce("Reduce", nwork, functor_type(nwork, Count),
                                result[i]);
      }
    }

    for (unsigned i = 0; i < Repeat; ++i) {
      for (unsigned j = 0; j < Count; ++j) {
        if (nwork == 0) {
          ScalarType amin(std::numeric_limits<ScalarType>::min());
          ScalarType amax(std::numeric_limits<ScalarType>::max());
          const ScalarType correct = (j % 2) ? amax : amin;
          ASSERT_EQ((ScalarType)correct, result[i][j]);
        } else {
          const uint64_t correct = j % 2 ? 1 : nwork;
          ASSERT_EQ((ScalarType)correct, result[i][j]);
        }
      }
    }
  }

  void run_test_dynamic_final(const size_type nwork) {
    using functor_type = Test::RuntimeReduceFunctorFinal<execution_space>;

    enum { Count = 3 };
    enum { Repeat = 100 };

    typename functor_type::scalar_type result[Repeat][Count];

    const uint64_t nw   = nwork;
    const uint64_t nsum = nw % 2 ? nw * ((nw + 1) / 2) : (nw / 2) * (nw + 1);

    for (unsigned i = 0; i < Repeat; ++i) {
      if (i % 2 == 0) {
        Kokkos::parallel_reduce(nwork, functor_type(nwork, Count), result[i]);
      } else {
        Kokkos::parallel_reduce("TestKernelReduce", nwork,
                                functor_type(nwork, Count), result[i]);
      }
    }

    for (unsigned i = 0; i < Repeat; ++i) {
      for (unsigned j = 0; j < Count; ++j) {
        const uint64_t correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ((ScalarType)correct, -result[i][j]);
      }
    }
  }
};

template <typename ScalarType, class DeviceType>
class TestReduceDynamicView {
 public:
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  TestReduceDynamicView(const size_type nwork) { run_test_dynamic_view(nwork); }

  void run_test_dynamic_view(const size_type nwork) {
    using functor_type =
        Test::RuntimeReduceFunctor<ScalarType, execution_space>;

    using result_type      = Kokkos::View<ScalarType*, DeviceType>;
    using result_host_type = typename result_type::HostMirror;

    const unsigned CountLimit = 23;

    const uint64_t nw   = nwork;
    const uint64_t nsum = nw % 2 ? nw * ((nw + 1) / 2) : (nw / 2) * (nw + 1);

    for (unsigned count = 0; count < CountLimit; ++count) {
      result_type result("result", count);
      result_host_type host_result = Kokkos::create_mirror(result);

      // Test result to host pointer:

      std::string str("TestKernelReduce");
      if (count % 2 == 0) {
        Kokkos::parallel_reduce(nw, functor_type(nw, count), host_result);
      } else {
        Kokkos::parallel_reduce(str, nw, functor_type(nw, count), host_result);
      }
      Kokkos::fence("Fence before accessing result on the host");

      for (unsigned j = 0; j < count; ++j) {
        const uint64_t correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ(host_result(j), (ScalarType)correct);
        host_result(j) = 0;
      }
    }
  }
};

}  // namespace

// FIXME_SYCL
// FIXME_OPENMPTARGET : The feature works with LLVM/13 on NVIDIA
// architectures. The jenkins currently tests with LLVM/12.
#if !defined(KOKKOS_ENABLE_SYCL) &&          \
    (!defined(KOKKOS_ENABLE_OPENMPTARGET) || \
     defined(KOKKOS_COMPILER_CLANG) && (KOKKOS_COMPILER_CLANG >= 1300))
TEST(TEST_CATEGORY, int64_t_reduce) {
  TestReduce<int64_t, TEST_EXECSPACE>(0);
  TestReduce<int64_t, TEST_EXECSPACE>(1000000);
}

TEST(TEST_CATEGORY, double_reduce) {
  TestReduce<double, TEST_EXECSPACE>(0);
  TestReduce<double, TEST_EXECSPACE>(1000000);
}

TEST(TEST_CATEGORY, int64_t_reduce_dynamic) {
  TestReduceDynamic<int64_t, TEST_EXECSPACE>(0);
  TestReduceDynamic<int64_t, TEST_EXECSPACE>(1000000);
}

TEST(TEST_CATEGORY, double_reduce_dynamic) {
  TestReduceDynamic<double, TEST_EXECSPACE>(0);
  TestReduceDynamic<double, TEST_EXECSPACE>(1000000);
}

TEST(TEST_CATEGORY, int64_t_reduce_dynamic_view) {
  TestReduceDynamicView<int64_t, TEST_EXECSPACE>(0);
  TestReduceDynamicView<int64_t, TEST_EXECSPACE>(1000000);
}
#endif

// FIXME_OPENMPTARGET: Not yet implemented.
#ifndef KOKKOS_ENABLE_OPENMPTARGET
// FIXME_OPENACC: Not yet implemented.
#ifndef KOKKOS_ENABLE_OPENACC
TEST(TEST_CATEGORY, int_combined_reduce) {
  using functor_type = CombinedReduceFunctorSameType<int64_t, TEST_EXECSPACE>;
  constexpr uint64_t nw = 1000;

  uint64_t nsum = (nw / 2) * (nw + 1);

  int64_t result1 = 0;
  int64_t result2 = 0;
  int64_t result3 = 0;

  Kokkos::parallel_reduce("int_combined_reduce",
                          Kokkos::RangePolicy<TEST_EXECSPACE>(0, nw),
                          functor_type(nw), result1, result2, result3);

  ASSERT_EQ(nw, uint64_t(result1));
  ASSERT_EQ(nsum, uint64_t(result2));
  ASSERT_EQ(nsum, uint64_t(result3));
}

TEST(TEST_CATEGORY, mdrange_combined_reduce) {
  using functor_type = CombinedReduceFunctorSameType<int64_t, TEST_EXECSPACE>;
  constexpr uint64_t nw = 1000;

  uint64_t nsum = (nw / 2) * (nw + 1);

  int64_t result1 = 0;
  int64_t result2 = 0;
  int64_t result3 = 0;

  Kokkos::parallel_reduce(
      "int_combined_reduce_mdrange",
      Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<3>>({{0, 0, 0}},
                                                             {{nw, 1, 1}}),
      functor_type(nw), result1, result2, result3);

  ASSERT_EQ(nw, uint64_t(result1));
  ASSERT_EQ(nsum, uint64_t(result2));
  ASSERT_EQ(nsum, uint64_t(result3));
}

TEST(TEST_CATEGORY, int_combined_reduce_mixed) {
  using functor_type = CombinedReduceFunctorSameType<int64_t, TEST_EXECSPACE>;

  constexpr uint64_t nw = 1000;

  uint64_t nsum = (nw / 2) * (nw + 1);
  {
    auto result1_v  = Kokkos::View<int64_t, Kokkos::HostSpace>{"result1_v"};
    int64_t result2 = 0;
    auto result3_v  = Kokkos::View<int64_t, Kokkos::HostSpace>{"result3_v"};
    Kokkos::parallel_reduce("int_combined-reduce_mixed",
                            Kokkos::RangePolicy<TEST_EXECSPACE>(0, nw),
                            functor_type(nw), result1_v, result2,
                            Kokkos::Sum<int64_t, Kokkos::HostSpace>{result3_v});
    ASSERT_EQ(int64_t(nw), result1_v());
    ASSERT_EQ(int64_t(nsum), result2);
    ASSERT_EQ(int64_t(nsum), result3_v());
  }
  {
    using MemorySpace = typename TEST_EXECSPACE::memory_space;
    auto result1_v    = Kokkos::View<int64_t, MemorySpace>{"result1_v"};
    int64_t result2   = 0;
    auto result3_v    = Kokkos::View<int64_t, MemorySpace>{"result3_v"};
    Kokkos::parallel_reduce("int_combined-reduce_mixed",
                            Kokkos::RangePolicy<TEST_EXECSPACE>(0, nw),
                            functor_type(nw), result1_v, result2,
                            Kokkos::Sum<int64_t, MemorySpace>{result3_v});
    int64_t result1;
    Kokkos::deep_copy(result1, result1_v);
    ASSERT_EQ(int64_t(nw), result1);
    ASSERT_EQ(int64_t(nsum), result2);
    int64_t result3;
    Kokkos::deep_copy(result3, result3_v);
    ASSERT_EQ(int64_t(nsum), result3);
  }
}
#endif
#endif
}  // namespace Test
