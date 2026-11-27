// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif

constexpr int N = 3;

namespace {

struct ArrayLike {
  double data[N];
  KOKKOS_FUNCTION operator double*() { return data; }
  KOKKOS_FUNCTION operator const double*() const { return data; }
};

// clang-format off
// clang-format gets confused, and thinks the lambda closure brace is the namespace
void from_array_like() {
  Kokkos::View<ArrayLike, TEST_EXECSPACE> a("A");
  int errors = 0;
  Kokkos::parallel_reduce(Kokkos::RangePolicy(TEST_EXECSPACE(), 0, 1), KOKKOS_LAMBDA(int, int& err) {
    {
      Kokkos::View<double*, TEST_EXECSPACE> b(a(), N);
      if(b.data() != &a().data[0]) err++;
    }
    {
      Kokkos::View<double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged> b(a(), N);
      if (b.data() != &a().data[0]) err++;
    }
    {
      Kokkos::View<const double*, TEST_EXECSPACE> b(a(), N);
      if (b.data() != &a().data[0]) err++;
    }
    {
      Kokkos::View<const double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged> b(a(), N);
      if (b.data() != &a().data[0]) err++;
    }
  }, errors);
  ASSERT_EQ(errors, 0);
}
// clang-format on

void from_carray() {
  int errors = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy(TEST_EXECSPACE(), 0, 1),
      KOKKOS_LAMBDA(int, int& err) {
        double a[N];
        {
          Kokkos::View<double*, TEST_EXECSPACE> b(a, N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged> b(a,
                                                                           N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE> b(a, N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged>
              b(a, N);
          if (b.data() != &a[0]) err++;
        }
      },
      errors);
  ASSERT_EQ(errors, 0);
}

struct CustomPtr {
  double* ptr;
  KOKKOS_FUNCTION operator double*() const { return ptr; }
};

void from_custom_ptr() {
  int errors = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy(TEST_EXECSPACE(), 0, 1),
      KOKKOS_LAMBDA(int, int& err) {
        double data[N];
        CustomPtr a{&data[0]};
        {
          Kokkos::View<double*, TEST_EXECSPACE> b(a, N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged> b(a,
                                                                           N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE> b(a, N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged>
              b(a, N);
          if (b.data() != &a[0]) err++;
        }
      },
      errors);
  ASSERT_EQ(errors, 0);
}

void from_ptr() {
  int errors = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy(TEST_EXECSPACE(), 0, 1),
      KOKKOS_LAMBDA(int, int& err) {
        double a[N];
        {
          Kokkos::View<double*, TEST_EXECSPACE> b(&a[0], N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged> b(
              &a[0], N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE> b(&a[0], N);
          if (b.data() != &a[0]) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged>
              b(&a[0], N);
          if (b.data() != &a[0]) err++;
        }
      },
      errors);
  ASSERT_EQ(errors, 0);
}

KOKKOS_FUNCTION std::nullptr_t& nullptr_ref(std::nullptr_t& val) { return val; }

void from_nullptr_and_0() {
  int errors = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy(TEST_EXECSPACE(), 0, 1),
      KOKKOS_LAMBDA(int, int& err) {
        {
          Kokkos::View<double*, TEST_EXECSPACE> b(nullptr, 0);
          if (b.data() != nullptr) err++;
        }
        {
          Kokkos::View<double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged> b(
              nullptr, 0);
          if (b.data() != nullptr) err++;
        }
        // need a temporary so we can test reference to nullptr_t
        // you can't get a non-const ref to nullptr (the variable)
        std::nullptr_t nullptrtmp = {};
        {
          Kokkos::View<const double*, TEST_EXECSPACE> b(nullptr_ref(nullptrtmp),
                                                        0);
          if (b.data() != nullptr) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged>
              b(nullptr_ref(nullptrtmp), 0);
          if (b.data() != nullptr) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE> b(
              0, 0);  // NOLINT(modernize-use-nullptr)
          if (b.data() != nullptr) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged>
              b(0, 0);  // NOLINT(modernize-use-nullptr)
          if (b.data() != nullptr) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE> b(
              NULL, 0);  // NOLINT(modernize-use-nullptr)
          if (b.data() != nullptr) err++;
        }
        {
          Kokkos::View<const double*, TEST_EXECSPACE, Kokkos::MemoryUnmanaged>
              b(NULL, 0);  // NOLINT(modernize-use-nullptr)
          if (b.data() != nullptr) err++;
        }
      },
      errors);
  ASSERT_EQ(errors, 0);
}

TEST(TEST_CATEGORY, view_ctor_ptr_convertible) {
  // FIXME_OPENACC
#ifdef KOKKOS_ENABLE_OPENACC
  GTEST_SKIP() << "Known to fail with OpenACC";
  KOKKOS_IMPL_UNREACHABLE();
#endif
  from_array_like();
  from_carray();
  from_custom_ptr();
  from_ptr();
  from_nullptr_and_0();
}

}  // namespace
