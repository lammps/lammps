// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOSTEST_UTILS_HPP
#define KOKKOSTEST_UTILS_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <type_traits>

namespace KokkosTest {

struct FloatingPointComparison {
 private:
  template <class T>
  KOKKOS_FUNCTION static double eps(T) {
    return DBL_EPSILON;
  }

#if !KOKKOS_HALF_T_IS_FLOAT
  KOKKOS_FUNCTION static Kokkos::Experimental::half_t eps(
      Kokkos::Experimental::half_t) {
// FIXME_NVHPC compile-time error
#ifdef KOKKOS_COMPILER_NVHPC
    return 0.0009765625F;
#else
    return Kokkos::Experimental::epsilon<Kokkos::Experimental::half_t>::value;
#endif
  }
#endif

#if !KOKKOS_BHALF_T_IS_FLOAT
  KOKKOS_FUNCTION
  Kokkos::Experimental::bhalf_t static eps(Kokkos::Experimental::bhalf_t) {
// FIXME_NVHPC compile-time error
#ifdef KOKKOS_COMPILER_NVHPC
    return 0.0078125;
#else
    return Kokkos::Experimental::epsilon<Kokkos::Experimental::bhalf_t>::value;
#endif
  }
#endif

  KOKKOS_FUNCTION static double eps(float) { return FLT_EPSILON; }

// POWER9 gives unexpected values with LDBL_EPSILON issues
// https://stackoverflow.com/questions/68960416/ppc64-long-doubles-machine-epsilon-calculation
#if defined(KOKKOS_ARCH_POWER9) || defined(KOKKOS_ARCH_POWER8)
  KOKKOS_FUNCTION static double eps(long double) { return DBL_EPSILON; }
#else
  KOKKOS_FUNCTION static double eps(long double) { return LDBL_EPSILON; }
#endif

  // Using absolute here instead of abs, since we actually test abs ...
  template <class T>
  KOKKOS_FUNCTION static T absolute(T val) {
    if constexpr (std::is_signed_v<T>) {
      return val < T(0) ? -val : val;
    }
    return val;
  }

 public:
  template <class FPT>
  KOKKOS_FUNCTION static bool compare_near_zero(FPT const& fpv, int ulp) {
    auto abs_tol = eps(fpv) * ulp;

    bool ar = absolute(fpv) <= abs_tol;
    if (!ar) {
      Kokkos::printf("absolute value exceeds tolerance [|%e| > %e]\n",
                     (double)fpv, (double)abs_tol);
    }

    return ar;
  }

  template <class Lhs, class Rhs>
  KOKKOS_FUNCTION static bool compare(Lhs const& lhs, Rhs const& rhs, int ulp) {
    if (lhs == 0) {
      return compare_near_zero(rhs, ulp);
    } else if (rhs == 0) {
      return compare_near_zero(lhs, ulp);
    } else {
      auto rel_tol     = (eps(lhs) < eps(rhs) ? eps(lhs) : eps(rhs)) * ulp;
      double abs_diff  = static_cast<double>(rhs > lhs ? rhs - lhs : lhs - rhs);
      double min_denom = static_cast<double>(
          absolute(rhs) < absolute(lhs) ? absolute(rhs) : absolute(lhs));
      double rel_diff = abs_diff / min_denom;
      bool ar         = rel_diff <= rel_tol;
      if (!ar) {
        Kokkos::printf("relative difference exceeds tolerance [%e > %e]\n",
                       (double)rel_diff, (double)rel_tol);
      }

      return ar;
    }
  }
};

struct IntegerComparison {
  template <class Lhs, class Rhs>
  KOKKOS_FUNCTION bool compare(Lhs const& lhs, Rhs const& rhs) const {
    static_assert(std::is_integral_v<Lhs>);
    static_assert(std::is_integral_v<Rhs>);
    return lhs == rhs;
  }
};

}  // namespace KokkosTest
#endif
