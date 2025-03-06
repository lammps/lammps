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

#include <Kokkos_Core.hpp>
#include <benchmark/benchmark.h>
#include "PerfTest_Category.hpp"

#include <cmath>
#include "PerfTestBlasKernels.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Test {

// Reduction   : result = dot( Q(:,j) , Q(:,j) );
// PostProcess : R(j,j) = result ; inv = 1 / result ;
template <class VectorView, class ValueView>
struct InvNorm2 : public Kokkos::DotSingle<VectorView> {
  using value_type = typename Kokkos::DotSingle<VectorView>::value_type;

  ValueView Rjj;
  ValueView inv;

  InvNorm2(const VectorView& argX, const ValueView& argR,
           const ValueView& argInv)
      : Kokkos::DotSingle<VectorView>(argX), Rjj(argR), inv(argInv) {}

  KOKKOS_INLINE_FUNCTION
  void final(value_type& result) const {
    result = Kokkos::sqrt(result);
    Rjj()  = result;
    inv()  = (0 < result) ? 1.0 / result : 0;
  }
};

template <class VectorView, class ValueView>
inline void invnorm2(const VectorView& x, const ValueView& r,
                     const ValueView& r_inv) {
  Kokkos::parallel_reduce(x.extent(0),
                          InvNorm2<VectorView, ValueView>(x, r, r_inv));
}

// PostProcess : tmp = - ( R(j,k) = result );
template <class VectorView, class ValueView>
struct DotM : public Kokkos::Dot<VectorView> {
  using value_type = typename Kokkos::Dot<VectorView>::value_type;

  ValueView Rjk;
  ValueView tmp;

  DotM(const VectorView& argX, const VectorView& argY, const ValueView& argR,
       const ValueView& argTmp)
      : Kokkos::Dot<VectorView>(argX, argY), Rjk(argR), tmp(argTmp) {}

  KOKKOS_INLINE_FUNCTION
  void final(value_type& result) const {
    Rjk() = result;
    tmp() = -result;
  }
};

template <class VectorView, class ValueView>
inline void dot_neg(const VectorView& x, const VectorView& y,
                    const ValueView& r, const ValueView& r_neg) {
  Kokkos::parallel_reduce(x.extent(0),
                          DotM<VectorView, ValueView>(x, y, r, r_neg));
}

template <typename Scalar, class DeviceType>
struct ModifiedGramSchmidt {
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  using multivector_type =
      Kokkos::View<Scalar**, Kokkos::LayoutLeft, execution_space>;

  using vector_type =
      Kokkos::View<Scalar*, Kokkos::LayoutLeft, execution_space>;

  using value_view = Kokkos::View<Scalar, Kokkos::LayoutLeft, execution_space>;

  multivector_type Q;
  multivector_type R;

  static double factorization(const multivector_type Q_,
                              const multivector_type R_) {
    const size_type count = Q_.extent(1);
    value_view tmp("tmp");
    value_view one("one");

    Kokkos::deep_copy(one, (Scalar)1);

    Kokkos::Timer timer;

    for (size_type j = 0; j < count; ++j) {
      // Reduction   : tmp = dot( Q(:,j) , Q(:,j) );
      // PostProcess : tmp = std::sqrt( tmp ); R(j,j) = tmp ; tmp = 1 / tmp ;
      const vector_type Qj = Kokkos::subview(Q_, Kokkos::ALL(), j);
      const value_view Rjj = Kokkos::subview(R_, j, j);

      invnorm2(Qj, Rjj, tmp);

      // Q(:,j) *= ( 1 / R(j,j) ); => Q(:,j) *= tmp ;
      Kokkos::scale(tmp, Qj);

      for (size_type k = j + 1; k < count; ++k) {
        const vector_type Qk = Kokkos::subview(Q_, Kokkos::ALL(), k);
        const value_view Rjk = Kokkos::subview(R_, j, k);

        // Reduction   : R(j,k) = dot( Q(:,j) , Q(:,k) );
        // PostProcess : tmp = - R(j,k);
        dot_neg(Qj, Qk, Rjk, tmp);

        // Q(:,k) -= R(j,k) * Q(:,j); => Q(:,k) += tmp * Q(:,j)
        Kokkos::axpby(tmp, Qj, one, Qk);
      }
    }

    execution_space().fence();

    return timer.seconds();
  }

  //--------------------------------------------------------------------------

  static double test(const size_type length, const size_type count) {
    multivector_type Q_("Q", length, count);
    multivector_type R_("R", count, count);

    typename multivector_type::HostMirror A = Kokkos::create_mirror(Q_);

    // Create and fill A on the host
    for (size_type j = 0; j < count; ++j) {
      for (size_type i = 0; i < length; ++i) {
        A(i, j) = (i + 1) * (j + 1);
      }
    }

    Kokkos::deep_copy(Q_, A);

    // A = Q * R
    const double dt = factorization(Q_, R_);

    return dt;
  }
};

template <class Scalar>
static void GramSchmidt(benchmark::State& state) {
  const int parallel_work_length = state.range(0);

  for (auto _ : state) {
    const double seconds =
        ModifiedGramSchmidt<Scalar, Kokkos::DefaultExecutionSpace>::test(
            parallel_work_length, 32);

    state.SetIterationTime(seconds);
    state.counters["Count"] = benchmark::Counter(parallel_work_length);
    state.counters["Time normalized"] =
        benchmark::Counter(seconds / parallel_work_length);
  }
}

// FIXME_SYCL SYCL+Cuda reports "an illegal memory access was encountered"
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
BENCHMARK(GramSchmidt<double>)
    ->ArgName("Count")
    ->ArgsProduct({
        benchmark::CreateRange(1 << 10, 1 << 18, 2),
    })
    ->UseManualTime()
    ->Iterations(5);
#else
BENCHMARK(GramSchmidt<double>)
    ->ArgName("Count")
    ->ArgsProduct({
        benchmark::CreateRange(1 << 10, 1 << 19, 2),
    })
    ->UseManualTime()
    ->Iterations(5);
#endif

}  // namespace Test
