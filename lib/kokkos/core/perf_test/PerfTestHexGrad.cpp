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
#include "Benchmark_Context.hpp"
#include "PerfTest_Category.hpp"

namespace Test {

template <class DeviceType, typename CoordScalarType = double,
          typename GradScalarType = float>
struct HexGrad {
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  using self_type = HexGrad<DeviceType, CoordScalarType, GradScalarType>;

  // 3D array : ( ParallelWork , Space , Node )

  enum { NSpace = 3, NNode = 8 };

  using elem_coord_type =
      Kokkos::View<CoordScalarType * [NSpace][NNode], execution_space>;

  using elem_grad_type =
      Kokkos::View<GradScalarType * [NSpace][NNode], execution_space>;

  elem_coord_type coords;
  elem_grad_type grad_op;

  enum { FLOPS = 318 };  // = 3 * ( 18 + 8 * 11 ) };
  enum { READS = 18 };
  enum { WRITES = 18 };

  HexGrad(const elem_coord_type& arg_coords, const elem_grad_type& arg_grad_op)
      : coords(arg_coords), grad_op(arg_grad_op) {}

  KOKKOS_INLINE_FUNCTION static void grad(const CoordScalarType x[],
                                          const CoordScalarType z[],
                                          GradScalarType grad_y[]) {
    const GradScalarType R42 = (x[3] - x[1]);
    const GradScalarType R52 = (x[4] - x[1]);
    const GradScalarType R54 = (x[4] - x[3]);

    const GradScalarType R63 = (x[5] - x[2]);
    const GradScalarType R83 = (x[7] - x[2]);
    const GradScalarType R86 = (x[7] - x[5]);

    const GradScalarType R31 = (x[2] - x[0]);
    const GradScalarType R61 = (x[5] - x[0]);
    const GradScalarType R74 = (x[6] - x[3]);

    const GradScalarType R72 = (x[6] - x[1]);
    const GradScalarType R75 = (x[6] - x[4]);
    const GradScalarType R81 = (x[7] - x[0]);

    const GradScalarType t1 = (R63 + R54);
    const GradScalarType t2 = (R61 + R74);
    const GradScalarType t3 = (R72 + R81);

    const GradScalarType t4 = (R86 + R42);
    const GradScalarType t5 = (R83 + R52);
    const GradScalarType t6 = (R75 + R31);

    //  Calculate Y gradient from X and Z data

    grad_y[0] = (z[1] * t1) - (z[2] * R42) - (z[3] * t5) + (z[4] * t4) +
                (z[5] * R52) - (z[7] * R54);
    grad_y[1] = (z[2] * t2) + (z[3] * R31) - (z[0] * t1) - (z[5] * t6) +
                (z[6] * R63) - (z[4] * R61);
    grad_y[2] = (z[3] * t3) + (z[0] * R42) - (z[1] * t2) - (z[6] * t4) +
                (z[7] * R74) - (z[5] * R72);
    grad_y[3] = (z[0] * t5) - (z[1] * R31) - (z[2] * t3) + (z[7] * t6) +
                (z[4] * R81) - (z[6] * R83);
    grad_y[4] = (z[5] * t3) + (z[6] * R86) - (z[7] * t2) - (z[0] * t4) -
                (z[3] * R81) + (z[1] * R61);
    grad_y[5] = (z[6] * t5) - (z[4] * t3) - (z[7] * R75) + (z[1] * t6) -
                (z[0] * R52) + (z[2] * R72);
    grad_y[6] = (z[7] * t1) - (z[5] * t5) - (z[4] * R86) + (z[2] * t4) -
                (z[1] * R63) + (z[3] * R83);
    grad_y[7] = (z[4] * t2) - (z[6] * t1) + (z[5] * R75) - (z[3] * t6) -
                (z[2] * R74) + (z[0] * R54);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type ielem) const {
    GradScalarType g[NNode];

    const CoordScalarType x[NNode] = {coords(ielem, 0, 0), coords(ielem, 0, 1),
                                      coords(ielem, 0, 2), coords(ielem, 0, 3),
                                      coords(ielem, 0, 4), coords(ielem, 0, 5),
                                      coords(ielem, 0, 6), coords(ielem, 0, 7)};

    const CoordScalarType y[NNode] = {coords(ielem, 1, 0), coords(ielem, 1, 1),
                                      coords(ielem, 1, 2), coords(ielem, 1, 3),
                                      coords(ielem, 1, 4), coords(ielem, 1, 5),
                                      coords(ielem, 1, 6), coords(ielem, 1, 7)};

    const CoordScalarType z[NNode] = {coords(ielem, 2, 0), coords(ielem, 2, 1),
                                      coords(ielem, 2, 2), coords(ielem, 2, 3),
                                      coords(ielem, 2, 4), coords(ielem, 2, 5),
                                      coords(ielem, 2, 6), coords(ielem, 2, 7)};

    grad(z, y, g);

    grad_op(ielem, 0, 0) = g[0];
    grad_op(ielem, 0, 1) = g[1];
    grad_op(ielem, 0, 2) = g[2];
    grad_op(ielem, 0, 3) = g[3];
    grad_op(ielem, 0, 4) = g[4];
    grad_op(ielem, 0, 5) = g[5];
    grad_op(ielem, 0, 6) = g[6];
    grad_op(ielem, 0, 7) = g[7];

    grad(x, z, g);

    grad_op(ielem, 1, 0) = g[0];
    grad_op(ielem, 1, 1) = g[1];
    grad_op(ielem, 1, 2) = g[2];
    grad_op(ielem, 1, 3) = g[3];
    grad_op(ielem, 1, 4) = g[4];
    grad_op(ielem, 1, 5) = g[5];
    grad_op(ielem, 1, 6) = g[6];
    grad_op(ielem, 1, 7) = g[7];

    grad(y, x, g);

    grad_op(ielem, 2, 0) = g[0];
    grad_op(ielem, 2, 1) = g[1];
    grad_op(ielem, 2, 2) = g[2];
    grad_op(ielem, 2, 3) = g[3];
    grad_op(ielem, 2, 4) = g[4];
    grad_op(ielem, 2, 5) = g[5];
    grad_op(ielem, 2, 6) = g[6];
    grad_op(ielem, 2, 7) = g[7];
  }

  //--------------------------------------------------------------------------

  struct Init {
    using execution_space = typename self_type::execution_space;

    elem_coord_type coords;

    Init(const elem_coord_type& arg_coords) : coords(arg_coords) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(size_type ielem) const {
      coords(ielem, 0, 0) = 0.;
      coords(ielem, 1, 0) = 0.;
      coords(ielem, 2, 0) = 0.;

      coords(ielem, 0, 1) = 1.;
      coords(ielem, 1, 1) = 0.;
      coords(ielem, 2, 1) = 0.;

      coords(ielem, 0, 2) = 1.;
      coords(ielem, 1, 2) = 1.;
      coords(ielem, 2, 2) = 0.;

      coords(ielem, 0, 3) = 0.;
      coords(ielem, 1, 3) = 1.;
      coords(ielem, 2, 3) = 0.;

      coords(ielem, 0, 4) = 0.;
      coords(ielem, 1, 4) = 0.;
      coords(ielem, 2, 4) = 1.;

      coords(ielem, 0, 5) = 1.;
      coords(ielem, 1, 5) = 0.;
      coords(ielem, 2, 5) = 1.;

      coords(ielem, 0, 6) = 1.;
      coords(ielem, 1, 6) = 1.;
      coords(ielem, 2, 6) = 1.;

      coords(ielem, 0, 7) = 0.;
      coords(ielem, 1, 7) = 1.;
      coords(ielem, 2, 7) = 1.;
    }
  };

  //--------------------------------------------------------------------------

  static double test(const int count) {
    elem_coord_type coord("coord", count);
    elem_grad_type grad("grad", count);

    // Execute the parallel kernels on the arrays:
    Kokkos::parallel_for(count, Init(coord));
    execution_space().fence();

    Kokkos::Timer timer;
    Kokkos::parallel_for(count, HexGrad<execution_space>(coord, grad));
    execution_space().fence();
    return timer.seconds();
  }
};

template <class CoordScalarType>
static void HexGrad_Benchmark(benchmark::State& state) {
  const auto parallel_work_length = state.range(0);

  for (auto _ : state) {
    const auto time =
        HexGrad<Kokkos::DefaultExecutionSpace, CoordScalarType>::test(
            parallel_work_length);

    state.SetIterationTime(time);
    state.counters["Count"] = benchmark::Counter(parallel_work_length);
    state.counters["Time normalized"] =
        benchmark::Counter(time / parallel_work_length);
  }
}

BENCHMARK(HexGrad_Benchmark<double>)
    ->ArgName("count")
    ->ArgsProduct({
        benchmark::CreateRange(1 << 10, 1 << 19, 2),
    })
    ->UseManualTime()
    ->Iterations(5);

}  // namespace Test
