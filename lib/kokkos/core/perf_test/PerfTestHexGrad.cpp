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
#include <gtest/gtest.h>
#include <PerfTest_Category.hpp>

namespace Test {

template <class DeviceType, typename CoordScalarType = double,
          typename GradScalarType = float>
struct HexGrad {
  typedef DeviceType execution_space;
  typedef typename execution_space::size_type size_type;

  typedef HexGrad<DeviceType, CoordScalarType, GradScalarType> self_type;

  // 3D array : ( ParallelWork , Space , Node )

  enum { NSpace = 3, NNode = 8 };

  typedef Kokkos::View<CoordScalarType * [NSpace][NNode], execution_space>
      elem_coord_type;

  typedef Kokkos::View<GradScalarType * [NSpace][NNode], execution_space>
      elem_grad_type;

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
    typedef typename self_type::execution_space execution_space;

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

  static double test(const int count, const int iter = 1) {
    elem_coord_type coord("coord", count);
    elem_grad_type grad("grad", count);

    // Execute the parallel kernels on the arrays:

    double dt_min = 0;

    Kokkos::parallel_for(count, Init(coord));
    execution_space().fence();

    for (int i = 0; i < iter; ++i) {
      Kokkos::Timer timer;
      Kokkos::parallel_for(count, HexGrad<execution_space>(coord, grad));
      execution_space().fence();
      const double dt = timer.seconds();
      if (0 == i)
        dt_min = dt;
      else
        dt_min = dt < dt_min ? dt : dt_min;
    }

    return dt_min;
  }
};

template <class DeviceType>
void run_test_hexgrad(int exp_beg, int exp_end, int num_trials,
                      const char deviceTypeName[]) {
  std::string label_hexgrad;
  label_hexgrad.append("\"HexGrad< double , ");
  label_hexgrad.append(deviceTypeName);
  label_hexgrad.append(" >\"");

  for (int i = exp_beg; i < exp_end; ++i) {
    double min_seconds = 0.0;
    double max_seconds = 0.0;
    double avg_seconds = 0.0;

    const int parallel_work_length = 1 << i;

    for (int j = 0; j < num_trials; ++j) {
      const double seconds = HexGrad<DeviceType>::test(parallel_work_length);

      if (0 == j) {
        min_seconds = seconds;
        max_seconds = seconds;
      } else {
        if (seconds < min_seconds) min_seconds = seconds;
        if (seconds > max_seconds) max_seconds = seconds;
      }
      avg_seconds += seconds;
    }
    avg_seconds /= num_trials;

    std::cout << label_hexgrad << " , " << parallel_work_length << " , "
              << min_seconds << " , " << (min_seconds / parallel_work_length)
              << std::endl;
  }
}

TEST(default_exec, hexgrad) {
  int exp_beg    = 10;
  int exp_end    = 20;
  int num_trials = 5;

  if (command_line_num_args() > 1) exp_beg = atoi(command_line_arg(1));
  if (command_line_num_args() > 2) exp_end = atoi(command_line_arg(2));
  if (command_line_num_args() > 3) num_trials = atoi(command_line_arg(3));

  EXPECT_NO_THROW(run_test_hexgrad<Kokkos::DefaultExecutionSpace>(
      exp_beg, exp_end, num_trials, Kokkos::DefaultExecutionSpace::name()));
}

}  // namespace Test
