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

// The function "test_view_collection" exposes the copy constructor
// and destructor overheads in Kokkos View objects
// Please see the lines marked by "NOTE".

#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>
#include <Kokkos_Core.hpp>
#include <iostream>

// NVIEWS is the number of Kokkos View objects in our ViewCollection object
// We have chosen a large value of 40 to make it easier to see performance
// differences when using the likelihood attribute
#define NVIEWS 40

class ViewCollection {
 public:
  Kokkos::View<double*> v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13,
      v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28,
      v29, v30, v31, v32, v33, v34, v35, v36, v37, v38, v39, v40;
  double m_expected_sum;
  double m_side_effect;
  int m_N;

  ViewCollection(int N)
      : v1("v1", N),
        v2("v2", N),
        v3("v3", N),
        v4("v4", N),
        v5("v5", N),
        v6("v6", N),
        v7("v7", N),
        v8("v8", N),
        v9("v9", N),
        v10("v10", N),
        v11("v11", N),
        v12("v12", N),
        v13("v13", N),
        v14("v14", N),
        v15("v15", N),
        v16("v16", N),
        v17("v17", N),
        v18("v18", N),
        v19("v19", N),
        v20("v20", N),
        v21("v21", N),
        v22("v22", N),
        v23("v23", N),
        v24("v24", N),
        v25("v25", N),
        v26("v26", N),
        v27("v27", N),
        v28("v28", N),
        v29("v29", N),
        v30("v30", N),
        v31("v31", N),
        v32("v32", N),
        v33("v33", N),
        v34("v34", N),
        v35("v35", N),
        v36("v36", N),
        v37("v37", N),
        v38("v38", N),
        v39("v39", N),
        v40("v40", N),
        m_expected_sum(N * NVIEWS),
        m_side_effect(0.0),
        m_N(N) {
    for (int i = 0; i < N; ++i) {
      v1(i)  = 1;
      v2(i)  = 1;
      v3(i)  = 1;
      v4(i)  = 1;
      v5(i)  = 1;
      v6(i)  = 1;
      v7(i)  = 1;
      v8(i)  = 1;
      v9(i)  = 1;
      v10(i) = 1;
      v11(i) = 1;
      v12(i) = 1;
      v13(i) = 1;
      v14(i) = 1;
      v15(i) = 1;
      v16(i) = 1;
      v17(i) = 1;
      v18(i) = 1;
      v19(i) = 1;
      v20(i) = 1;
      v21(i) = 1;
      v22(i) = 1;
      v23(i) = 1;
      v24(i) = 1;
      v25(i) = 1;
      v26(i) = 1;
      v27(i) = 1;
      v28(i) = 1;
      v29(i) = 1;
      v30(i) = 1;
      v31(i) = 1;
      v32(i) = 1;
      v33(i) = 1;
      v34(i) = 1;
      v35(i) = 1;
      v36(i) = 1;
      v37(i) = 1;
      v38(i) = 1;
      v39(i) = 1;
      v40(i) = 1;
    }
  }

// The ADD_COPY_CONSTRUCTOR macro is helpful to compare time in the copy
// constructor between compilers. We have found that the GNU compiler
// is sometimes able to inline the default copy constructor.
#ifdef ADD_COPY_CONSTRUCTOR
  __attribute__((noinline)) ViewCollection(const ViewCollection& other)
      : v1(other.v1),
        v2(other.v2),
        v3(other.v3),
        v4(other.v4),
        v5(other.v5),
        v6(other.v6),
        v7(other.v7),
        v8(other.v8),
        v9(other.v9),
        v10(other.v10),
        v11(other.v11),
        v12(other.v12),
        v13(other.v13),
        v14(other.v14),
        v15(other.v15),
        v16(other.v16),
        v17(other.v17),
        v18(other.v18),
        v19(other.v19),
        v20(other.v20),
        v21(other.v21),
        v22(other.v22),
        v23(other.v23),
        v24(other.v24),
        v25(other.v25),
        v26(other.v26),
        v27(other.v27),
        v28(other.v28),
        v29(other.v29),
        v30(other.v30),
        v31(other.v31),
        v32(other.v32),
        v33(other.v33),
        v34(other.v34),
        v35(other.v35),
        v36(other.v36),
        v37(other.v37),
        v38(other.v38),
        v39(other.v39),
        v40(other.v40),
        m_expected_sum(other.m_expected_sum),
        m_side_effect(other.m_side_effect),
        m_N(other.m_N) {}
#endif

  KOKKOS_INLINE_FUNCTION
  double sum_views(int ii, bool execute_kernel) {
    double result = 0.0;
    if (execute_kernel) {
      // This code is only executed when using the command line option -k
      // The computation references all Kokkos views. This may help our
      // effort to stop compilers from optimizing away the Kokkos views
      for (int i = 0; i < m_N; ++i) {
        result += v1(i) + v2(i) + v3(i) + v4(i) + v5(i) + v6(i) + v7(i) +
                  v8(i) + v9(i) + v10(i) + v11(i) + v12(i) + v13(i) + v14(i) +
                  v15(i) + v16(i) + v17(i) + v18(i) + v19(i) + v20(i) + v21(i) +
                  v22(i) + v23(i) + v24(i) + v25(i) + v26(i) + v27(i) + v28(i) +
                  v29(i) + v30(i) + v31(i) + v32(i) + v33(i) + v34(i) + v35(i) +
                  v36(i) + v37(i) + v38(i) + v39(i) + v40(i);
      }
    } else {
      result = m_expected_sum;
    }
    // This statement introduces a side effect that may help our effort to
    // stop compilers from optimizing away the temporary ViewCollection object
    m_side_effect = result * (ii + 1);
    return result;
  }
};

void test_view_collection_kk(int N, int num_iter, bool execute_kernel) {
  ViewCollection view_collection(N);

  Kokkos::Timer view_collection_timer;
  double max_value = 0.0;
  // Max Reduction boilerplate code taken from slide 53 of
  // kokkos-tutorials/LectureSeries/KokkosTutorial_02_ViewsAndSpaces.pdf
  Kokkos::parallel_reduce(
      "collection-reduction", num_iter,
      KOKKOS_LAMBDA(int i, double& valueToUpdate) {
        // NOTE: The following lines expose the Kokkos View overheads
        ViewCollection tmp_view_collection = view_collection;
        double my_value = tmp_view_collection.sum_views(i, execute_kernel);
        if (my_value > valueToUpdate) valueToUpdate = my_value;
      },
      Kokkos::Max<double>(max_value));
  double view_collection_time = view_collection_timer.seconds();

  bool success = std::fabs(max_value - N * NVIEWS) < 1.E-6;
  std::cout << "View Time = " << view_collection_time << " seconds"
            << std::endl;
  if (success) {
    std::cout << "Kokkos run:" << std::endl;
    std::cout << "SUCCESS" << std::endl;
  } else {
    std::cout << "FAILURE" << std::endl;
  }
}

void test_view_collection_serial(int N, int num_iter, bool execute_kernel) {
  ViewCollection view_collection(N);

  Kokkos::Timer view_collection_timer;
  double max_value = 0.0;
  // Max Reduction boilerplate code taken from slide 53 of
  // kokkos-tutorials/LectureSeries/KokkosTutorial_02_ViewsAndSpaces.pdf
  for (int i = 0; i < num_iter; ++i) {
    // NOTE: The following lines expose the Kokkos View overheads
    ViewCollection tmp_view_collection = view_collection;
    double my_value = tmp_view_collection.sum_views(i, execute_kernel);
    if (my_value > max_value) max_value = my_value;
  }
  double view_collection_time = view_collection_timer.seconds();

  bool success = std::fabs(max_value - N * NVIEWS) < 1.E-6;
  std::cout << "View Time 2 = " << view_collection_time << " seconds"
            << std::endl;
  if (success) {
    std::cout << "Serial run:" << std::endl;
    std::cout << "SUCCESS" << std::endl;
  } else {
    std::cout << "FAILURE" << std::endl;
  }
}

int main(int argc, char* argv[]) {
  // The benchmark is only testing reference counting for views on host.
#if defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_SERIAL) || \
    defined(KOKKOS_ENABLE_THREADS) || defined(KOKKOS_ENABLE_HPX)
  int N               = 1;
  int num_iter        = 1 << 27;
  bool execute_kernel = false;

  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-N") == 0)) {
      N = atoi(argv[++i]);
      if (N < 1) {
        std::cout << "Array extent must be >= 1" << std::endl;
        exit(1);
      }
    } else if (strcmp(argv[i], "-i") == 0) {
      num_iter = atoi(argv[++i]);
      if (num_iter < 1) {
        std::cout << "Number of iterations must be >= 1" << std::endl;
        exit(1);
      }
    } else if (strcmp(argv[i], "-k") == 0) {
      execute_kernel = true;
    } else if ((strcmp(argv[i], "-h") == 0)) {
      printf("  Options:\n");
      printf("  -N <int>: Array extent\n");
      printf("  -i <int>: Number of iterations\n");
      printf("  -k:       Execute the summation kernel\n");
      printf("  -h:       Print this message\n\n");
      exit(1);
    }
  }

  std::cout << "Array extent = " << N << std::endl;
  std::cout << "Iterations = " << num_iter << std::endl;
  std::cout << "Execute summation kernel = " << std::boolalpha << execute_kernel
            << std::noboolalpha << std::endl;

  // Test inside a Kokkos kernel.
  Kokkos::initialize(argc, argv);
  { test_view_collection_kk(N, num_iter, execute_kernel); }

  // Test outside Kokkos kernel.
  test_view_collection_serial(N, num_iter, execute_kernel);

  Kokkos::finalize();
#endif

  return 0;
}
