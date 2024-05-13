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
#include <cstdio>

//
// MDRangePolicy example with parallel_for and parallel_reduce:
//   1. Start up Kokkos
//   2. Execute a parallel_for loop in the default execution space,
//      using a functor to define the loop body
//   3. Shut down Kokkos
//
// Two examples are provided:
// Example 1: Rank 2 case with minimal default parameters and arguments used
//            in the MDRangePolicy
//
// Example 2: Rank 3 case with additional outer/inner iterate pattern parameters
//            and tile dims passed to the ctor

// Simple functor for computing/storing the product of indices in a View v
template <class ViewType>
struct MDFunctor2D {
  using value_type = long;

  ViewType v;
  size_t size;

  MDFunctor2D(const ViewType& v_, const size_t size_) : v(v_), size(size_) {}

  // 2D case - used by parallel_for
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {
    v(i, j) = i * j;  // compute the product of indices
  }

  // 2D case - reduction
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, value_type& incorrect_count) const {
    if (v(i, j) != i * j) {
      incorrect_count += 1;
    }
  }
};

template <class ViewType>
struct MDFunctor3D {
  using value_type = long;

  ViewType v;
  size_t size;

  MDFunctor3D(const ViewType& v_, const size_t size_) : v(v_), size(size_) {}

  // 3D case - used by parallel_for
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    v(i, j, k) = i * j * k;  // compute the product of indices
  }

  // 3D case - reduction
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k,
                  value_type& incorrect_count) const {
    if (v(i, j, k) != i * j * k) {
      incorrect_count += 1;
    }
  }
};

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);

  // Bound(s) for MDRangePolicy
  const int n = 100;

  // ViewType aliases for Rank<2>, Rank<3> for example usage
  using ScalarType  = double;
  using ViewType_2D = Kokkos::View<ScalarType**>;
  using ViewType_3D = Kokkos::View<ScalarType***>;

  /////////////////////////////////////////////////////////////////////////////
  // Explanation of MDRangePolicy usage, template parameters, constructor
  // arguments
  //
  // MDRangePolicy aliases for Rank<2>, Rank<3> cases
  // Required template parameters:
  //   Kokkos::Rank<N>: where N=rank
  //
  // Optional template parameters to Rank<...>:
  //   Kokkos::Iterate::{Default,Left,Right}: Outer iteration pattern across
  //   tiles;
  //     defaults based on the execution space similar to Kokkos::Layout
  //   Kokkos::Iterate::{Default,Left,Right}: Inner iteration pattern within
  //   tiles;
  //     defaults based on the execution space similar to Kokkos::Layout
  //
  //   e.g. using rank2ll = Rank<2, Iterate::Left, Iterate::Left>;
  //
  //
  // Optional template parameters to MDRangePolicy:
  //   ExecutionSpace: Kokkos::Serial, Kokkos::OpenMP, Kokkos::Cuda, etc.
  //
  //   Kokkos::IndexType< T >: where T = int, long, unsigned int, etc.
  //
  //   struct Tag{}: A user-provided tag for tagging functor operators
  //
  //   e.g. 1:  MDRangePolicy< Kokkos::Serial, Rank<2, Iterate::Left,
  //   Iterate::Left>, IndexType<int>, Tag > mdpolicy; e.g. 2:  MDRangePolicy<
  //   Kokkos::Serial, rank2ll, IndexType<int>, Tag > mdpolicy;
  //
  //
  // Required arguments to ctor:
  //   {{ l0, l1, ... }}: Lower bounds, provided as Kokkos::Array or
  //   std::initializer_list
  //   {{ u0, u1, ... }}: Upper bounds, provided as Kokkos::Array or
  //   std::initializer_list
  //
  // Optional arguments to ctor:
  //   {{ t0, t1, ... }}: Tile dimensions, provided as Kokkos::Array or
  //   std::initializer_list
  //                      defaults based on the execution space
  //
  //  e.g. mdpolicy( {{0,0}}, {{u0,u1}}, {{t0,t1}};
  //
  /////////////////////////////////////////////////////////////////////////////

  // Example 1:
  long incorrect_count_2d = 0;
  {
    // Rank<2> Case: Rank is provided, all other parameters are default
    using MDPolicyType_2D = Kokkos::MDRangePolicy<Kokkos::Rank<2> >;

    // Construct 2D MDRangePolicy: lower and upper bounds provided, tile dims
    // defaulted
    MDPolicyType_2D mdpolicy_2d({{0, 0}}, {{n, n}});

    // Construct a 2D view to store result of product of indices
    ViewType_2D v2("v2", n, n);

    // Execute parallel_for with rank 2 MDRangePolicy
    Kokkos::parallel_for("md2d", mdpolicy_2d, MDFunctor2D<ViewType_2D>(v2, n));

    // Check results with a parallel_reduce using the MDRangePolicy
    Kokkos::parallel_reduce("md2dredux", mdpolicy_2d,
                            MDFunctor2D<ViewType_2D>(v2, n),
                            incorrect_count_2d);

    printf("Rank 2 MDRangePolicy incorrect count: %ld\n",
           incorrect_count_2d);  // should be 0
  }

  // Example 2:
  long incorrect_count_3d = 0;
  {
    // Rank<3> Case: Rank, inner iterate pattern, outer iterate pattern provided
    using MDPolicyType_3D = Kokkos::MDRangePolicy<
        Kokkos::Rank<3, Kokkos::Iterate::Left, Kokkos::Iterate::Left> >;

    // Construct 3D MDRangePolicy: lower, upper bounds, tile dims provided
    MDPolicyType_3D mdpolicy_3d({{0, 0, 0}}, {{n, n, n}}, {{4, 4, 4}});

    // Construct a 3D view to store result of product of indices
    ViewType_3D v3("v3", n, n, n);

    // Execute parallel_for with rank 3 MDRangePolicy
    Kokkos::parallel_for("md3d", mdpolicy_3d, MDFunctor3D<ViewType_3D>(v3, n));

    // Check results with a parallel_reduce using the MDRangePolicy
    Kokkos::parallel_reduce("md3dredux", mdpolicy_3d,
                            MDFunctor3D<ViewType_3D>(v3, n),
                            incorrect_count_3d);

    printf("Rank 3 MDRangePolicy incorrect count: %ld\n",
           incorrect_count_3d);  // should be 0
  }

  Kokkos::finalize();

  return (incorrect_count_2d == long(0) && incorrect_count_3d == long(0)) ? 0
                                                                          : -1;
}
