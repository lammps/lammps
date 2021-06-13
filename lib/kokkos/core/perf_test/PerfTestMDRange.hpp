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

namespace Test {
template <class DeviceType, typename ScalarType = double,
          typename TestLayout = Kokkos::LayoutRight>
struct MultiDimRangePerf3D {
  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;

  using iterate_type = Kokkos::Iterate;

  using view_type      = Kokkos::View<ScalarType ***, TestLayout, DeviceType>;
  using host_view_type = typename view_type::HostMirror;

  view_type A;
  view_type B;
  const long irange;
  const long jrange;
  const long krange;

  MultiDimRangePerf3D(const view_type &A_, const view_type &B_,
                      const long &irange_, const long &jrange_,
                      const long &krange_)
      : A(A_), B(B_), irange(irange_), jrange(jrange_), krange(krange_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const long i, const long j, const long k) const {
    A(i, j, k) =
        0.25 * (ScalarType)(B(i + 2, j, k) + B(i + 1, j, k) + B(i, j + 2, k) +
                            B(i, j + 1, k) + B(i, j, k + 2) + B(i, j, k + 1) +
                            B(i, j, k));
  }

  struct InitZeroTag {};
  //  struct InitViewTag {};

  struct Init {
    Init(const view_type &input_, const long &irange_, const long &jrange_,
         const long &krange_)
        : input(input_), irange(irange_), jrange(jrange_), krange(krange_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const long i, const long j, const long k) const {
      input(i, j, k) = 1.0;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const InitZeroTag &, const long i, const long j,
                    const long k) const {
      input(i, j, k) = 0;
    }

    view_type input;
    const long irange;
    const long jrange;
    const long krange;
  };

  static double test_multi_index(const unsigned int icount,
                                 const unsigned int jcount,
                                 const unsigned int kcount,
                                 const unsigned int Ti = 1,
                                 const unsigned int Tj = 1,
                                 const unsigned int Tk = 1,
                                 const long iter       = 1) {
    // This test performs multidim range over all dims
    view_type Atest("Atest", icount, jcount, kcount);
    view_type Btest("Btest", icount + 2, jcount + 2, kcount + 2);
    using FunctorType =
        MultiDimRangePerf3D<execution_space, ScalarType, TestLayout>;

    double dt_min = 0;

    // LayoutRight
    if (std::is_same<TestLayout, Kokkos::LayoutRight>::value) {
      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, iterate_type::Right, iterate_type::Right>,
          execution_space>
          policy_initA({{0, 0, 0}}, {{icount, jcount, kcount}}, {{Ti, Tj, Tk}});
      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, iterate_type::Right, iterate_type::Right>,
          execution_space>
          policy_initB({{0, 0, 0}}, {{icount + 2, jcount + 2, kcount + 2}},
                       {{Ti, Tj, Tk}});

      using MDRangeType = typename Kokkos::MDRangePolicy<
          Kokkos::Rank<3, iterate_type::Right, iterate_type::Right>,
          execution_space>;
      using tile_type  = typename MDRangeType::tile_type;
      using point_type = typename MDRangeType::point_type;

      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, iterate_type::Right, iterate_type::Right>,
          execution_space>
          policy(point_type{{0, 0, 0}}, point_type{{icount, jcount, kcount}},
                 tile_type{{Ti, Tj, Tk}});

      Kokkos::parallel_for(policy_initA, Init(Atest, icount, jcount, kcount));
      execution_space().fence();
      Kokkos::parallel_for(policy_initB,
                           Init(Btest, icount + 2, jcount + 2, kcount + 2));
      execution_space().fence();

      for (int i = 0; i < iter; ++i) {
        Kokkos::Timer timer;
        Kokkos::parallel_for(policy,
                             FunctorType(Atest, Btest, icount, jcount, kcount));
        execution_space().fence();
        const double dt = timer.seconds();
        if (0 == i)
          dt_min = dt;
        else
          dt_min = dt < dt_min ? dt : dt_min;

        // Correctness check - only the first run
        if (0 == i) {
          long numErrors = 0;
          host_view_type Ahost("Ahost", icount, jcount, kcount);
          Kokkos::deep_copy(Ahost, Atest);
          host_view_type Bhost("Bhost", icount + 2, jcount + 2, kcount + 2);
          Kokkos::deep_copy(Bhost, Btest);

          // On KNL, this may vectorize - add print statement to prevent
          // Also, compare against epsilon, as vectorization can change bitwise
          // answer
          for (long l = 0; l < static_cast<long>(icount); ++l) {
            for (long j = 0; j < static_cast<long>(jcount); ++j) {
              for (long k = 0; k < static_cast<long>(kcount); ++k) {
                ScalarType check =
                    0.25 *
                    (ScalarType)(Bhost(l + 2, j, k) + Bhost(l + 1, j, k) +
                                 Bhost(l, j + 2, k) + Bhost(l, j + 1, k) +
                                 Bhost(l, j, k + 2) + Bhost(l, j, k + 1) +
                                 Bhost(l, j, k));
                if (Ahost(l, j, k) - check != 0) {
                  ++numErrors;
                  std::cout << "  Correctness error at index: " << l << "," << j
                            << "," << k << "\n"
                            << "  multi Ahost = " << Ahost(l, j, k)
                            << "  expected = " << check
                            << "  multi Bhost(ijk) = " << Bhost(l, j, k)
                            << "  multi Bhost(l+1jk) = " << Bhost(l + 1, j, k)
                            << "  multi Bhost(l+2jk) = " << Bhost(l + 2, j, k)
                            << "  multi Bhost(ij+1k) = " << Bhost(l, j + 1, k)
                            << "  multi Bhost(ij+2k) = " << Bhost(l, j + 2, k)
                            << "  multi Bhost(ijk+1) = " << Bhost(l, j, k + 1)
                            << "  multi Bhost(ijk+2) = " << Bhost(l, j, k + 2)
                            << std::endl;
                  // exit(-1);
                }
              }
            }
          }
          if (numErrors != 0) {
            std::cout << "LR multi: errors " << numErrors << "  range product "
                      << icount * jcount * kcount << "  LL " << jcount * kcount
                      << "  LR " << icount * jcount << std::endl;
          }
          // else { std::cout << " multi: No errors!" <<  std::endl; }
        }
      }  // end for

    }
    // LayoutLeft
    else {
      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, iterate_type::Left, iterate_type::Left>,
          execution_space>
          policy_initA({{0, 0, 0}}, {{icount, jcount, kcount}}, {{Ti, Tj, Tk}});
      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, iterate_type::Left, iterate_type::Left>,
          execution_space>
          policy_initB({{0, 0, 0}}, {{icount + 2, jcount + 2, kcount + 2}},
                       {{Ti, Tj, Tk}});

      // using MDRangeType =
      //     typename Kokkos::MDRangePolicy<
      //         Kokkos::Rank<3, iterate_type::Left, iterate_type::Left>,
      //         execution_space >;
      // using tile_type = typename MDRangeType::tile_type;
      // using point_type = typename MDRangeType::point_type;
      // MDRangeType policy(point_type{{0,0,0}},
      //                    point_type{{icount,jcount,kcount}},
      //                    tile_type{{Ti,Tj,Tk}});
      Kokkos::MDRangePolicy<
          Kokkos::Rank<3, iterate_type::Left, iterate_type::Left>,
          execution_space>
          policy({{0, 0, 0}}, {{icount, jcount, kcount}}, {{Ti, Tj, Tk}});

      Kokkos::parallel_for(policy_initA, Init(Atest, icount, jcount, kcount));
      execution_space().fence();
      Kokkos::parallel_for(policy_initB,
                           Init(Btest, icount + 2, jcount + 2, kcount + 2));
      execution_space().fence();

      for (int i = 0; i < iter; ++i) {
        Kokkos::Timer timer;
        Kokkos::parallel_for(policy,
                             FunctorType(Atest, Btest, icount, jcount, kcount));
        execution_space().fence();
        const double dt = timer.seconds();
        if (0 == i)
          dt_min = dt;
        else
          dt_min = dt < dt_min ? dt : dt_min;

        // Correctness check - only the first run
        if (0 == i) {
          long numErrors = 0;
          host_view_type Ahost("Ahost", icount, jcount, kcount);
          Kokkos::deep_copy(Ahost, Atest);
          host_view_type Bhost("Bhost", icount + 2, jcount + 2, kcount + 2);
          Kokkos::deep_copy(Bhost, Btest);

          // On KNL, this may vectorize - add print statement to prevent
          // Also, compare against epsilon, as vectorization can change bitwise
          // answer
          for (long l = 0; l < static_cast<long>(icount); ++l) {
            for (long j = 0; j < static_cast<long>(jcount); ++j) {
              for (long k = 0; k < static_cast<long>(kcount); ++k) {
                ScalarType check =
                    0.25 *
                    (ScalarType)(Bhost(l + 2, j, k) + Bhost(l + 1, j, k) +
                                 Bhost(l, j + 2, k) + Bhost(l, j + 1, k) +
                                 Bhost(l, j, k + 2) + Bhost(l, j, k + 1) +
                                 Bhost(l, j, k));
                if (Ahost(l, j, k) - check != 0) {
                  ++numErrors;
                  std::cout << "  Correctness error at index: " << l << "," << j
                            << "," << k << "\n"
                            << "  multi Ahost = " << Ahost(l, j, k)
                            << "  expected = " << check
                            << "  multi Bhost(ijk) = " << Bhost(l, j, k)
                            << "  multi Bhost(l+1jk) = " << Bhost(l + 1, j, k)
                            << "  multi Bhost(l+2jk) = " << Bhost(l + 2, j, k)
                            << "  multi Bhost(ij+1k) = " << Bhost(l, j + 1, k)
                            << "  multi Bhost(ij+2k) = " << Bhost(l, j + 2, k)
                            << "  multi Bhost(ijk+1) = " << Bhost(l, j, k + 1)
                            << "  multi Bhost(ijk+2) = " << Bhost(l, j, k + 2)
                            << std::endl;
                  // exit(-1);
                }
              }
            }
          }
          if (numErrors != 0) {
            std::cout << " LL multi run: errors " << numErrors
                      << "  range product " << icount * jcount * kcount
                      << "  LL " << jcount * kcount << "  LR "
                      << icount * jcount << std::endl;
          }
          // else { std::cout << " multi: No errors!" <<  std::endl; }
        }
      }  // end for
    }

    return dt_min;
  }
};

template <class DeviceType, typename ScalarType = double,
          typename TestLayout = Kokkos::LayoutRight>
struct RangePolicyCollapseTwo {
  // RangePolicy for 3D range, but will collapse only 2 dims => like Rank<2> for
  // multi-dim; unroll 2 dims in one-dim

  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;
  using layout          = TestLayout;

  using iterate_type = Kokkos::Iterate;

  using view_type      = Kokkos::View<ScalarType ***, TestLayout, DeviceType>;
  using host_view_type = typename view_type::HostMirror;

  view_type A;
  view_type B;
  const long irange;
  const long jrange;
  const long krange;

  RangePolicyCollapseTwo(view_type &A_, const view_type &B_,
                         const long &irange_, const long &jrange_,
                         const long &krange_)
      : A(A_), B(B_), irange(irange_), jrange(jrange_), krange(krange_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const long r) const {
    if (std::is_same<TestLayout, Kokkos::LayoutRight>::value) {
      // id(i,j,k) = k + j*Nk + i*Nk*Nj = k + Nk*(j + i*Nj) = k + Nk*r
      // r = j + i*Nj
      long i = int(r / jrange);
      long j = int(r - i * jrange);
      for (int k = 0; k < krange; ++k) {
        A(i, j, k) =
            0.25 * (ScalarType)(B(i + 2, j, k) + B(i + 1, j, k) +
                                B(i, j + 2, k) + B(i, j + 1, k) +
                                B(i, j, k + 2) + B(i, j, k + 1) + B(i, j, k));
      }
    } else if (std::is_same<TestLayout, Kokkos::LayoutLeft>::value) {
      // id(i,j,k) = i + j*Ni + k*Ni*Nj = i + Ni*(j + k*Nj) = i + Ni*r
      // r = j + k*Nj
      long k = int(r / jrange);
      long j = int(r - k * jrange);
      for (int i = 0; i < irange; ++i) {
        A(i, j, k) =
            0.25 * (ScalarType)(B(i + 2, j, k) + B(i + 1, j, k) +
                                B(i, j + 2, k) + B(i, j + 1, k) +
                                B(i, j, k + 2) + B(i, j, k + 1) + B(i, j, k));
      }
    }
  }

  struct Init {
    view_type input;
    const long irange;
    const long jrange;
    const long krange;

    Init(const view_type &input_, const long &irange_, const long &jrange_,
         const long &krange_)
        : input(input_), irange(irange_), jrange(jrange_), krange(krange_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const long r) const {
      if (std::is_same<TestLayout, Kokkos::LayoutRight>::value) {
        long i = int(r / jrange);
        long j = int(r - i * jrange);
        for (int k = 0; k < krange; ++k) {
          input(i, j, k) = 1;
        }
      } else if (std::is_same<TestLayout, Kokkos::LayoutLeft>::value) {
        long k = int(r / jrange);
        long j = int(r - k * jrange);
        for (int i = 0; i < irange; ++i) {
          input(i, j, k) = 1;
        }
      }
    }
  };

  static double test_index_collapse_two(const unsigned int icount,
                                        const unsigned int jcount,
                                        const unsigned int kcount,
                                        const long iter = 1) {
    // This test refers to collapsing two dims while using the RangePolicy
    view_type Atest("Atest", icount, jcount, kcount);
    view_type Btest("Btest", icount + 2, jcount + 2, kcount + 2);
    using FunctorType =
        RangePolicyCollapseTwo<execution_space, ScalarType, TestLayout>;

    long collapse_index_rangeA = 0;
    long collapse_index_rangeB = 0;
    if (std::is_same<TestLayout, Kokkos::LayoutRight>::value) {
      collapse_index_rangeA = icount * jcount;
      collapse_index_rangeB = (icount + 2) * (jcount + 2);
      //      std::cout << "   LayoutRight " << std::endl;
    } else if (std::is_same<TestLayout, Kokkos::LayoutLeft>::value) {
      collapse_index_rangeA = kcount * jcount;
      collapse_index_rangeB = (kcount + 2) * (jcount + 2);
      //      std::cout << "   LayoutLeft " << std::endl;
    } else {
      std::cout << "  LayoutRight or LayoutLeft required - will pass 0 as "
                   "range instead "
                << std::endl;
      exit(-1);
    }

    Kokkos::RangePolicy<execution_space> policy(0, (collapse_index_rangeA));
    Kokkos::RangePolicy<execution_space> policy_initB(0,
                                                      (collapse_index_rangeB));

    double dt_min = 0;

    Kokkos::parallel_for(policy, Init(Atest, icount, jcount, kcount));
    execution_space().fence();
    Kokkos::parallel_for(policy_initB,
                         Init(Btest, icount + 2, jcount + 2, kcount + 2));
    execution_space().fence();

    for (int i = 0; i < iter; ++i) {
      Kokkos::Timer timer;
      Kokkos::parallel_for(policy,
                           FunctorType(Atest, Btest, icount, jcount, kcount));
      execution_space().fence();
      const double dt = timer.seconds();
      if (0 == i)
        dt_min = dt;
      else
        dt_min = dt < dt_min ? dt : dt_min;

      // Correctness check - first iteration only
      if (0 == i) {
        long numErrors = 0;
        host_view_type Ahost("Ahost", icount, jcount, kcount);
        Kokkos::deep_copy(Ahost, Atest);
        host_view_type Bhost("Bhost", icount + 2, jcount + 2, kcount + 2);
        Kokkos::deep_copy(Bhost, Btest);

        // On KNL, this may vectorize - add print statement to prevent
        // Also, compare against epsilon, as vectorization can change bitwise
        // answer
        for (long l = 0; l < static_cast<long>(icount); ++l) {
          for (long j = 0; j < static_cast<long>(jcount); ++j) {
            for (long k = 0; k < static_cast<long>(kcount); ++k) {
              ScalarType check =
                  0.25 * (ScalarType)(Bhost(l + 2, j, k) + Bhost(l + 1, j, k) +
                                      Bhost(l, j + 2, k) + Bhost(l, j + 1, k) +
                                      Bhost(l, j, k + 2) + Bhost(l, j, k + 1) +
                                      Bhost(l, j, k));
              if (Ahost(l, j, k) - check != 0) {
                ++numErrors;
                std::cout << "  Correctness error at index: " << l << "," << j
                          << "," << k << "\n"
                          << "  flat Ahost = " << Ahost(l, j, k)
                          << "  expected = " << check << std::endl;
                // exit(-1);
              }
            }
          }
        }
        if (numErrors != 0) {
          std::cout << " RP collapse2: errors " << numErrors
                    << "  range product " << icount * jcount * kcount << "  LL "
                    << jcount * kcount << "  LR " << icount * jcount
                    << std::endl;
        }
        // else { std::cout << " RP collapse2: Pass! " << std::endl; }
      }
    }

    return dt_min;
  }
};

template <class DeviceType, typename ScalarType = double,
          typename TestLayout = Kokkos::LayoutRight>
struct RangePolicyCollapseAll {
  // RangePolicy for 3D range, but will collapse all dims

  using execution_space = DeviceType;
  using size_type       = typename execution_space::size_type;
  using layout          = TestLayout;

  using view_type      = Kokkos::View<ScalarType ***, TestLayout, DeviceType>;
  using host_view_type = typename view_type::HostMirror;

  view_type A;
  view_type B;
  const long irange;
  const long jrange;
  const long krange;

  RangePolicyCollapseAll(view_type &A_, const view_type &B_,
                         const long &irange_, const long &jrange_,
                         const long &krange_)
      : A(A_), B(B_), irange(irange_), jrange(jrange_), krange(krange_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const long r) const {
    if (std::is_same<TestLayout, Kokkos::LayoutRight>::value) {
      long i = int(r / (jrange * krange));
      long j = int((r - i * jrange * krange) / krange);
      long k = int(r - i * jrange * krange - j * krange);
      A(i, j, k) =
          0.25 * (ScalarType)(B(i + 2, j, k) + B(i + 1, j, k) + B(i, j + 2, k) +
                              B(i, j + 1, k) + B(i, j, k + 2) + B(i, j, k + 1) +
                              B(i, j, k));
    } else if (std::is_same<TestLayout, Kokkos::LayoutLeft>::value) {
      long k = int(r / (irange * jrange));
      long j = int((r - k * irange * jrange) / irange);
      long i = int(r - k * irange * jrange - j * irange);
      A(i, j, k) =
          0.25 * (ScalarType)(B(i + 2, j, k) + B(i + 1, j, k) + B(i, j + 2, k) +
                              B(i, j + 1, k) + B(i, j, k + 2) + B(i, j, k + 1) +
                              B(i, j, k));
    }
  }

  struct Init {
    view_type input;
    const long irange;
    const long jrange;
    const long krange;

    Init(const view_type &input_, const long &irange_, const long &jrange_,
         const long &krange_)
        : input(input_), irange(irange_), jrange(jrange_), krange(krange_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const long r) const {
      if (std::is_same<TestLayout, Kokkos::LayoutRight>::value) {
        long i         = int(r / (jrange * krange));
        long j         = int((r - i * jrange * krange) / krange);
        long k         = int(r - i * jrange * krange - j * krange);
        input(i, j, k) = 1;
      } else if (std::is_same<TestLayout, Kokkos::LayoutLeft>::value) {
        long k         = int(r / (irange * jrange));
        long j         = int((r - k * irange * jrange) / irange);
        long i         = int(r - k * irange * jrange - j * irange);
        input(i, j, k) = 1;
      }
    }
  };

  static double test_collapse_all(const unsigned int icount,
                                  const unsigned int jcount,
                                  const unsigned int kcount,
                                  const long iter = 1) {
    // This test refers to collapsing all dims using the RangePolicy
    view_type Atest("Atest", icount, jcount, kcount);
    view_type Btest("Btest", icount + 2, jcount + 2, kcount + 2);
    using FunctorType =
        RangePolicyCollapseAll<execution_space, ScalarType, TestLayout>;

    const long flat_index_range = icount * jcount * kcount;
    Kokkos::RangePolicy<execution_space> policy(0, flat_index_range);
    Kokkos::RangePolicy<execution_space> policy_initB(
        0, (icount + 2) * (jcount + 2) * (kcount + 2));

    double dt_min = 0;

    Kokkos::parallel_for(policy, Init(Atest, icount, jcount, kcount));
    execution_space().fence();
    Kokkos::parallel_for(policy_initB,
                         Init(Btest, icount + 2, jcount + 2, kcount + 2));
    execution_space().fence();

    for (int i = 0; i < iter; ++i) {
      Kokkos::Timer timer;
      Kokkos::parallel_for(policy,
                           FunctorType(Atest, Btest, icount, jcount, kcount));
      execution_space().fence();
      const double dt = timer.seconds();
      if (0 == i)
        dt_min = dt;
      else
        dt_min = dt < dt_min ? dt : dt_min;

      // Correctness check - first iteration only
      if (0 == i) {
        long numErrors = 0;
        host_view_type Ahost("Ahost", icount, jcount, kcount);
        Kokkos::deep_copy(Ahost, Atest);
        host_view_type Bhost("Bhost", icount + 2, jcount + 2, kcount + 2);
        Kokkos::deep_copy(Bhost, Btest);

        // On KNL, this may vectorize - add print statement to prevent
        // Also, compare against epsilon, as vectorization can change bitwise
        // answer
        for (long l = 0; l < static_cast<long>(icount); ++l) {
          for (long j = 0; j < static_cast<long>(jcount); ++j) {
            for (long k = 0; k < static_cast<long>(kcount); ++k) {
              ScalarType check =
                  0.25 * (ScalarType)(Bhost(l + 2, j, k) + Bhost(l + 1, j, k) +
                                      Bhost(l, j + 2, k) + Bhost(l, j + 1, k) +
                                      Bhost(l, j, k + 2) + Bhost(l, j, k + 1) +
                                      Bhost(l, j, k));
              if (Ahost(l, j, k) - check != 0) {
                ++numErrors;
                std::cout << "  Callapse ALL Correctness error at index: " << l
                          << "," << j << "," << k << "\n"
                          << "  flat Ahost = " << Ahost(l, j, k)
                          << "  expected = " << check << std::endl;
                // exit(-1);
              }
            }
          }
        }
        if (numErrors != 0) {
          std::cout << " RP collapse all: errors " << numErrors
                    << "  range product " << icount * jcount * kcount << "  LL "
                    << jcount * kcount << "  LR " << icount * jcount
                    << std::endl;
        }
        // else { std::cout << " RP collapse all: Pass! " << std::endl; }
      }
    }

    return dt_min;
  }
};

}  // end namespace Test
