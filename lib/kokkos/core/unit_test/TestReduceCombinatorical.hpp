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

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <limits>

#include <Kokkos_Core.hpp>

namespace Test {

namespace ReduceCombinatorical {

template <class Scalar, class Space = Kokkos::HostSpace>
struct AddPlus {
 public:
  // Required.
  using reducer    = AddPlus;
  using value_type = Scalar;

  using result_view_type =
      Kokkos::View<value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

 private:
  result_view_type result;

 public:
  AddPlus(value_type& result_) : result(&result_) {}

  // Required.
  KOKKOS_INLINE_FUNCTION
  void join(value_type& dest, const value_type& src) const { dest += src + 1; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile value_type& dest, const volatile value_type& src) const {
    dest += src + 1;
  }

  // Optional.
  KOKKOS_INLINE_FUNCTION
  void init(value_type& val) const { val = value_type(); }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const { return result(); }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const { return result; }
};

template <int ISTEAM>
struct FunctorScalar;

template <>
struct FunctorScalar<0> {
  Kokkos::View<double> result;

  FunctorScalar(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }
};

// FIXME_SYCL requires TeamPolicy
#ifndef KOKKOS_ENABLE_SYCL
template <>
struct FunctorScalar<1> {
  using team_type = Kokkos::TeamPolicy<>::member_type;

  Kokkos::View<double> result;

  FunctorScalar(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_type& team, double& update) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }
};
#endif

template <int ISTEAM>
struct FunctorScalarInit;

template <>
struct FunctorScalarInit<0> {
  Kokkos::View<double> result;

  FunctorScalarInit(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }

  KOKKOS_INLINE_FUNCTION
  void init(double& update) const { update = 0.0; }
};

// FIXME_SYCL requires TeamPolicy
#ifndef KOKKOS_ENABLE_SYCL
template <>
struct FunctorScalarInit<1> {
  using team_type = Kokkos::TeamPolicy<>::member_type;

  Kokkos::View<double> result;

  FunctorScalarInit(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_type& team, double& update) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void init(double& update) const { update = 0.0; }
};
#endif

template <int ISTEAM>
struct FunctorScalarFinal;

template <>
struct FunctorScalarFinal<0> {
  Kokkos::View<double> result;

  FunctorScalarFinal(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }

  KOKKOS_INLINE_FUNCTION
  void final(double& update) const { result() = update; }
};

// FIXME_SYCL requires TeamPolicy
#ifndef KOKKOS_ENABLE_SYCL
template <>
struct FunctorScalarFinal<1> {
  using team_type = Kokkos::TeamPolicy<>::member_type;

  Kokkos::View<double> result;

  FunctorScalarFinal(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_type& team, double& update) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void final(double& update) const { result() = update; }
};
#endif

template <int ISTEAM>
struct FunctorScalarJoin;

template <>
struct FunctorScalarJoin<0> {
  Kokkos::View<double> result;

  FunctorScalarJoin(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }
};

// FIXME_SYCL requires TeamPolicy
#ifndef KOKKOS_ENABLE_SYCL
template <>
struct FunctorScalarJoin<1> {
  using team_type = Kokkos::TeamPolicy<>::member_type;

  Kokkos::View<double> result;

  FunctorScalarJoin(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_type& team, double& update) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }
};
#endif

template <int ISTEAM>
struct FunctorScalarJoinFinal;

template <>
struct FunctorScalarJoinFinal<0> {
  Kokkos::View<double> result;

  FunctorScalarJoinFinal(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final(double& update) const { result() = update; }
};

// FIXME_SYCL requires TeamPolicy
#ifndef KOKKOS_ENABLE_SYCL
template <>
struct FunctorScalarJoinFinal<1> {
  using team_type = Kokkos::TeamPolicy<>::member_type;

  Kokkos::View<double> result;

  FunctorScalarJoinFinal(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_type& team, double& update) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final(double& update) const { result() = update; }
};
#endif

template <int ISTEAM>
struct FunctorScalarJoinInit;

template <>
struct FunctorScalarJoinInit<0> {
  Kokkos::View<double> result;

  FunctorScalarJoinInit(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void init(double& update) const { update = 0.0; }
};

// FIXME_SYCL requires TeamPolicy
#ifndef KOKKOS_ENABLE_SYCL
template <>
struct FunctorScalarJoinInit<1> {
  using team_type = Kokkos::TeamPolicy<>::member_type;

  Kokkos::View<double> result;

  FunctorScalarJoinInit(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_type& team, double& update) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void init(double& update) const { update = 0.0; }
};
#endif

template <int ISTEAM>
struct FunctorScalarJoinFinalInit;

template <>
struct FunctorScalarJoinFinalInit<0> {
  Kokkos::View<double> result;

  FunctorScalarJoinFinalInit(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final(double& update) const { result() = update; }

  KOKKOS_INLINE_FUNCTION
  void init(double& update) const { update = 0.0; }
};

// FIXME_SYCL requires TeamPolicy
#ifndef KOKKOS_ENABLE_SYCL
template <>
struct FunctorScalarJoinFinalInit<1> {
  using team_type = Kokkos::TeamPolicy<>::member_type;

  Kokkos::View<double> result;

  FunctorScalarJoinFinalInit(Kokkos::View<double> r) : result(r) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_type& team, double& update) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double& dst, const volatile double& update) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final(double& update) const { result() = update; }

  KOKKOS_INLINE_FUNCTION
  void init(double& update) const { update = 0.0; }
};
#endif

struct Functor1 {
  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, double& update) const { update += i; }
};

struct Functor2 {
  using value_type = double[];

  const unsigned value_count;

  Functor2(unsigned n) : value_count(n) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned& i, double update[]) const {
    for (unsigned j = 0; j < value_count; j++) {
      update[j] += i;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(double dst[]) const {
    for (unsigned i = 0; i < value_count; ++i) dst[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join(volatile double dst[], const volatile double src[]) const {
    for (unsigned i = 0; i < value_count; ++i) dst[i] += src[i];
  }
};

}  // namespace ReduceCombinatorical

template <class ExecSpace = Kokkos::DefaultExecutionSpace>
struct TestReduceCombinatoricalInstantiation {
  template <class... Args>
  static void CallParallelReduce(Args... args) {
    Kokkos::parallel_reduce(args...);
  }

  template <class... Args>
  static void AddReturnArgument(int N, Args... args) {
    Kokkos::View<double, Kokkos::HostSpace> result_view("ResultViewHost");
    Kokkos::View<double, ExecSpace> result_view_device("ResultViewDevice");
    double expected_result = (1.0 * N) * (1.0 * N - 1.0) / 2.0;

    double value = 99;
    Kokkos::parallel_reduce(args..., value);
    ASSERT_EQ(expected_result, value);

    result_view() = 99;
    CallParallelReduce(args..., result_view);
    Kokkos::fence();
    ASSERT_EQ(expected_result, result_view());

#ifndef KOKKOS_ENABLE_OPENMPTARGET
    result_view() = 99;
    CallParallelReduce(args..., result_view_device);
    Kokkos::fence();
    Kokkos::deep_copy(result_view, result_view_device);
    ASSERT_EQ(expected_result, result_view());
#endif

    value = 99;
    CallParallelReduce(
        args...,
        Kokkos::View<double, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >(&value));
    Kokkos::fence();
    ASSERT_EQ(expected_result, value);

    result_view() = 99;
    const Kokkos::View<double, Kokkos::HostSpace,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> >
        result_view_const_um = result_view;
    CallParallelReduce(args..., result_view_const_um);
    Kokkos::fence();
    ASSERT_EQ(expected_result, result_view_const_um());

    value = 99;
// WORKAROUND OPENMPTARGET Custom Reducers not implemented
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    CallParallelReduce(args...,
                       Test::ReduceCombinatorical::AddPlus<double>(value));
    if ((Kokkos::DefaultExecutionSpace::concurrency() > 1) &&
        (ExecSpace::concurrency() > 1) && (expected_result > 0)) {
      ASSERT_TRUE(expected_result < value);
    } else if (((Kokkos::DefaultExecutionSpace::concurrency() > 1) ||
                (ExecSpace::concurrency() > 1)) &&
               (expected_result > 0)) {
      ASSERT_TRUE(expected_result <= value);
    } else {
      ASSERT_EQ(expected_result, value);
    }

    value = 99;
    Test::ReduceCombinatorical::AddPlus<double> add(value);
    CallParallelReduce(args..., add);
    if ((Kokkos::DefaultExecutionSpace::concurrency() > 1) &&
        (ExecSpace::concurrency() > 1) && (expected_result > 0)) {
      ASSERT_TRUE(expected_result < value);
    } else if (((Kokkos::DefaultExecutionSpace::concurrency() > 1) ||
                (ExecSpace::concurrency() > 1)) &&
               (expected_result > 0)) {
      ASSERT_TRUE(expected_result <= value);
    } else {
      ASSERT_EQ(expected_result, value);
    }
#endif
  }

  template <class... Args>
  static void AddLambdaRange(int N, void*, Args... args) {
    AddReturnArgument(
        N, args..., KOKKOS_LAMBDA(const int& i, double& lsum) { lsum += i; });
  }

  template <class... Args>
  static void AddLambdaTeam(int N, void*, Args... args) {
    AddReturnArgument(
        N, args...,
        KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team,
                      double& update) {
          update += 1.0 / team.team_size() * team.league_rank();
        });
  }

  template <class... Args>
  static void AddLambdaRange(int, Kokkos::InvalidType, Args... /*args*/) {}

  template <class... Args>
  static void AddLambdaTeam(int, Kokkos::InvalidType, Args... /*args*/) {}

  template <int ISTEAM, class... Args>
  static void AddFunctor(int N, Args... args) {
    Kokkos::View<double, ExecSpace> result_view("FunctorView");
    auto h_r = Kokkos::create_mirror_view(result_view);
    Test::ReduceCombinatorical::FunctorScalar<ISTEAM> functor(result_view);

    AddReturnArgument(N, args..., functor);
    AddReturnArgument(
        N, args...,
        Test::ReduceCombinatorical::FunctorScalar<ISTEAM>(result_view));
// WORKAROUND OPENMPTARGET: reductions with functor join/init/final
// not implemented
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
    AddReturnArgument(
        N, args...,
        Test::ReduceCombinatorical::FunctorScalarInit<ISTEAM>(result_view));
    AddReturnArgument(
        N, args...,
        Test::ReduceCombinatorical::FunctorScalarJoin<ISTEAM>(result_view));
    AddReturnArgument(
        N, args...,
        Test::ReduceCombinatorical::FunctorScalarJoinInit<ISTEAM>(result_view));
    double expected_result = (1.0 * N) * (1.0 * N - 1.0) / 2.0;

    h_r() = 0;
    Kokkos::deep_copy(result_view, h_r);
    CallParallelReduce(
        args...,
        Test::ReduceCombinatorical::FunctorScalarFinal<ISTEAM>(result_view));
    Kokkos::fence();
    Kokkos::deep_copy(h_r, result_view);
    ASSERT_EQ(expected_result, h_r());

    h_r() = 0;
    Kokkos::deep_copy(result_view, h_r);
    CallParallelReduce(
        args..., Test::ReduceCombinatorical::FunctorScalarJoinFinal<ISTEAM>(
                     result_view));
    Kokkos::fence();
    Kokkos::deep_copy(h_r, result_view);
    ASSERT_EQ(expected_result, h_r());

    h_r() = 0;
    Kokkos::deep_copy(result_view, h_r);
    CallParallelReduce(
        args..., Test::ReduceCombinatorical::FunctorScalarJoinFinalInit<ISTEAM>(
                     result_view));
    Kokkos::fence();
    Kokkos::deep_copy(h_r, result_view);
    ASSERT_EQ(expected_result, h_r());
#endif
  }

  template <class... Args>
  static void AddFunctorLambdaRange(int N, Args... args) {
    AddFunctor<0, Args...>(N, args...);
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
    AddLambdaRange(
        N,
        typename std::conditional<
            std::is_same<ExecSpace, Kokkos::DefaultExecutionSpace>::value,
            void*, Kokkos::InvalidType>::type(),
        args...);
#endif
  }

  template <class... Args>
  static void AddFunctorLambdaTeam(int N, Args... args) {
    AddFunctor<1, Args...>(N, args...);
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
    AddLambdaTeam(
        N,
        typename std::conditional<
            std::is_same<ExecSpace, Kokkos::DefaultExecutionSpace>::value,
            void*, Kokkos::InvalidType>::type(),
        args...);
#endif
  }

  template <class... Args>
  static void AddPolicy_1(int N, Args... args) {
    Kokkos::RangePolicy<ExecSpace> policy(0, N);

    AddFunctorLambdaRange(1000, args..., 1000);
    AddFunctorLambdaRange(N, args..., N);
    AddFunctorLambdaRange(N, args..., policy);
  }

  template <class... Args>
  static void AddPolicy_2(int N, Args... args) {
    AddFunctorLambdaRange(N, args..., Kokkos::RangePolicy<ExecSpace>(0, N));
    AddFunctorLambdaRange(
        N, args...,
        Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >(0,
                                                                           N));
    AddFunctorLambdaRange(
        N, args...,
        Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Static> >(0, N)
            .set_chunk_size(16));
    AddFunctorLambdaRange(
        N, args...,
        Kokkos::RangePolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >(0, N)
            .set_chunk_size(16));
  }

  template <class... Args>
  static void AddPolicy_3(int N, Args... args) {
    AddFunctorLambdaTeam(N, args...,
                         Kokkos::TeamPolicy<ExecSpace>(N, Kokkos::AUTO));
    AddFunctorLambdaTeam(
        N, args...,
        Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >(
            N, Kokkos::AUTO));
    AddFunctorLambdaTeam(
        N, args...,
        Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Static> >(
            N, Kokkos::AUTO)
            .set_chunk_size(16));
    AddFunctorLambdaTeam(
        N, args...,
        Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >(
            N, Kokkos::AUTO)
            .set_chunk_size(16));
  }

  static void execute_a1() { AddPolicy_1(1000); }

  static void execute_b1() {
    std::string s("Std::String");
    AddPolicy_1(1000, s.c_str());
    AddPolicy_1(1000, "Char Constant");
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    AddPolicy_1(0, "Char Constant");
#endif
  }

  static void execute_c1() {
    std::string s("Std::String");
    AddPolicy_1(1000, s);
  }

  static void execute_a2() { AddPolicy_2(1000); }

  static void execute_b2() {
    std::string s("Std::String");
    AddPolicy_2(1000, s.c_str());
    AddPolicy_2(1000, "Char Constant");
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    AddPolicy_2(0, "Char Constant");
#endif
  }

  static void execute_c2() {
    std::string s("Std::String");
    AddPolicy_2(1000, s);
  }

  static void execute_a3() {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    AddPolicy_3(1000);
#endif
  }

  static void execute_b3() {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    std::string s("Std::String");
    AddPolicy_3(1000, s.c_str());
    AddPolicy_3(1000, "Char Constant");
    AddPolicy_3(0, "Char Constant");
#endif
  }

  static void execute_c3() {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    std::string s("Std::String");
    AddPolicy_3(1000, s);
#endif
  }
};

}  // namespace Test
