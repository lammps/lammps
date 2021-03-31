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

#ifndef KOKKOS_BLAS_KERNELS_HPP
#define KOKKOS_BLAS_KERNELS_HPP

#include <type_traits>

namespace Kokkos {

template <class Type>
struct Dot {
  using execution_space = typename Type::execution_space;

  static_assert(static_cast<unsigned>(Type::Rank) == static_cast<unsigned>(1),
                "Dot static_assert Fail: Rank != 1");

  using value_type = double;

#if 1
  typename Type::const_type X;
  typename Type::const_type Y;
#else
  Type X;
  Type Y;
#endif

  Dot(const Type& arg_x, const Type& arg_y) : X(arg_x), Y(arg_y) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, value_type& update) const { update += X[i] * Y[i]; }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& update,
                   const volatile value_type& source) {
    update += source;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& update) { update = 0; }
};

template <class Type>
struct DotSingle {
  using execution_space = typename Type::execution_space;

  static_assert(static_cast<unsigned>(Type::Rank) == static_cast<unsigned>(1),
                "DotSingle static_assert Fail: Rank != 1");

  using value_type = double;

#if 1
  typename Type::const_type X;
#else
  Type X;
#endif

  DotSingle(const Type& arg_x) : X(arg_x) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, value_type& update) const {
    const typename Type::value_type& x = X[i];
    update += x * x;
  }

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& update,
                   const volatile value_type& source) {
    update += source;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type& update) { update = 0; }
};

template <class ScalarType, class VectorType>
struct Scale {
  using execution_space = typename VectorType::execution_space;

  static_assert(static_cast<unsigned>(ScalarType::Rank) ==
                    static_cast<unsigned>(0),
                "Scale static_assert Fail: ScalarType::Rank != 0");

  static_assert(static_cast<unsigned>(VectorType::Rank) ==
                    static_cast<unsigned>(1),
                "Scale static_assert Fail: VectorType::Rank != 1");

#if 1
  typename ScalarType::const_type alpha;
#else
  ScalarType alpha;
#endif

  VectorType Y;

  Scale(const ScalarType& arg_alpha, const VectorType& arg_Y)
      : alpha(arg_alpha), Y(arg_Y) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { Y[i] *= alpha(); }
};

template <class ScalarType, class ConstVectorType, class VectorType>
struct AXPBY {
  using execution_space = typename VectorType::execution_space;

  static_assert(static_cast<unsigned>(ScalarType::Rank) ==
                    static_cast<unsigned>(0),
                "AXPBY static_assert Fail: ScalarType::Rank != 0");

  static_assert(static_cast<unsigned>(ConstVectorType::Rank) ==
                    static_cast<unsigned>(1),
                "AXPBY static_assert Fail: ConstVectorType::Rank != 1");

  static_assert(static_cast<unsigned>(VectorType::Rank) ==
                    static_cast<unsigned>(1),
                "AXPBY static_assert Fail: VectorType::Rank != 1");

#if 1
  typename ScalarType::const_type alpha, beta;
  typename ConstVectorType::const_type X;
#else
  ScalarType alpha, beta;
  ConstVectorType X;
#endif

  VectorType Y;

  AXPBY(const ScalarType& arg_alpha, const ConstVectorType& arg_X,
        const ScalarType& arg_beta, const VectorType& arg_Y)
      : alpha(arg_alpha), beta(arg_beta), X(arg_X), Y(arg_Y) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { Y[i] = alpha() * X[i] + beta() * Y[i]; }
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
/** \brief  Y = alpha * X + beta * Y */
template <class ConstScalarType, class ConstVectorType, class VectorType>
void axpby(const ConstScalarType& alpha, const ConstVectorType& X,
           const ConstScalarType& beta, const VectorType& Y) {
  using functor = AXPBY<ConstScalarType, ConstVectorType, VectorType>;

  parallel_for(Y.extent(0), functor(alpha, X, beta, Y));
}

/** \brief  Y *= alpha */
template <class ConstScalarType, class VectorType>
void scale(const ConstScalarType& alpha, const VectorType& Y) {
  using functor = Scale<ConstScalarType, VectorType>;

  parallel_for(Y.extent(0), functor(alpha, Y));
}

template <class ConstVectorType, class Finalize>
void dot(const ConstVectorType& X, const ConstVectorType& Y,
         const Finalize& finalize) {
  using functor = Dot<ConstVectorType>;

  parallel_reduce(X.extent(0), functor(X, Y), finalize);
}

template <class ConstVectorType, class Finalize>
void dot(const ConstVectorType& X, const Finalize& finalize) {
  using functor = DotSingle<ConstVectorType>;

  parallel_reduce(X.extent(0), functor(X), finalize);
}

} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_BLAS_KERNELS_HPP */
