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

#ifndef KOKKOS_HIP_VECTORIZATION_HPP
#define KOKKOS_HIP_VECTORIZATION_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_HIP.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

// Include all lanes
constexpr unsigned shfl_all_mask = 0xffffffff;

//----------------------------------------------------------------------------
// Shuffle operations require input to be a register (stack) variable

// Derived implements do_shfl_op(unsigned mask, T& in, int lane, int width),
// which turns in to one of KOKKOS_IMPL_HIP_SHFL(_UP_|_DOWN_|_)MASK
// Since the logic with respect to value sizes, etc., is the same everywhere,
// put it all in one place.
template <class Derived>
struct in_place_shfl_op {
  // CRTP boilerplate
  __device__ KOKKOS_IMPL_FORCEINLINE const Derived& self() const noexcept {
    return *static_cast<Derived const*>(this);
  }

  // sizeof(Scalar) == sizeof(int) case
  template <class Scalar>
  // requires _assignable_from_bits<Scalar>
  __device__ inline typename std::enable_if<sizeof(Scalar) == sizeof(int)>::type
  operator()(Scalar& out, Scalar const& in, int lane_or_delta, int width,
             unsigned mask = shfl_all_mask) const noexcept {
    //------------------------------------------------
    reinterpret_cast<int&>(out) = self().do_shfl_op(
        mask, reinterpret_cast<int const&>(in), lane_or_delta, width);
    //------------------------------------------------
  }

  template <class Scalar>
  __device__ inline
      typename std::enable_if<sizeof(Scalar) == sizeof(double)>::type
      operator()(Scalar& out, Scalar const& in, int lane_or_delta, int width,
                 unsigned mask = shfl_all_mask) const noexcept {
    //------------------------------------------------
    reinterpret_cast<double&>(out) = self().do_shfl_op(
        mask, *reinterpret_cast<double const*>(&in), lane_or_delta, width);
    //------------------------------------------------
  }

  // sizeof(Scalar) > sizeof(double) case
  template <typename Scalar>
  __device__ inline
      typename std::enable_if<(sizeof(Scalar) > sizeof(double))>::type
      operator()(Scalar& out, const Scalar& val, int lane_or_delta, int width,
                 unsigned mask = shfl_all_mask) const noexcept {
    using shuffle_as_t = int;
    enum : int { N = sizeof(Scalar) / sizeof(shuffle_as_t) };

    for (int i = 0; i < N; ++i) {
      reinterpret_cast<shuffle_as_t*>(&out)[i] = self().do_shfl_op(
          mask, reinterpret_cast<shuffle_as_t const*>(&val)[i], lane_or_delta,
          width);
    }
  }
};

struct in_place_shfl_fn : in_place_shfl_op<in_place_shfl_fn> {
  template <class T>
  __device__ KOKKOS_IMPL_FORCEINLINE T do_shfl_op(unsigned mask, T& val,
                                                  int lane, int width) const
      noexcept {
    return KOKKOS_IMPL_HIP_SHFL_MASK(mask, val, lane, width);
  }
};

template <class... Args>
__device__ KOKKOS_IMPL_FORCEINLINE void in_place_shfl(Args&&... args) noexcept {
  in_place_shfl_fn{}((Args &&) args...);
}

struct in_place_shfl_up_fn : in_place_shfl_op<in_place_shfl_up_fn> {
  template <class T>
  __device__ KOKKOS_IMPL_FORCEINLINE T do_shfl_op(unsigned mask, T& val,
                                                  int lane, int width) const
      noexcept {
    return KOKKOS_IMPL_HIP_SHFL_UP_MASK(mask, val, lane, width);
  }
};

template <class... Args>
__device__ KOKKOS_IMPL_FORCEINLINE void in_place_shfl_up(
    Args&&... args) noexcept {
  in_place_shfl_up_fn{}((Args &&) args...);
}

struct in_place_shfl_down_fn : in_place_shfl_op<in_place_shfl_down_fn> {
  template <class T>
  __device__ KOKKOS_IMPL_FORCEINLINE T do_shfl_op(unsigned mask, T& val,
                                                  int lane, int width) const
      noexcept {
    return KOKKOS_IMPL_HIP_SHFL_DOWN_MASK(mask, val, lane, width);
  }
};

template <class... Args>
__device__ KOKKOS_IMPL_FORCEINLINE void in_place_shfl_down(
    Args&&... args) noexcept {
  in_place_shfl_down_fn{}((Args &&) args...);
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
