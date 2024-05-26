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
#ifndef KOKKOS_CUDA_VECTORIZATION_HPP
#define KOKKOS_CUDA_VECTORIZATION_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <type_traits>

#if !defined(KOKKOS_COMPILER_CLANG)
#define KOKKOS_IMPL_CUDA_MAX_SHFL_SIZEOF sizeof(long long)
#else
#define KOKKOS_IMPL_CUDA_MAX_SHFL_SIZEOF sizeof(int)
#endif

namespace Kokkos {

namespace Impl {

// Include all lanes
constexpr unsigned shfl_all_mask = 0xffffffffu;

//----------------------------------------------------------------------------
// Shuffle operations require input to be a register (stack) variable

// Derived implements do_shfl_op(unsigned mask, T& in, int lane, int width),
// which turns in to one of __shfl_sync(_up|_down)
// Since the logic with respect to value sizes, etc., is the same everywhere,
// put it all in one place.
template <class Derived>
struct in_place_shfl_op {
  // CRTP boilerplate
  __device__ KOKKOS_IMPL_FORCEINLINE const Derived& self() const noexcept {
    return *static_cast<Derived const*>(this);
  }

  // sizeof(Scalar) <= sizeof(int) case
  template <class Scalar>
  // requires _assignable_from_bits<Scalar>
  __device__ inline std::enable_if_t<sizeof(Scalar) <= sizeof(int)> operator()(
      Scalar& out, Scalar const& in, int lane_or_delta, int width,
      unsigned mask = shfl_all_mask) const noexcept {
    using shfl_type = int;
    union conv_type {
      Scalar orig;
      shfl_type conv;
      // This should be fine, members get explicitly reset, which changes the
      // active member
      KOKKOS_FUNCTION conv_type() { conv = 0; }
    };
    conv_type tmp_in;
    tmp_in.orig = in;
    shfl_type tmp_out;
    tmp_out = reinterpret_cast<shfl_type&>(tmp_in.orig);
    conv_type res;
    //------------------------------------------------
    res.conv = self().do_shfl_op(mask, tmp_out, lane_or_delta, width);
    //------------------------------------------------
    out = reinterpret_cast<Scalar&>(res.conv);
  }

// TODO: figure out why 64-bit shfl fails in Clang
#if !defined(KOKKOS_COMPILER_CLANG)
  // sizeof(Scalar) == sizeof(double) case
  // requires _assignable_from_bits<Scalar>
  template <class Scalar>
  __device__ inline std::enable_if_t<sizeof(Scalar) == sizeof(double)>
  operator()(Scalar& out, Scalar const& in, int lane_or_delta, int width,
             unsigned mask = shfl_all_mask) const noexcept {
    //------------------------------------------------
    reinterpret_cast<double&>(out) = self().do_shfl_op(
        mask, *reinterpret_cast<double const*>(&in), lane_or_delta, width);
    //------------------------------------------------
  }
#else
  // sizeof(Scalar) == sizeof(double) case
  // requires _assignable_from_bits<Scalar>
  template <typename Scalar>
  __device__ inline std::enable_if_t<sizeof(Scalar) == sizeof(double)>
  operator()(Scalar& out, const Scalar& val, int lane_or_delta, int width,
             unsigned mask = shfl_all_mask) const noexcept {
    //------------------------------------------------
    int lo   = __double2loint(*reinterpret_cast<const double*>(&val));
    int hi   = __double2hiint(*reinterpret_cast<const double*>(&val));
    lo       = self().do_shfl_op(mask, lo, lane_or_delta, width);
    hi       = self().do_shfl_op(mask, hi, lane_or_delta, width);
    auto tmp = __hiloint2double(hi, lo);
    out      = reinterpret_cast<Scalar&>(tmp);
    //------------------------------------------------
  }
#endif

  // sizeof(Scalar) > sizeof(double) case
  template <typename Scalar>
  __device__ inline std::enable_if_t<(sizeof(Scalar) > sizeof(double))>
  operator()(Scalar& out, const Scalar& val, int lane_or_delta, int width,
             unsigned mask = shfl_all_mask) const noexcept {
    // TODO DSH shouldn't this be KOKKOS_IMPL_CUDA_MAX_SHFL_SIZEOF instead of
    //      sizeof(int)? (Need benchmarks to decide which is faster)
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
    (void)mask;
    (void)val;
    (void)lane;
    (void)width;
    return __shfl_sync(mask, val, lane, width);
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
    return __shfl_up_sync(mask, val, lane, width);
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
    (void)mask;
    (void)val;
    (void)lane;
    (void)width;
    return __shfl_down_sync(mask, val, lane, width);
  }
};
template <class... Args>
__device__ KOKKOS_IMPL_FORCEINLINE void in_place_shfl_down(
    Args&&... args) noexcept {
  in_place_shfl_down_fn{}((Args &&) args...);
}

}  // namespace Impl

template <class T>
// requires default_constructible<T> && _assignable_from_bits<T>
__device__ inline T shfl(const T& val, const int& srcLane, const int& width,
                         unsigned mask = Impl::shfl_all_mask) {
  T rv = {};
  Impl::in_place_shfl(rv, val, srcLane, width, mask);
  return rv;
}

template <class T>
// requires default_constructible<T> && _assignable_from_bits<T>
__device__ inline T shfl_down(const T& val, int delta, int width,
                              unsigned mask = Impl::shfl_all_mask) {
  T rv = {};
  Impl::in_place_shfl_down(rv, val, delta, width, mask);
  return rv;
}

template <class T>
// requires default_constructible<T> && _assignable_from_bits<T>
__device__ inline T shfl_up(const T& val, int delta, int width,
                            unsigned mask = Impl::shfl_all_mask) {
  T rv = {};
  Impl::in_place_shfl_up(rv, val, delta, width, mask);
  return rv;
}

}  // end namespace Kokkos

#undef KOKKOS_IMPL_CUDA_MAX_SHFL_SIZEOF

#endif  // defined( KOKKOS_ENABLE_CUDA )
#endif  // !defined( KOKKOS_CUDA_VECTORIZATION_HPP )
