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
#ifndef KOKKOS_ROCM_VECTORIZATION_HPP
#define KOKKOS_ROCM_VECTORIZATION_HPP

#include <Kokkos_Macros.hpp>

/* only compile this file if ROCM is enabled for Kokkos */
#ifdef KOKKOS_ENABLE_ROCM

#include <Kokkos_ROCm.hpp>

namespace Kokkos {
using namespace hc;

// Shuffle only makes sense on >= Fiji GPUs; it doesn't work on CPUs
// or other GPUs.  We provide a generic definition (which is trivial
// and doesn't do what it claims to do) because we don't actually use
// this function unless we are on a suitable GPU, with a suitable
// Scalar type.  (For example, in the mat-vec, the "ThreadsPerRow"
// internal parameter depends both on the ExecutionSpace and the Scalar type,
// and it controls whether shfl_down() gets called.)
namespace Impl {

template <typename Scalar>
struct shfl_union {
  enum { n = sizeof(Scalar) / 4 };
  float fval[n];
  KOKKOS_INLINE_FUNCTION
  Scalar value() { return *(Scalar*)fval; }
  KOKKOS_INLINE_FUNCTION
  void operator=(Scalar& value_) {
    float* const val_ptr = (float*)&value_;
    for (int i = 0; i < n; i++) {
      fval[i] = val_ptr[i];
    }
  }
  KOKKOS_INLINE_FUNCTION
  void operator=(const Scalar& value_) {
    float* const val_ptr = (float*)&value_;
    for (int i = 0; i < n; i++) {
      fval[i] = val_ptr[i];
    }
  }
};
}  // namespace Impl

#ifdef __HCC_ACCELERATOR__

KOKKOS_INLINE_FUNCTION
int __long2loint(const long val) {
  union {
    long l;
    int i[2];
  } u;
  u.l = val;
  return u.i[0];
}

KOKKOS_INLINE_FUNCTION
int __long2hiint(const long val) {
  union {
    long l;
    int i[2];
  } u;
  u.l = val;
  return u.i[1];
}

KOKKOS_INLINE_FUNCTION
int __double2loint(const double val) {
  union {
    double d;
    int i[2];
  } u;
  u.d = val;
  return u.i[0];
}

KOKKOS_INLINE_FUNCTION
int __double2hiint(const double val) {
  union {
    double d;
    int i[2];
  } u;
  u.d = val;
  return u.i[1];
}

KOKKOS_INLINE_FUNCTION
long __hiloint2long(const int hi, const int lo) {
  union {
    long l;
    int i[2];
  } u;
  u.i[0] = lo;
  u.i[1] = hi;
  return u.l;
}

KOKKOS_INLINE_FUNCTION
double __hiloint2double(const int hi, const int lo) {
  union {
    double d;
    int i[2];
  } u;
  u.i[0] = lo;
  u.i[1] = hi;
  return u.d;
}

KOKKOS_INLINE_FUNCTION
int shfl(const int& val, const int& srcLane, const int& width) {
  return __shfl(val, srcLane, width);
}

KOKKOS_INLINE_FUNCTION
float shfl(const float& val, const int& srcLane, const int& width) {
  return __shfl(val, srcLane, width);
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
shfl(const Scalar& val, const int& srcLane,
     const typename Impl::enable_if<(sizeof(Scalar) == 4), int>::type& width) {
  Scalar tmp1 = val;
  float tmp   = *reinterpret_cast<float*>(&tmp1);
  tmp         = __shfl(tmp, srcLane, width);
  return *reinterpret_cast<Scalar*>(&tmp);
}

KOKKOS_INLINE_FUNCTION
double shfl(const double& val, const int& srcLane, const int& width) {
  int lo = __double2loint(val);
  int hi = __double2hiint(val);
  lo     = __shfl(lo, srcLane, width);
  hi     = __shfl(hi, srcLane, width);
  return __hiloint2double(hi, lo);
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
shfl(const Scalar& val, const int& srcLane,
     const typename Impl::enable_if<(sizeof(Scalar) == 8), int>::type& width) {
  int lo           = __double2loint(*reinterpret_cast<const double*>(&val));
  int hi           = __double2hiint(*reinterpret_cast<const double*>(&val));
  lo               = __shfl(lo, srcLane, width);
  hi               = __shfl(hi, srcLane, width);
  const double tmp = __hiloint2double(hi, lo);
  return *(reinterpret_cast<const Scalar*>(&tmp));
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar
shfl(const Scalar& val, const int& srcLane,
     const typename Impl::enable_if<(sizeof(Scalar) > 8), int>::type& width) {
  Impl::shfl_union<Scalar> s_val;
  Impl::shfl_union<Scalar> r_val;
  s_val = val;

  for (int i = 0; i < s_val.n; i++)
    r_val.fval[i] = __shfl(s_val.fval[i], srcLane, width);
  return r_val.value();
}

KOKKOS_INLINE_FUNCTION
int shfl_down(const int& val, const int& delta, const int& width) {
  return __shfl_down(val, delta, width);
}

KOKKOS_INLINE_FUNCTION
float shfl_down(const float& val, const int& delta, const int& width) {
  return __shfl_down(val, delta, width);
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar shfl_down(
    const Scalar& val, const int& delta,
    const typename Impl::enable_if<(sizeof(Scalar) == 4), int>::type& width) {
  Scalar tmp1 = val;
  float tmp   = *reinterpret_cast<float*>(&tmp1);
  tmp         = __shfl_down(tmp, delta, width);
  return *reinterpret_cast<Scalar*>(&tmp);
}

KOKKOS_INLINE_FUNCTION
long shfl_down(const long& val, const int& delta, const int& width) {
  int lo = __long2loint(val);
  int hi = __long2hiint(val);
  lo     = __shfl_down(lo, delta, width);
  hi     = __shfl_down(hi, delta, width);
  return __hiloint2long(hi, lo);
}

KOKKOS_INLINE_FUNCTION
double shfl_down(const double& val, const int& delta, const int& width) {
  int lo = __double2loint(val);
  int hi = __double2hiint(val);
  lo     = __shfl_down(lo, delta, width);
  hi     = __shfl_down(hi, delta, width);
  return __hiloint2double(hi, lo);
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar shfl_down(
    const Scalar& val, const int& delta,
    const typename Impl::enable_if<(sizeof(Scalar) == 8), int>::type& width) {
  int lo           = __double2loint(*reinterpret_cast<const double*>(&val));
  int hi           = __double2hiint(*reinterpret_cast<const double*>(&val));
  lo               = __shfl_down(lo, delta, width);
  hi               = __shfl_down(hi, delta, width);
  const double tmp = __hiloint2double(hi, lo);
  return *(reinterpret_cast<const Scalar*>(&tmp));
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar shfl_down(
    const Scalar& val, const int& delta,
    const typename Impl::enable_if<(sizeof(Scalar) > 8), int>::type& width) {
  Impl::shfl_union<Scalar> s_val;
  Impl::shfl_union<Scalar> r_val;
  s_val = val;

  for (int i = 0; i < s_val.n; i++)
    r_val.fval[i] = __shfl_down(s_val.fval[i], delta, width);
  return r_val.value();
}

KOKKOS_INLINE_FUNCTION
int shfl_up(const int& val, const int& delta, const int& width) {
  return __shfl_up(val, delta, width);
}

KOKKOS_INLINE_FUNCTION
float shfl_up(const float& val, const int& delta, const int& width) {
  return __shfl_up(val, delta, width);
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar shfl_up(
    const Scalar& val, const int& delta,
    const typename Impl::enable_if<(sizeof(Scalar) == 4), int>::type& width) {
  Scalar tmp1 = val;
  float tmp   = *reinterpret_cast<float*>(&tmp1);
  tmp         = __shfl_up(tmp, delta, width);
  return *reinterpret_cast<Scalar*>(&tmp);
}

KOKKOS_INLINE_FUNCTION
double shfl_up(const double& val, const int& delta, const int& width) {
  int lo = __double2loint(val);
  int hi = __double2hiint(val);
  lo     = __shfl_up(lo, delta, width);
  hi     = __shfl_up(hi, delta, width);
  return __hiloint2double(hi, lo);
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar shfl_up(
    const Scalar& val, const int& delta,
    const typename Impl::enable_if<(sizeof(Scalar) == 8), int>::type& width) {
  int lo           = __double2loint(*reinterpret_cast<const double*>(&val));
  int hi           = __double2hiint(*reinterpret_cast<const double*>(&val));
  lo               = __shfl_up(lo, delta, width);
  hi               = __shfl_up(hi, delta, width);
  const double tmp = __hiloint2double(hi, lo);
  return *(reinterpret_cast<const Scalar*>(&tmp));
}

template <typename Scalar>
KOKKOS_INLINE_FUNCTION Scalar shfl_up(
    const Scalar& val, const int& delta,
    const typename Impl::enable_if<(sizeof(Scalar) > 8), int>::type& width) {
  Impl::shfl_union<Scalar> s_val;
  Impl::shfl_union<Scalar> r_val;
  s_val = val;

  for (int i = 0; i < s_val.n; i++)
    r_val.fval[i] = __shfl_up(s_val.fval[i], delta, width);
  return r_val.value();
}

#else
template <typename Scalar>
inline Scalar shfl(const Scalar& val, const int& srcLane, const int& width) {
  if (width > 1)
    Kokkos::abort("Error: calling shfl from a device with CC<8.0.");
  return val;
}

template <typename Scalar>
inline Scalar shfl_down(const Scalar& val, const int& delta, const int& width) {
  if (width > 1)
    Kokkos::abort("Error: calling shfl_down from a device with CC<8.0.");
  return val;
}

template <typename Scalar>
inline Scalar shfl_up(const Scalar& val, const int& delta, const int& width) {
  if (width > 1)
    Kokkos::abort("Error: calling shfl_down from a device with CC<8.0.");
  return val;
}
#endif

}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_ROCM
#endif
