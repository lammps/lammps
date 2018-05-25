/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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
#ifndef KOKKOS_CUDA_VECTORIZATION_HPP
#define KOKKOS_CUDA_VECTORIZATION_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA

#include <Kokkos_Cuda.hpp>
#include <Cuda/Kokkos_Cuda_Version_9_8_Compatibility.hpp>
namespace Kokkos {


// Shuffle only makes sense on >= Kepler GPUs; it doesn't work on CPUs
// or other GPUs.  We provide a generic definition (which is trivial
// and doesn't do what it claims to do) because we don't actually use
// this function unless we are on a suitable GPU, with a suitable
// Scalar type.  (For example, in the mat-vec, the "ThreadsPerRow"
// internal parameter depends both on the ExecutionSpace and the Scalar type,
// and it controls whether shfl_down() gets called.)
namespace Impl {

  template< typename Scalar >
  struct shfl_union {
    enum {n = sizeof(Scalar)/4};
    float fval[n];
    KOKKOS_INLINE_FUNCTION
    Scalar value() {
      return *(Scalar*) fval;
    }
    KOKKOS_INLINE_FUNCTION
    void operator= (Scalar& value_) {
      float* const val_ptr = (float*) &value_;
      for(int i=0; i<n ; i++) {
        fval[i] = val_ptr[i];
      }
    }
    KOKKOS_INLINE_FUNCTION
    void operator= (const Scalar& value_) {
      float* const val_ptr = (float*) &value_;
      for(int i=0; i<n ; i++) {
        fval[i] = val_ptr[i];
      }
    }

  };
}

#ifdef __CUDA_ARCH__
#if (__CUDA_ARCH__ >= 300)

    KOKKOS_INLINE_FUNCTION
    int shfl(const int &val, const int& srcLane, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL(val,srcLane,width);
    }

    KOKKOS_INLINE_FUNCTION
    float shfl(const float &val, const int& srcLane, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL(val,srcLane,width);
    }

// TODO: figure out why 64-bit shfl fails with Clang
#if ( CUDA_VERSION >= 9000 ) && (!defined(KOKKOS_COMPILER_CLANG))

    KOKKOS_INLINE_FUNCTION
    long shfl(const long &val, const int& srcLane, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL(val,srcLane,width);
    }

    KOKKOS_INLINE_FUNCTION
    long long shfl(const long long &val, const int& srcLane, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL(val,srcLane,width);
    }

    KOKKOS_INLINE_FUNCTION
    double shfl(const double &val, const int& srcLane, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL(val,srcLane,width);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl(const Scalar &val, const int& srcLane, const typename Impl::enable_if< (sizeof(Scalar) == 8) , int >::type& width
        ) {
      Scalar tmp1 = val;
      double tmp = *reinterpret_cast<double*>(&tmp1);
      tmp = KOKKOS_IMPL_CUDA_SHFL(tmp,srcLane,width);
      return *reinterpret_cast<Scalar*>(&tmp);
    }

#else // ( CUDA_VERSION < 9000 )

    KOKKOS_INLINE_FUNCTION
    double shfl(const double &val, const int& srcLane, const int& width) {
      int lo = __double2loint(val);
      int hi = __double2hiint(val);
      lo = KOKKOS_IMPL_CUDA_SHFL(lo,srcLane,width);
      hi = KOKKOS_IMPL_CUDA_SHFL(hi,srcLane,width);
      return __hiloint2double(hi,lo);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl(const Scalar &val, const int& srcLane, const typename Impl::enable_if< (sizeof(Scalar) == 8) ,int>::type& width) {
      int lo = __double2loint(*reinterpret_cast<const double*>(&val));
      int hi = __double2hiint(*reinterpret_cast<const double*>(&val));
      lo = KOKKOS_IMPL_CUDA_SHFL(lo,srcLane,width);
      hi = KOKKOS_IMPL_CUDA_SHFL(hi,srcLane,width);
      const double tmp = __hiloint2double(hi,lo);
      return *(reinterpret_cast<const Scalar*>(&tmp));
    }

#endif // ( CUDA_VERSION < 9000 )

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl(const Scalar &val, const int& srcLane, const typename Impl::enable_if< (sizeof(Scalar) == 4) , int >::type& width
        ) {
      Scalar tmp1 = val;
      float tmp = *reinterpret_cast<float*>(&tmp1);
      tmp = KOKKOS_IMPL_CUDA_SHFL(tmp,srcLane,width);
      return *reinterpret_cast<Scalar*>(&tmp);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl(const Scalar &val, const int& srcLane, const typename Impl::enable_if< (sizeof(Scalar) > 8) ,int>::type& width) {
      Impl::shfl_union<Scalar> s_val;
      Impl::shfl_union<Scalar> r_val;
      s_val = val;

      for(int i = 0; i<s_val.n; i++)
        r_val.fval[i] = KOKKOS_IMPL_CUDA_SHFL(s_val.fval[i],srcLane,width);
      return r_val.value();
    }

    KOKKOS_INLINE_FUNCTION
    int shfl_down(const int &val, const int& delta, const int& width) {
      return KOKKOS_IMPL_CUDA_SHFL_DOWN(val,delta,width);
    }

    KOKKOS_INLINE_FUNCTION
    float shfl_down(const float &val, const int& delta, const int& width) {
      return KOKKOS_IMPL_CUDA_SHFL_DOWN(val,delta,width);
    }

// TODO: figure out why 64-bit shfl fails with Clang
#if ( CUDA_VERSION >= 9000 ) && (!defined(KOKKOS_COMPILER_CLANG))

    KOKKOS_INLINE_FUNCTION
    long shfl_down(const long &val, const int& delta, const int& width) {
      return KOKKOS_IMPL_CUDA_SHFL_DOWN(val,delta,width);
    }

    KOKKOS_INLINE_FUNCTION
    long long shfl_down(const long long &val, const int& delta, const int& width) {
      return KOKKOS_IMPL_CUDA_SHFL_DOWN(val,delta,width);
    }

    KOKKOS_INLINE_FUNCTION
    double shfl_down(const double &val, const int& delta, const int& width) {
      return KOKKOS_IMPL_CUDA_SHFL_DOWN(val,delta,width);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_down(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) == 8) , int >::type & width) {
      Scalar tmp1 = val;
      double tmp = *reinterpret_cast<double*>(&tmp1);
      tmp = KOKKOS_IMPL_CUDA_SHFL_DOWN(tmp,delta,width);
      return *reinterpret_cast<Scalar*>(&tmp);
    }

#else // ( CUDA_VERSION < 9000 )

    KOKKOS_INLINE_FUNCTION
    double shfl_down(const double &val, const int& delta, const int& width) {
      int lo = __double2loint(val);
      int hi = __double2hiint(val);
      lo = KOKKOS_IMPL_CUDA_SHFL_DOWN(lo,delta,width);
      hi = KOKKOS_IMPL_CUDA_SHFL_DOWN(hi,delta,width);
      return __hiloint2double(hi,lo);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_down(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) == 8) , int >::type & width) {
      int lo = __double2loint(*reinterpret_cast<const double*>(&val));
      int hi = __double2hiint(*reinterpret_cast<const double*>(&val));
      lo = KOKKOS_IMPL_CUDA_SHFL_DOWN(lo,delta,width);
      hi = KOKKOS_IMPL_CUDA_SHFL_DOWN(hi,delta,width);
      const double tmp = __hiloint2double(hi,lo);
      return *(reinterpret_cast<const Scalar*>(&tmp));
    }

#endif // ( CUDA_VERSION < 9000 )

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_down(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) == 4) , int >::type & width) {
      Scalar tmp1 = val;
      float tmp = *reinterpret_cast<float*>(&tmp1);
      tmp = KOKKOS_IMPL_CUDA_SHFL_DOWN(tmp,delta,width);
      return *reinterpret_cast<Scalar*>(&tmp);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_down(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) > 8) , int >::type & width) {
      Impl::shfl_union<Scalar> s_val;
      Impl::shfl_union<Scalar> r_val;
      s_val = val;

      for(int i = 0; i<s_val.n; i++)
        r_val.fval[i] = KOKKOS_IMPL_CUDA_SHFL_DOWN(s_val.fval[i],delta,width);
      return r_val.value();
    }

    KOKKOS_INLINE_FUNCTION
    int shfl_up(const int &val, const int& delta, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL_UP(val,delta,width);
    }

    KOKKOS_INLINE_FUNCTION
    float shfl_up(const float &val, const int& delta, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL_UP(val,delta,width);
    }

// TODO: figure out why 64-bit shfl fails with Clang
#if ( CUDA_VERSION >= 9000 ) && (!defined(KOKKOS_COMPILER_CLANG))

    KOKKOS_INLINE_FUNCTION
    long shfl_up(const long &val, const int& delta, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL_UP(val,delta,width);
    }

    KOKKOS_INLINE_FUNCTION
    long long shfl_up(const long long &val, const int& delta, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL_UP(val,delta,width);
    }

    KOKKOS_INLINE_FUNCTION
    double shfl_up(const double &val, const int& delta, const int& width ) {
      return KOKKOS_IMPL_CUDA_SHFL_UP(val,delta,width);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_up(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) == 8) , int >::type & width) {
      Scalar tmp1 = val;
      double tmp = *reinterpret_cast<double*>(&tmp1);
      tmp = KOKKOS_IMPL_CUDA_SHFL_UP(tmp,delta,width);
      return *reinterpret_cast<Scalar*>(&tmp);
    }

#else // ( CUDA_VERSION < 9000 )

    KOKKOS_INLINE_FUNCTION
    double shfl_up(const double &val, const int& delta, const int& width ) {
      int lo = __double2loint(val);
      int hi = __double2hiint(val);
      lo = KOKKOS_IMPL_CUDA_SHFL_UP(lo,delta,width);
      hi = KOKKOS_IMPL_CUDA_SHFL_UP(hi,delta,width);
      return __hiloint2double(hi,lo);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_up(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) == 8) , int >::type & width) {
      int lo = __double2loint(*reinterpret_cast<const double*>(&val));
      int hi = __double2hiint(*reinterpret_cast<const double*>(&val));
      lo = KOKKOS_IMPL_CUDA_SHFL_UP(lo,delta,width);
      hi = KOKKOS_IMPL_CUDA_SHFL_UP(hi,delta,width);
      const double tmp = __hiloint2double(hi,lo);
      return *(reinterpret_cast<const Scalar*>(&tmp));
    }

#endif // ( CUDA_VERSION < 9000 )

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_up(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) == 4) , int >::type & width) {
      Scalar tmp1 = val;
      float tmp = *reinterpret_cast<float*>(&tmp1);
      tmp = KOKKOS_IMPL_CUDA_SHFL_UP(tmp,delta,width);
      return *reinterpret_cast<Scalar*>(&tmp);
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_up(const Scalar &val, const int& delta, const typename Impl::enable_if< (sizeof(Scalar) > 8) , int >::type & width) {
      Impl::shfl_union<Scalar> s_val;
      Impl::shfl_union<Scalar> r_val;
      s_val = val;

      for(int i = 0; i<s_val.n; i++)
        r_val.fval[i] = KOKKOS_IMPL_CUDA_SHFL_UP(s_val.fval[i],delta,width);
      return r_val.value();
    }

#else // (__CUDA_ARCH__ < 300)
    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl(const Scalar &val, const int& srcLane, const int& width) {
      if(width > 1) Kokkos::abort("Error: calling shfl from a device with CC<3.0.");
      return val;
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_down(const Scalar &val, const int& delta, const int& width) {
      if(width > 1) Kokkos::abort("Error: calling shfl_down from a device with CC<3.0.");
      return val;
    }

    template<typename Scalar>
    KOKKOS_INLINE_FUNCTION
    Scalar shfl_up(const Scalar &val, const int& delta, const int& width) {
      if(width > 1) Kokkos::abort("Error: calling shfl_down from a device with CC<3.0.");
      return val;
    }
#endif // (__CUDA_ARCH__ < 300)
#else // !defined( __CUDA_ARCH__ )
    template<typename Scalar>
    inline
    Scalar shfl(const Scalar &val, const int& srcLane, const int& width) {
      if(width > 1) Kokkos::abort("Error: calling shfl outside __CUDA_ARCH__.");
      return val;
    }

    template<typename Scalar>
    inline
    Scalar shfl_down(const Scalar &val, const int& delta, const int& width) {
      if(width > 1) Kokkos::abort("Error: calling shfl_down outside __CUDA_ARCH__.");
      return val;
    }

    template<typename Scalar>
    inline
    Scalar shfl_up(const Scalar &val, const int& delta, const int& width) {
      if(width > 1) Kokkos::abort("Error: calling shfl_down outside __CUDA_ARCH__.");
      return val;
    }
#endif // !defined( __CUDA_ARCH__ )

} // end namespace Kokkos

#endif // defined( KOKKOS_ENABLE_CUDA )
#endif // !defined( KOKKOS_CUDA_VECTORIZATION_HPP )

