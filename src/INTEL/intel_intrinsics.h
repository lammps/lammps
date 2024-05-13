// clang-format off
/* *- c++ -*- -----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Markus Hoehnerbach (RWTH)
------------------------------------------------------------------------- */

// This file provides an intrinsics abstraction that allows access to the
// underlying SIMD registers.
// It allows for algorithms that are templated over the used vector length
// The final interface is provided by vector_routines, which provides the
// support for different precision modes.
// The vector_ops interface provides routines specific to one floating point
// data type, and is specialized for various architectures.
// The routines work best with AVX-512 and AVX2, as both support the gather
// instructions.
// For both AVX and SSE we miss some optimization opportunities in the gather
// implementations.

// Vector classes provided with the intel compiler
#if defined(__MIC__) && !defined(__AVX512F__)
#include <mic/micvec.h>
#else
#include <dvec.h> // icc-mmic hates generating movq
#include <fvec.h>
#endif

namespace lmp_intel {

// Self explanatory mostly, KNC=IMCI and AVX-512, NONE=Scalar.
enum CalculationMode {KNC, AVX, AVX2, SSE, NONE};
#ifdef __MIC__
  #ifdef LMP_INTEL_VECTOR_MIC
  static const CalculationMode mode = LMP_INTEL_VECTOR_MIC;
  #else
  static const CalculationMode mode = KNC;
  #endif
#else
  #ifdef LMP_INTEL_VECTOR_HOST
  static const CalculationMode mode = LMP_INTEL_VECTOR_HOST;
  #else
    #ifdef __AVX512F__
    static const CalculationMode mode = KNC;
    #else
      #ifdef __AVX2__
      static const CalculationMode mode = AVX2;
      #else
        #ifdef __AVX__
        static const CalculationMode mode = AVX;
        #else
        static const CalculationMode mode = SSE;
        #endif
      #endif
    #endif
  #endif
#endif

// This is used in the selection logic
template<CalculationMode mode>
struct vector_traits {
    static const bool support_integer_and_gather_ops = true;
};

template<>
struct vector_traits<AVX> {
    static const bool support_integer_and_gather_ops = false;
};

// This is the base template for all the different architectures
// It will get specialized
template<class flt_t, CalculationMode mode>
struct vector_ops {};

// Intrinsic routines for IMCI and AVX-512
#if defined(__MIC__) || defined(__AVX512F__)
// Integer vector class
#ifdef __INTEL_LLVM_COMPILER
#pragma pack(push,16)
#else
#pragma pack(push,64)
#endif
struct ivec32x16 {
  __m512i vec;
  ivec32x16() {}
  ivec32x16(__m512i m) { vec = m; }
  ivec32x16(const int * a) {
    vec = _mm512_load_epi32(reinterpret_cast<const int *>(a));
  }
  explicit ivec32x16(int i) { vec = _mm512_set1_epi32(i); }
  operator __m512i() const { return vec; }
  friend ivec32x16 operator &(const ivec32x16 &a, const ivec32x16 &b) {
    return _mm512_and_epi32(a, b);
  }
  friend ivec32x16 operator |(const ivec32x16 &a, const ivec32x16 &b) {
    return _mm512_or_epi32(a, b);
  }
  friend ivec32x16 operator +(const ivec32x16 &a, const ivec32x16 &b) {
    return _mm512_add_epi32(a, b);
  }
};
#pragma pack(pop)
// Double precision routines
template<>
struct vector_ops<double, KNC> {
    static const int VL = 8;
    typedef double fscal;
    typedef F64vec8 fvec;
    typedef ivec32x16 ivec;
    typedef __mmask16 bvec;
    typedef double farr[8] __attribute__((aligned(64)));
    typedef int iarr[16] __attribute__((aligned(64)));
    static fvec recip(const fvec &a) { return _mm512_recip_pd(a); }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, bvec mask, const void *base) {
      _mm512_mask_prefetch_i32gather_ps(idx, mask, base, scale, _MM_HINT_T0);
    }
    template<int scale>
    static fvec gather(const fvec &from, bvec mask, const ivec &idx, const void *base) {
      return _mm512_mask_i32gather_pd(from, mask, _mm512_castsi512_si256(idx),
                                      base, scale);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return _mm512_mask_blend_pd(mask, a, b);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return _mm512_fmadd_pd(a, b, c);
    }
    static fvec zero() {
      return _mm512_setzero_pd();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm512_cmp_pd_mask(a, b, _CMP_EQ_OQ);
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm512_cmp_pd_mask(a, b, _CMP_NLE_US);
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm512_cmp_pd_mask(a, b, _CMP_LE_OS);
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm512_cmp_pd_mask(a, b, _CMP_LT_OS);
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      return _mm512_cmpneq_epi32_mask(a, b);
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return _mm512_cmplt_epi32_mask(a, b);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm512_invsqrt_pd(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      #if __INTEL_COMPILER+0 < 1500
      *reinterpret_cast<__m512d *>(cos) = _mm512_cos_pd(a);
      return _mm512_sin_pd(a);
      #else
      return _mm512_sincos_pd(reinterpret_cast<__m512d *>(cos), a);
      #endif
    }
    static fscal reduce_add(const fvec &a) {
      return _mm512_reduce_add_pd(a);
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return _mm512_mullo_epi32(a, b);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return _mm512_mask_add_epi32(src, mask, a, b);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      return _mm512_mask_i32gather_epi32(from, mask, idx, base, scale);
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return _mm512_mask_add_pd(src, mask, a, b);
    }
    static void store(void *at, const fvec &a) {
      _mm512_store_pd(at, a);
    }
    static void int_store(void *at, const ivec &a) {
      _mm512_store_epi32(at, a);
    }
    static void mask_store(int *at, const bvec &a) {
      for (int i = 0; i < 8; i++) {
        at[i] = (a >> i) & 1;
      }
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm512_min_pd(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return mask & (1 << at);
    }
    static bool mask_testz(const bvec &mask) {
      return mask == 0;
    }

    static bvec mask_enable_lower(int n) {
      return 0xFF >> (VL - n);
    }

    static ivec int_load_vl(const int *a) {
      return _mm512_load_epi32(a);
    }
    static void int_clear_arr(int *a) {
      _mm512_store_epi32(a, ivec(0));
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm512_store_epi32(tmp, a);
      for (int i = 0; i < 8; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      *x = gather<1>(*x, mask, idxs, &base->x);
      *y = gather<1>(*y, mask, idxs, &base->y);
      *z = gather<1>(*z, mask, idxs, &base->z);
      *w = int_gather<1>(*w, mask, idxs, &base->w);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char *>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char *>(base) +  8);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char *>(base) + 16);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char *>(base) + 24);
      *r4 = gather<4>(*r4, mask, idxs, reinterpret_cast<const char *>(base) + 32);
      *r5 = gather<4>(*r5, mask, idxs, reinterpret_cast<const char *>(base) + 40);
      *r6 = gather<4>(*r6, mask, idxs, reinterpret_cast<const char *>(base) + 48);
      *r7 = gather<4>(*r7, mask, idxs, reinterpret_cast<const char *>(base) + 56);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char *>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char *>(base) +  8);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char *>(base) + 16);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char *>(base) + 24);
    }
};

template<>
struct vector_ops<float, KNC> {
    static const int VL = 16;
    static const int ALIGN = 64;
    typedef float fscal;
    typedef F32vec16 fvec;
    typedef ivec32x16 ivec;
    typedef __mmask16 bvec;
    typedef float farr[16] __attribute__((aligned(64)));
    typedef int iarr[16] __attribute__((aligned(64)));
    static const bvec full_mask = 0xFFFF;
    static fvec recip(const fvec &a) { return _mm512_recip_ps(a); }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, bvec mask, const void *base) {
      _mm512_mask_prefetch_i32gather_ps(idx, mask, base, scale, _MM_HINT_T0);
    }
    template<int scale>
    static fvec gather(const fvec &from, bvec mask, const ivec &idx, const void *base) {
      return _mm512_mask_i32gather_ps(from, mask, idx, base, scale);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return _mm512_mask_blend_ps(mask, a, b);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return _mm512_fmadd_ps(a, b, c);
    }
    static fvec zero() {
      return _mm512_setzero_ps();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm512_cmpeq_ps_mask(a, b);
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm512_cmpnle_ps_mask(a, b);
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm512_cmple_ps_mask(a, b);
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm512_cmplt_ps_mask(a, b);
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      return _mm512_cmpneq_epi32_mask(a, b);
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return _mm512_cmplt_epi32_mask(a, b);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm512_invsqrt_ps(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      #if __INTEL_COMPILER+0 < 1500
      *reinterpret_cast<__m512 *>(cos) = _mm512_cos_ps(a);
      return _mm512_sin_ps(a);
      #else
      return _mm512_sincos_ps(reinterpret_cast<__m512 *>(cos), a);
      #endif
    }
    static fscal reduce_add(const fvec &a) {
      return _mm512_reduce_add_ps(a);
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return _mm512_mullo_epi32(a, b);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return _mm512_mask_add_epi32(src, mask, a, b);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      return _mm512_mask_i32gather_epi32(from, mask, idx, base, scale);
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return _mm512_mask_add_ps(src, mask, a, b);
    }
    static void store(void *at, const fvec &a) {
      _mm512_store_ps(at, a);
    }
    static void int_store(void *at, const ivec &a) {
      _mm512_store_epi32(at, a);
    }
    static void mask_store(int *at, const bvec &a) {
      for (int i = 0; i < 16; i++) {
        at[i] = (a >> i) & 1;
      }
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm512_min_ps(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return mask & (1 << at);
    }
    static bool mask_testz(const bvec &mask) {
      return mask == 0;
    }

    static bvec mask_enable_lower(int n) {
      return 0xFFFF >> (VL - n);
    }

    static ivec int_load_vl(const int *a) {
      return _mm512_load_epi32(a);
    }
    static void int_clear_arr(int *a) {
      _mm512_store_epi32(a, ivec(0));
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm512_store_epi32(tmp, a);
      for (int i = 0; i < 16; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      *x = gather<1>(*x, mask, idxs, &base->x);
      *y = gather<1>(*y, mask, idxs, &base->y);
      *z = gather<1>(*z, mask, idxs, &base->z);
      *w = int_gather<1>(*w, mask, idxs, &base->w);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char *>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char *>(base) +  4);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char *>(base) +  8);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char *>(base) + 12);
      *r4 = gather<4>(*r4, mask, idxs, reinterpret_cast<const char *>(base) + 16);
      *r5 = gather<4>(*r5, mask, idxs, reinterpret_cast<const char *>(base) + 20);
      *r6 = gather<4>(*r6, mask, idxs, reinterpret_cast<const char *>(base) + 24);
      *r7 = gather<4>(*r7, mask, idxs, reinterpret_cast<const char *>(base) + 28);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char *>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char *>(base) +  4);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char *>(base) +  8);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char *>(base) + 12);
    }
    // Additional routines needed for the implementation of mixed precision
    static fvec cvtdown(const vector_ops<double,KNC>::fvec &lo,
                        const vector_ops<double,KNC>::fvec &hi) {
      __m512 t1 = _mm512_cvtpd_pslo(lo);
      __m512 t2 = _mm512_cvtpd_pslo(hi);
      return _mm512_mask_shuffle_f32x4(_mm512_undefined_ps(), 0xFF00, t2, t2,
                                       0x4E);
    }
    static vector_ops<double,KNC>::fvec cvtup_lo(const fvec &a) {
      return _mm512_cvtpslo_pd(a);
    }
    static vector_ops<double,KNC>::fvec cvtup_hi(const fvec &a) {
      return _mm512_cvtpslo_pd(_mm512_shuffle_f32x4(a, a, 0x4E));
    }
    static void mask_cvtup(const bvec &a, vector_ops<double,KNC>::bvec *blo, vector_ops<double,KNC>::bvec *bhi) {
      *blo = a & 0xFF;
      *bhi = a >> 8;
    }
};
#endif

//////////////////////////////////////////////////////////////////////////////
// AVX/SSE
//////////////////////////////////////////////////////////////////////////////

#ifndef __MIC__
// class definitions for integer and masks for AVX
// Note that we have to lower a number of operations to SSE, notably comparison
// and integer operations.
// The gather operations are emulated.
struct ivec32x8 {
  __m256i vec;
  ivec32x8() {}
  ivec32x8(__m256i m) { vec = m; }
  ivec32x8(const int * a) {
    vec = _mm256_load_si256(reinterpret_cast<const __m256i *>(a));
  }
  explicit ivec32x8(int i) { vec = _mm256_set1_epi32(i); }
  operator __m256i() const { return vec; }
  friend ivec32x8 operator &(const ivec32x8 &a, const ivec32x8 &b) {
    return _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(b)));
  }
  friend ivec32x8 operator |(const ivec32x8 &a, const ivec32x8 &b) {
    return _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(b)));
  }
  friend ivec32x8 operator +(const ivec32x8 &a, const ivec32x8 &b) {
    __m128i alo = _mm256_castsi256_si128(a);
    __m128i ahi = _mm256_extractf128_si256(a, 1);
    __m128i blo = _mm256_castsi256_si128(b);
    __m128i bhi = _mm256_extractf128_si256(b, 1);
    __m128i rlo = _mm_add_epi32(alo, blo);
    __m128i rhi = _mm_add_epi32(ahi, bhi);
    return _mm256_setr_m128i(rlo, rhi);
  }
};

struct avx_bvec {
  __m256i vec;
  avx_bvec() {}
  avx_bvec(__m256i m) { vec = m; }
  explicit avx_bvec(int i) { vec = _mm256_set1_epi32(i); }
  operator __m256i() const { return vec; }
  operator F64vec4() const { return _mm256_castsi256_pd(vec); }
  operator F32vec8() const { return _mm256_castsi256_ps(vec); }
  operator ivec32x8() const { return vec; }
  friend avx_bvec operator &(const avx_bvec &a, const avx_bvec &b) {
    return _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(b)));
  }
  friend avx_bvec operator |(const avx_bvec &a, const avx_bvec &b) {
    return _mm256_castpd_si256(_mm256_or_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(b)));
  }
  friend avx_bvec operator ~(const avx_bvec &a) { return _mm256_castpd_si256(_mm256_andnot_pd(_mm256_castsi256_pd(a), _mm256_castsi256_pd(avx_bvec(0xFFFFFFFF)))); }
  avx_bvec& operator &=(const avx_bvec &a) { return *this = _mm256_castpd_si256(_mm256_and_pd(_mm256_castsi256_pd(vec), _mm256_castsi256_pd(a))); }
};

template<>
struct vector_ops<double, AVX> {
    static const int VL = 4;
    typedef double fscal;
    typedef F64vec4 fvec;
    typedef ivec32x8 ivec;
    typedef avx_bvec bvec;
    typedef double farr[4] __attribute__((aligned(32)));
    typedef int iarr[8] __attribute__((aligned(32)));
    static fvec recip(const fvec &a) {
      // newton-raphson
      fvec b = _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(a)));
      b = b + b - a * b * b;
      return  b + b - a * b * b;
    }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      farr result;
      farr src;
      iarr idxs;
      _mm256_store_si256(reinterpret_cast<__m256i*>(idxs), idx);
      _mm256_store_pd(reinterpret_cast<double*>(src), from);
      for (int i = 0; i < VL; i++) {
        result[i] = mask_test_at(mask, i)
            ? *reinterpret_cast<const double*>(reinterpret_cast<const char*>(base) + scale * idxs[2*i])
            : src[i];
      }
      return _mm256_load_pd(reinterpret_cast<double*>(result));
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      iarr i, m;
      int_store(i, idxs);
      mask_store(m, mask);
      __m256d a0 = m[0] ? _mm256_load_pd(reinterpret_cast<const double*>(&base[i[0]/32])) : _mm256_setzero_pd();
      __m256d a1 = m[2] ? _mm256_load_pd(reinterpret_cast<const double*>(&base[i[1]/32])) : _mm256_setzero_pd();
      __m256d b0 = _mm256_unpacklo_pd(a0, a1);
      __m256d b1 = _mm256_unpackhi_pd(a0, a1);
      __m256d a2 = m[4] ? _mm256_load_pd(reinterpret_cast<const double*>(&base[i[2]/32])) : _mm256_setzero_pd();
      __m256d a3 = m[6] ? _mm256_load_pd(reinterpret_cast<const double*>(&base[i[3]/32])) : _mm256_setzero_pd();
      __m256d b2 = _mm256_unpacklo_pd(a2, a3);
      __m256d b3 = _mm256_unpackhi_pd(a2, a3);
      __m256d c0 = _mm256_permute2f128_pd(b0, b2, 0x20);
      __m256d c1 = _mm256_permute2f128_pd(b1, b3, 0x20);
      __m256d c2 = _mm256_permute2f128_pd(b0, b2, 0x31);
      __m256d c3 = _mm256_permute2f128_pd(b1, b3, 0x31);
      *x = blend(mask, *x, c0);
      *y = blend(mask, *y, c1);
      *z = blend(mask, *z, c2);
      *w = int_blend(mask, *w, _mm256_castps_si256(_mm256_permute_ps(_mm256_castpd_ps(c3), 0xA0)));
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      fvec a = zero(), b = zero(), c = zero(), d = zero();
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 32, r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      iarr i, m;
      _mm256_store_si256(reinterpret_cast<__m256i*>(i), idxs);
      mask_store(m, mask);
      __m256d z = _mm256_setzero_pd();
      const char * m0 = reinterpret_cast<const char*>(base) + 4 * i[0];
      const char * m1 = reinterpret_cast<const char*>(base) + 4 * i[2];
      const char * m2 = reinterpret_cast<const char*>(base) + 4 * i[4];
      const char * m3 = reinterpret_cast<const char*>(base) + 4 * i[6];
      const double * e0 = reinterpret_cast<const double*>(m[0] == 0 ? reinterpret_cast<const char*>(&z) : m0);
      const double * e1 = reinterpret_cast<const double*>(m[2] == 0 ? reinterpret_cast<const char*>(&z) : m1);
      const double * e2 = reinterpret_cast<const double*>(m[4] == 0 ? reinterpret_cast<const char*>(&z) : m2);
      const double * e3 = reinterpret_cast<const double*>(m[6] == 0 ? reinterpret_cast<const char*>(&z) : m3);
      __m256d a0 = _mm256_load_pd(e0);
      __m256d a1 = _mm256_load_pd(e1);
      __m256d b0 = _mm256_unpacklo_pd(a0, a1);
      __m256d b1 = _mm256_unpackhi_pd(a0, a1);
      __m256d a2 = _mm256_load_pd(e2);
      __m256d a3 = _mm256_load_pd(e3);
      __m256d b2 = _mm256_unpacklo_pd(a2, a3);
      __m256d b3 = _mm256_unpackhi_pd(a2, a3);
      __m256d c0 = _mm256_permute2f128_pd(b0, b2, 0x20);
      __m256d c1 = _mm256_permute2f128_pd(b1, b3, 0x20);
      __m256d c2 = _mm256_permute2f128_pd(b0, b2, 0x31);
      __m256d c3 = _mm256_permute2f128_pd(b1, b3, 0x31);
      *r0 = blend(mask, *r0, c0);
      *r1 = blend(mask, *r1, c1);
      *r2 = blend(mask, *r2, c2);
      *r3 = blend(mask, *r3, c3);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return _mm256_setzero_pd();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_EQ_OQ));
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_NLE_US));
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_LE_OS));
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_LT_OS));
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      __m128i alo = _mm256_castsi256_si128(a);
      __m128i ahi = _mm256_extractf128_si256(a, 1);
      __m128i blo = _mm256_castsi256_si128(b);
      __m128i bhi = _mm256_extractf128_si256(b, 1);
      __m128i rlo = _mm_andnot_si128(_mm_cmpeq_epi32(alo, blo), _mm_cmpeq_epi8(alo, alo));
      __m128i rhi = _mm_andnot_si128(_mm_cmpeq_epi32(ahi, bhi), _mm_cmpeq_epi8(alo, alo));
      return _mm256_setr_m128i(rlo, rhi);
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      __m128i alo = _mm256_castsi256_si128(a);
      __m128i ahi = _mm256_extractf128_si256(a, 1);
      __m128i blo = _mm256_castsi256_si128(b);
      __m128i bhi = _mm256_extractf128_si256(b, 1);
      __m128i rlo = _mm_cmplt_epi32(alo, blo);
      __m128i rhi = _mm_cmplt_epi32(ahi, bhi);
      return _mm256_setr_m128i(rlo, rhi);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm256_invsqrt_pd(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      return _mm256_sincos_pd(reinterpret_cast<__m256d *>(cos), a);
    }
    static fscal reduce_add(const fvec &a) {
      __m256d t1 = _mm256_hadd_pd(a, a);
      __m128d t2 = _mm256_extractf128_pd(t1, 1);
      __m128d t3 = _mm256_castpd256_pd128(t1);
      return _mm_cvtsd_f64(_mm_add_pd(t2, t3));
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      __m128i alo = _mm256_castsi256_si128(a);
      __m128i ahi = _mm256_extractf128_si256(a, 1);
      __m128i blo = _mm256_castsi256_si128(b);
      __m128i bhi = _mm256_extractf128_si256(b, 1);
      __m128i rlo = _mm_mullo_epi32(alo, blo);
      __m128i rhi = _mm_mullo_epi32(ahi, bhi);
      return _mm256_setr_m128i(rlo, rhi);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      iarr result;
      iarr src;
      iarr idxs;
      _mm256_store_si256(reinterpret_cast<__m256i*>(idxs), idx);
      _mm256_store_si256(reinterpret_cast<__m256i*>(src), from);
      for (int i = 0; i < VL; i++) {
        int tmp;
        if (mask_test_at(mask, i)) {
          tmp = *reinterpret_cast<const int*>(reinterpret_cast<const char*>(base) + scale * idxs[2*i]);
        } else {
          tmp = src[2*i];
        }
        result[2*i] = tmp;
        result[2*i+1] = tmp;
      }
      return _mm256_load_si256(reinterpret_cast<__m256i*>(result));
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    static void store(void *at, const fvec &a) {
      _mm256_store_pd(reinterpret_cast<double*>(at), a);
    }
    static void int_store(int *at, const ivec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
      at[1] = at[2];
      at[2] = at[4];
      at[3] = at[6];
    }
    static void mask_store(int *at, const bvec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm256_min_pd(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return reinterpret_cast<const int*>(&mask)[2*at];
    }
    static bool mask_testz(const bvec &mask) {
      return _mm256_testz_si256(mask, mask);
    }
    static bvec mask_enable_lower(int n) {
      static const int base[8] __attribute__((aligned(64))) = {0,0,1,1,2,2,3,3};
      return int_cmplt(ivec(base), ivec(n));
    }
    static ivec int_load_vl(const int *a) {
      __m128i b = _mm_load_si128(reinterpret_cast<const __m128i*>(a));
      __m128i c = _mm_unpacklo_epi32(b, b);
      __m128i d = _mm_unpackhi_epi32(b, b);
      return _mm256_setr_m128i(c, d);
      // stolen from http://stackoverflow.com/questions/6687535/unexpected-result-from-avx-m256-unpack-ps-unpack-intrinsic
    }
    static void int_clear_arr(int *a) {
      _mm256_store_si256(reinterpret_cast<__m256i *>(a), ivec(0));
    }
    static bvec full_mask() {
      return bvec(0xFFFFFFFF);
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm256_store_si256(reinterpret_cast<__m256i *>(tmp), a);
      for (int i = 0; i < 8; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
};

template<>
struct vector_ops<float, AVX> {
    static const int VL = 8;
    static const int ALIGN = 32;
    typedef float fscal;
    typedef F32vec8 fvec;
    typedef ivec32x8 ivec;
    typedef avx_bvec bvec;
    typedef float farr[8] __attribute__((aligned(32)));
    typedef int iarr[8] __attribute__((aligned(32)));
    static fvec recip(const fvec &a) {
      fvec b = _mm256_rcp_ps(a);
      b = b + b - a * b * b;
      return  b + b - a * b * b;
    }

    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, bvec mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, bvec mask, const ivec &idx, const void *base) {
      farr result;
      farr src;
      iarr idxs;
      _mm256_store_si256(reinterpret_cast<__m256i*>(idxs), idx);
      _mm256_store_ps(reinterpret_cast<float*>(src), from);
      for (int i = 0; i < VL; i++) {
        result[i] = mask_test_at(mask, i)
            ? *reinterpret_cast<const float*>(reinterpret_cast<const char*>(base) + scale * idxs[i])
            : src[i];
      }
      return _mm256_load_ps(reinterpret_cast<float*>(result));
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      iarr i, m;
      int_store(i, idxs);
      mask_store(m, mask);
      farr zero_mem;
      store(zero_mem, zero());
      const float *e0 = m[0] ? reinterpret_cast<const float*>(&base[i[0]/16]) : zero_mem;
      const float *e1 = m[1] ? reinterpret_cast<const float*>(&base[i[1]/16]) : zero_mem;
      const float *e2 = m[2] ? reinterpret_cast<const float*>(&base[i[2]/16]) : zero_mem;
      const float *e3 = m[3] ? reinterpret_cast<const float*>(&base[i[3]/16]) : zero_mem;
      const float *e4 = m[4] ? reinterpret_cast<const float*>(&base[i[4]/16]) : zero_mem;
      const float *e5 = m[5] ? reinterpret_cast<const float*>(&base[i[5]/16]) : zero_mem;
      const float *e6 = m[6] ? reinterpret_cast<const float*>(&base[i[6]/16]) : zero_mem;
      const float *e7 = m[7] ? reinterpret_cast<const float*>(&base[i[7]/16]) : zero_mem;
      __m256 a0 = _mm256_loadu2_m128(e4, e0);
      __m256 a1 = _mm256_loadu2_m128(e5, e1);
      __m256 b0 = _mm256_unpacklo_ps(a0, a1);
      __m256 b1 = _mm256_unpackhi_ps(a0, a1);
      __m256 a2 = _mm256_loadu2_m128(e6, e2);
      __m256 a3 = _mm256_loadu2_m128(e7, e3);
      __m256 b2 = _mm256_unpacklo_ps(a2, a3);
      __m256 b3 = _mm256_unpackhi_ps(a2, a3);
      __m256 c0 = _mm256_shuffle_ps(b0, b2, 0x44);
      __m256 c1 = _mm256_shuffle_ps(b0, b2, 0xEE);
      __m256 c2 = _mm256_shuffle_ps(b1, b3, 0x44);
      __m256 c3 = _mm256_shuffle_ps(b1, b3, 0xEE);
      *x = blend(mask, *x, c0);
      *y = blend(mask, *y, c1);
      *z = blend(mask, *z, c2);
      *w = int_blend(mask, *w, _mm256_castps_si256(c3));
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      fvec a = zero(), b = zero(), c = zero(), d = zero();
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 16, r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      iarr i, m;
      int_store(i, idxs);
      mask_store(m, mask);
      farr zero_mem;
      store(zero_mem, zero());
      const float *e0 = m[0] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[0]]) : zero_mem;
      const float *e1 = m[1] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[1]]) : zero_mem;
      const float *e2 = m[2] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[2]]) : zero_mem;
      const float *e3 = m[3] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[3]]) : zero_mem;
      const float *e4 = m[4] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[4]]) : zero_mem;
      const float *e5 = m[5] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[5]]) : zero_mem;
      const float *e6 = m[6] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[6]]) : zero_mem;
      const float *e7 = m[7] ? reinterpret_cast<const float*>(&reinterpret_cast<const char *>(base)[4*i[7]]) : zero_mem;
      __m256 a0 = _mm256_loadu2_m128(e4, e0);
      __m256 a1 = _mm256_loadu2_m128(e5, e1);
      __m256 b0 = _mm256_unpacklo_ps(a0, a1);
      __m256 b1 = _mm256_unpackhi_ps(a0, a1);
      __m256 a2 = _mm256_loadu2_m128(e6, e2);
      __m256 a3 = _mm256_loadu2_m128(e7, e3);
      __m256 b2 = _mm256_unpacklo_ps(a2, a3);
      __m256 b3 = _mm256_unpackhi_ps(a2, a3);
      __m256 c0 = _mm256_shuffle_ps(b0, b2, 0x44);
      __m256 c1 = _mm256_shuffle_ps(b0, b2, 0xEE);
      __m256 c2 = _mm256_shuffle_ps(b1, b3, 0x44);
      __m256 c3 = _mm256_shuffle_ps(b1, b3, 0xEE);
      *r0 = blend(mask, *r0, c0);
      *r1 = blend(mask, *r1, c1);
      *r2 = blend(mask, *r2, c2);
      *r3 = blend(mask, *r3, c3);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return _mm256_setzero_ps();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_EQ_OQ));
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_NLE_US));
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_LE_OS));
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_LT_OS));
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      __m128i alo = _mm256_castsi256_si128(a);
      __m128i ahi = _mm256_extractf128_si256(a, 1);
      __m128i blo = _mm256_castsi256_si128(b);
      __m128i bhi = _mm256_extractf128_si256(b, 1);
      __m128i rlo = _mm_andnot_si128(_mm_cmpeq_epi32(alo, blo), _mm_cmpeq_epi8(alo, alo));
      __m128i rhi = _mm_andnot_si128(_mm_cmpeq_epi32(ahi, bhi), _mm_cmpeq_epi8(alo, alo));
      return _mm256_setr_m128i(rlo, rhi);
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      __m128i alo = _mm256_castsi256_si128(a);
      __m128i ahi = _mm256_extractf128_si256(a, 1);
      __m128i blo = _mm256_castsi256_si128(b);
      __m128i bhi = _mm256_extractf128_si256(b, 1);
      __m128i rlo = _mm_cmplt_epi32(alo, blo);
      __m128i rhi = _mm_cmplt_epi32(ahi, bhi);
      return _mm256_setr_m128i(rlo, rhi);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm256_invsqrt_ps(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      return _mm256_sincos_ps(reinterpret_cast<__m256 *>(cos), a);
    }
    static fscal reduce_add(const fvec &a) {
      __m256 t1 = _mm256_hadd_ps(a, a);
      __m128 t2 = _mm256_extractf128_ps(t1, 1);
      __m128 t3 = _mm256_castps256_ps128(t1);
      __m128 t4 = _mm_add_ps(t2, t3);
     __m128 t5 = _mm_permute_ps(t4, 0x1B); // 0x1B = reverse
      return _mm_cvtss_f32(_mm_add_ps(t4, t5));
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      __m128i alo = _mm256_castsi256_si128(a);
      __m128i ahi = _mm256_extractf128_si256(a, 1);
      __m128i blo = _mm256_castsi256_si128(b);
      __m128i bhi = _mm256_extractf128_si256(b, 1);
      __m128i rlo = _mm_mullo_epi32(alo, blo);
      __m128i rhi = _mm_mullo_epi32(ahi, bhi);
      return _mm256_setr_m128i(rlo, rhi);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      iarr result;
      iarr src;
      iarr idxs;
      _mm256_store_si256(reinterpret_cast<__m256i*>(idxs), idx);
      _mm256_store_si256(reinterpret_cast<__m256i*>(src), from);
      for (int i = 0; i < VL; i++) {
        result[i] = mask_test_at(mask, i)
            ? *reinterpret_cast<const int*>(reinterpret_cast<const char*>(base) + scale * idxs[i])
            : src[i];
      }
      return _mm256_load_si256(reinterpret_cast<__m256i*>(result));
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    static void store(void *at, const fvec &a) {
      _mm256_store_ps(reinterpret_cast<float*>(at), a);
    }
    static void int_store(void *at, const ivec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
    }
    static void mask_store(int *at, const bvec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm256_min_ps(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return reinterpret_cast<const int*>(&mask)[at];
    }
    static bool mask_testz(const bvec &mask) {
      return _mm256_testz_si256(mask, mask);
    }
    static bvec mask_enable_lower(int n) {
      static const int base[8] __attribute__((aligned(64))) = {0,1,2,3,4,5,6,7};
      return int_cmplt(ivec(base), ivec(n));
    }
    static ivec int_load_vl(const int *a) {
     return _mm256_load_si256(reinterpret_cast<const __m256i*>(a));
    }
    static void int_clear_arr(int *a) {
      _mm256_store_si256(reinterpret_cast<__m256i *>(a), ivec(0));
    }
    static bvec full_mask() {
      return bvec(0xFFFFFFFF);
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm256_store_si256(reinterpret_cast<__m256i *>(tmp), a);
      for (int i = 0; i < 16; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
    static fvec cvtdown(const vector_ops<double,AVX>::fvec &lo, const vector_ops<double,AVX>::fvec &hi) {
      __m128 t1 = _mm256_cvtpd_ps(lo);
      __m128 t2 = _mm256_cvtpd_ps(hi);
      return _mm256_setr_m128(t1, t2);
    }
    static vector_ops<double,AVX>::fvec cvtup_lo(const fvec &a) {
      return _mm256_cvtps_pd(_mm256_castps256_ps128(a));
    }
    static vector_ops<double,AVX>::fvec cvtup_hi(const fvec &a) {
      return _mm256_cvtps_pd(_mm256_extractf128_ps(a, 1)); // permute DCBA -> BADC
    }
    static void mask_cvtup(const bvec &a, vector_ops<double,AVX>::bvec *blo, vector_ops<double,AVX>::bvec *bhi) {
      __m256i t1 = _mm256_castps_si256(_mm256_unpacklo_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(a)));
      __m256i t2 = _mm256_castps_si256(_mm256_unpackhi_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(a)));
      *blo = _mm256_permute2f128_si256(t1, t2, 0x20);
      *bhi = _mm256_permute2f128_si256(t1, t2, 0x31);
    }
};

// AVX2
//

struct avx2_ivec32 {
  __m256i vec;
  avx2_ivec32() {}
  avx2_ivec32(__m256i m) { vec = m; }
  avx2_ivec32(const int * a) {
    vec = _mm256_load_si256(reinterpret_cast<const __m256i *>(a));
  }
  explicit avx2_ivec32(int i) { vec = _mm256_set1_epi32(i); }
  operator __m256i() const { return vec; }
  friend avx2_ivec32 operator &(const avx2_ivec32 &a, const avx2_ivec32 &b) {
    return _mm256_and_si256(a, b);
  }
  friend avx2_ivec32 operator |(const avx2_ivec32 &a, const avx2_ivec32 &b) {
    return _mm256_or_si256(a, b);
  }
  friend avx2_ivec32 operator +(const avx2_ivec32 &a, const avx2_ivec32 &b) {
    return _mm256_add_epi32(a, b);
  }
};

struct avx2_bvec {
  __m256i vec;
  avx2_bvec() {}
  avx2_bvec(__m256i m) { vec = m; }
  explicit avx2_bvec(int i) { vec = _mm256_set1_epi32(i); }
  operator __m256i() const { return vec; }
  operator __m256d() const { return _mm256_castsi256_pd(vec); }
  operator __m256() const { return _mm256_castsi256_ps(vec); }
  operator F64vec4() const { return _mm256_castsi256_pd(vec); }
  operator F32vec8() const { return _mm256_castsi256_ps(vec); }
  operator avx2_ivec32() const { return vec; }
  friend avx2_bvec operator &(const avx2_bvec &a, const avx2_bvec &b) {
    return _mm256_and_si256(a, b);
  }
  friend avx2_bvec operator |(const avx2_bvec &a, const avx2_bvec &b) {
    return _mm256_or_si256(a, b);
  }
  friend avx2_bvec operator ~(const avx2_bvec &a) {
    return _mm256_andnot_si256(a, avx2_bvec(0xFFFFFFFF));
  }
  avx2_bvec& operator &=(const avx2_bvec &a) { return *this = _mm256_and_si256(vec,a); }
};

template<>
struct vector_ops<double, AVX2> {
    static const int VL = 4;
    typedef double fscal;
    typedef F64vec4 fvec;
    typedef avx2_ivec32 ivec;
    typedef avx2_bvec bvec;
    typedef double farr[4] __attribute__((aligned(32)));
    typedef int iarr[8] __attribute__((aligned(32)));
    static fvec recip(const fvec &a) {
      // newton-raphson
      fvec b = _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(a)));
      b = b + b - a * b * b;
      return  b + b - a * b * b;
    }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      ivec idx0 = _mm256_shuffle_epi32(idx, 0xD8); // 11011000 ->3120
      ivec idx1 = _mm256_permute4x64_epi64(idx0, 0xD8);
      return _mm256_mask_i32gather_pd(from, static_cast<const double*>(base), _mm256_castsi256_si128(idx1), mask, scale);
    }
    template<class T>
    static void gather_x(const ivec &idx, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      ivec idx0 = _mm256_shuffle_epi32(idx, 0xD8); // 11011000 ->3120
      ivec idx1 = _mm256_permute4x64_epi64(idx0, 0xD8);
      *x = _mm256_mask_i32gather_pd(*x, &base->x, _mm256_castsi256_si128(idx1), mask, 1);
      *y = _mm256_mask_i32gather_pd(*y, &base->y, _mm256_castsi256_si128(idx1), mask, 1);
      *z = _mm256_mask_i32gather_pd(*z, &base->z, _mm256_castsi256_si128(idx1), mask, 1);
      *w = _mm256_mask_i32gather_epi32(*w, &base->w, idx, mask, 1);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      fvec a = zero(), b = zero(), c = zero(), d = zero();
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 32, r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idx, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      ivec idx0 = _mm256_shuffle_epi32(idx, 0xD8); // 11011000 ->3120
      ivec idx1 = _mm256_permute4x64_epi64(idx0, 0xD8);
      *r0 = _mm256_mask_i32gather_pd(*r0, static_cast<const double*>(base) + 0, _mm256_castsi256_si128(idx1), mask, 4);
      *r1 = _mm256_mask_i32gather_pd(*r1, static_cast<const double*>(base) + 1, _mm256_castsi256_si128(idx1), mask, 4);
      *r2 = _mm256_mask_i32gather_pd(*r2, static_cast<const double*>(base) + 2, _mm256_castsi256_si128(idx1), mask, 4);
      *r3 = _mm256_mask_i32gather_pd(*r3, static_cast<const double*>(base) + 3, _mm256_castsi256_si128(idx1), mask, 4);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return _mm256_setzero_pd();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_EQ_OQ));
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_NLE_US));
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_LE_OS));
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm256_castpd_si256(_mm256_cmp_pd(a, b, _CMP_LT_OS));
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      return ~bvec(_mm256_cmpeq_epi32(a, b));
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return _mm256_cmpgt_epi32(b, a);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm256_invsqrt_pd(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      return _mm256_sincos_pd(reinterpret_cast<__m256d *>(cos), a);
    }
    static fscal reduce_add(const fvec &a) {
      __m256d t1 = _mm256_hadd_pd(a, a);
      __m128d t2 = _mm256_extractf128_pd(t1, 1);
      __m128d t3 = _mm256_castpd256_pd128(t1);
      return _mm_cvtsd_f64(_mm_add_pd(t2, t3));
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return _mm256_mullo_epi32(a, b);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      return _mm256_mask_i32gather_epi32(from, static_cast<const int*>(base), idx, mask, scale);
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    static void store(void *at, const fvec &a) {
      _mm256_store_pd(reinterpret_cast<double*>(at), a);
    }
    static void int_store(int *at, const ivec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
      at[1] = at[2];
      at[2] = at[4];
      at[3] = at[6];
    }
    static void mask_store(int *at, const bvec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm256_min_pd(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return reinterpret_cast<const int*>(&mask)[2*at];
    }
    static bool mask_testz(const bvec &mask) {
      return _mm256_testz_si256(mask, mask);
    }
    static bvec mask_enable_lower(int n) {
      static const int base[8] __attribute__((aligned(64))) = {0,0,1,1,2,2,3,3};
      return int_cmplt(ivec(base), ivec(n));
    }
    static ivec int_load_vl(const int *a) {
      __m128i b = _mm_load_si128(reinterpret_cast<const __m128i*>(a));
      __m128i c = _mm_unpacklo_epi32(b, b);
      __m128i d = _mm_unpackhi_epi32(b, b);
      return _mm256_setr_m128i(c, d);
      // stolen from http://stackoverflow.com/questions/6687535/unexpected-result-from-avx-m256-unpack-ps-unpack-intrinsic
    }
    static void int_clear_arr(int *a) {
      _mm256_store_si256(reinterpret_cast<__m256i *>(a), ivec(0));
    }
    static bvec full_mask() {
      return bvec(0xFFFFFFFF);
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm256_store_si256(reinterpret_cast<__m256i *>(tmp), a);
      for (int i = 0; i < 8; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
};

template<>
struct vector_ops<float, AVX2> {
    static const int VL = 8;
    static const int ALIGN = 32;
    typedef float fscal;
    typedef F32vec8 fvec;
    typedef avx2_ivec32 ivec;
    typedef avx2_bvec bvec;
    typedef float farr[8] __attribute__((aligned(32)));
    typedef int iarr[8] __attribute__((aligned(32)));
    static fvec recip(const fvec &a) {
      fvec b = _mm256_rcp_ps(a);
      b = b + b - a * b * b;
      return  b + b - a * b * b;
    }

    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, bvec mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, bvec mask, const ivec &idx, const void *base) {
      return _mm256_mask_i32gather_ps(from, static_cast<const float*>(base), idx, mask, scale);
    }
    template<class T>
    static void gather_x(const ivec &idx, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      *x = _mm256_mask_i32gather_ps(*x, reinterpret_cast<const float*>(base) + 0, idx, mask, 1);
      *y = _mm256_mask_i32gather_ps(*y, reinterpret_cast<const float*>(base) + 1, idx, mask, 1);
      *z = _mm256_mask_i32gather_ps(*z, reinterpret_cast<const float*>(base) + 2, idx, mask, 1);
      *w = _mm256_mask_i32gather_epi32(*w, reinterpret_cast<const int*>(base) + 3, idx, mask, 1);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char *>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char *>(base) +  4);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char *>(base) +  8);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char *>(base) + 12);
      *r4 = gather<4>(*r4, mask, idxs, reinterpret_cast<const char *>(base) + 16);
      *r5 = gather<4>(*r5, mask, idxs, reinterpret_cast<const char *>(base) + 20);
      *r6 = gather<4>(*r6, mask, idxs, reinterpret_cast<const char *>(base) + 24);
      *r7 = gather<4>(*r7, mask, idxs, reinterpret_cast<const char *>(base) + 28);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char *>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char *>(base) +  4);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char *>(base) +  8);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char *>(base) + 12);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return (b & mask) | (a & ~ mask);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return _mm256_setzero_ps();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_EQ_OQ));
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_NLE_US));
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_LE_OS));
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm256_castps_si256(_mm256_cmp_ps(a, b, _CMP_LT_OS));
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      return ~bvec(_mm256_cmpeq_epi32(a, b));
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return _mm256_cmpgt_epi32(b, a);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm256_invsqrt_ps(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      return _mm256_sincos_ps(reinterpret_cast<__m256 *>(cos), a);
    }
    static fscal reduce_add(const fvec &a) {
      __m256 t1 = _mm256_hadd_ps(a, a);
      __m128 t2 = _mm256_extractf128_ps(t1, 1);
      __m128 t3 = _mm256_castps256_ps128(t1);
      __m128 t4 = _mm_add_ps(t2, t3);
     __m128 t5 = _mm_permute_ps(t4, 0x1B); // 0x1B = reverse
      return _mm_cvtss_f32(_mm_add_ps(t4, t5));
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return _mm256_mullo_epi32(a, b);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      return _mm256_mask_i32gather_epi32(from, static_cast<const int*>(base), idx, mask, scale);
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return ((a + b) & mask) | (src & ~ mask);
    }
    static void store(void *at, const fvec &a) {
      _mm256_store_ps(reinterpret_cast<float*>(at), a);
    }
    static void int_store(void *at, const ivec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
    }
    static void mask_store(int *at, const bvec &a) {
      _mm256_store_si256(reinterpret_cast<__m256i*>(at), a);
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm256_min_ps(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return reinterpret_cast<const int*>(&mask)[at];
    }
    static bool mask_testz(const bvec &mask) {
      return _mm256_testz_si256(mask, mask);
    }
    static bvec mask_enable_lower(int n) {
      static const int base[8] __attribute__((aligned(64))) = {0,1,2,3,4,5,6,7};
      return int_cmplt(ivec(base), ivec(n));
    }
    static ivec int_load_vl(const int *a) {
     return _mm256_load_si256(reinterpret_cast<const __m256i*>(a));
    }
    static void int_clear_arr(int *a) {
      _mm256_store_si256(reinterpret_cast<__m256i *>(a), ivec(0));
    }
    static bvec full_mask() {
      return bvec(0xFFFFFFFF);
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm256_store_si256(reinterpret_cast<__m256i *>(tmp), a);
      for (int i = 0; i < 16; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
    static fvec cvtdown(const vector_ops<double,AVX2>::fvec &lo, const vector_ops<double,AVX2>::fvec &hi) {
      __m128 t1 = _mm256_cvtpd_ps(lo);
      __m128 t2 = _mm256_cvtpd_ps(hi);
      return _mm256_setr_m128(t1, t2);
    }
    static vector_ops<double,AVX2>::fvec cvtup_lo(const fvec &a) {
      return _mm256_cvtps_pd(_mm256_castps256_ps128(a));
    }
    static vector_ops<double,AVX2>::fvec cvtup_hi(const fvec &a) {
      return _mm256_cvtps_pd(_mm256_extractf128_ps(a, 1)); // permute DCBA -> BADC
    }
    static void mask_cvtup(const bvec &a, vector_ops<double,AVX2>::bvec *blo, vector_ops<double,AVX2>::bvec *bhi) {
      __m256i t1 = _mm256_castps_si256(_mm256_unpacklo_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(a)));
      __m256i t2 = _mm256_castps_si256(_mm256_unpackhi_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(a)));
      *blo = _mm256_permute2f128_si256(t1, t2, 0x20);
      *bhi = _mm256_permute2f128_si256(t1, t2, 0x31);
    }
};



//////////////////////////////////////////////////////////////////////////////
// SSE
//////////////////////////////////////////////////////////////////////////////

#pragma pack(push,16)

struct ivec32x4 {
  __m128i vec;
  ivec32x4() {}
  ivec32x4(__m128i m) { vec = m; }
  ivec32x4(const int * a) {
    vec = _mm_load_si128(reinterpret_cast<const __m128i *>(a));
  }
  explicit ivec32x4(int i) { vec = _mm_set1_epi32(i); }
  operator __m128i() const { return vec; }
  friend ivec32x4 operator &(const ivec32x4 &a, const ivec32x4 &b) {
    return _mm_castpd_si128(_mm_and_pd(_mm_castsi128_pd(a), _mm_castsi128_pd(b)));
  }
  friend ivec32x4 operator |(const ivec32x4 &a, const ivec32x4 &b) {
    return _mm_castpd_si128(_mm_or_pd(_mm_castsi128_pd(a), _mm_castsi128_pd(b)));
  }
  friend ivec32x4 operator +(const ivec32x4 &a, const ivec32x4 &b) {
    return _mm_add_epi32(a, b);
  }
};

struct sse_bvecx4 {
  __m128i vec;
  sse_bvecx4() {}
  sse_bvecx4(__m128i m) { vec = m; }
  explicit sse_bvecx4(int i) { vec = _mm_set1_epi32(i); }
  operator __m128i() const { return vec; }
  operator F64vec2() const { return _mm_castsi128_pd(vec); }
  operator ivec32x4() const { return vec; }
  friend sse_bvecx4 operator &(const sse_bvecx4 &a, const sse_bvecx4 &b) {
    return _mm_castpd_si128(_mm_and_pd(_mm_castsi128_pd(a), _mm_castsi128_pd(b)));
  }
  friend sse_bvecx4 operator |(const sse_bvecx4 &a, const sse_bvecx4 &b) {
    return _mm_castpd_si128(_mm_or_pd(_mm_castsi128_pd(a), _mm_castsi128_pd(b)));
  }
  friend sse_bvecx4 operator ~(const sse_bvecx4 &a) { return _mm_castpd_si128(_mm_andnot_pd(_mm_castsi128_pd(a), _mm_castsi128_pd(sse_bvecx4(0xFFFFFFFF)))); }
  sse_bvecx4& operator &=(const sse_bvecx4 &a) { return *this = _mm_and_si128(vec,a); }
};


#pragma pack(pop)

template<>
struct vector_ops<double, SSE> {
    static const int VL = 2;
    typedef double fscal;
    typedef F64vec2 fvec;
    typedef ivec32x4 ivec;
    typedef sse_bvecx4 bvec;
    typedef double farr[2] __attribute__((aligned(16)));
    typedef int iarr[4] __attribute__((aligned(16)));
    static fvec recip(const fvec &a) {
      fvec b = _mm_cvtps_pd(_mm_rcp_ps(_mm_cvtpd_ps(a)));
      b = b + b - a * b * b;
      return  b + b - a * b * b;
    }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      fvec res = from;
      if (_mm_extract_epi32(mask, 0)) {
          res = _mm_loadl_pd(res, reinterpret_cast<const double*>(reinterpret_cast<const char*>(base) + scale*_mm_extract_epi32(idx, 0)));
      }
      if (_mm_extract_epi32(mask, 2)) {
          res = _mm_loadh_pd(res, reinterpret_cast<const double*>(reinterpret_cast<const char*>(base) + scale*_mm_extract_epi32(idx, 2)));
      }
      return res;
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      __m128d a0lo, a0hi, a1lo, a1hi;
      if (_mm_extract_epi32(mask, 0)) {
        a0lo = _mm_load_pd(reinterpret_cast<const double*>(&base[_mm_extract_epi32(idxs, 0)/32]));
        a0hi = _mm_load_pd(reinterpret_cast<const double*>(&base[_mm_extract_epi32(idxs, 0)/32].z));
      }
      if (_mm_extract_epi32(mask, 2)) {
        a1lo = _mm_load_pd(reinterpret_cast<const double*>(&base[_mm_extract_epi32(idxs, 2)/32]));
        a1hi = _mm_load_pd(reinterpret_cast<const double*>(&base[_mm_extract_epi32(idxs, 2)/32].z));
      }
      __m128d c0 = _mm_unpacklo_pd(a0lo, a1lo);
      __m128d c1 = _mm_unpackhi_pd(a0lo, a1lo);
      __m128d c2 = _mm_unpacklo_pd(a0hi, a1hi);
      __m128d c3 = _mm_unpackhi_pd(a0hi, a1hi);
      *x = blend(mask, *x, c0);
      *y = blend(mask, *y, c1);
      *z = blend(mask, *z, c2);
      *w = int_blend(mask, *w, _mm_shuffle_epi32(_mm_castpd_si128(c3), 0xA0));
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      fvec a = zero(), b = zero(), c = zero(), d = zero();
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 32, r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char*>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char*>(base) +  8);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char*>(base) + 16);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char*>(base) + 24);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return (b & _mm_castsi128_pd(mask)) | _mm_andnot_pd(_mm_castsi128_pd(mask), a);
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return (b & mask) | _mm_andnot_si128(mask, a);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return _mm_setzero_pd();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm_castpd_si128(_mm_cmpeq_pd(a, b));
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm_castpd_si128(_mm_cmpnle_pd(a, b));
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm_castpd_si128(_mm_cmple_pd(a, b));
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm_castpd_si128(_mm_cmplt_pd(a, b));
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      __m128i we = _mm_undefined_si128();
      return _mm_andnot_si128(_mm_cmpeq_epi32(a, b), _mm_cmpeq_epi8(we, we));
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return _mm_cmplt_epi32(a, b);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm_invsqrt_pd(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      return _mm_sincos_pd(reinterpret_cast<__m128d *>(cos), a);
    }
    static fscal reduce_add(const fvec &a) {
      __m128d t1 = _mm_shuffle_pd(a, a, 1); // reverse vector
      return _mm_cvtsd_f64(_mm_add_pd(t1, a));
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return _mm_mullo_epi32(a, b);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return ((a + b) & mask) | _mm_andnot_si128(mask, src);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      iarr result;
      iarr src;
      iarr idxs, m;
      _mm_store_si128(reinterpret_cast<__m128i*>(idxs), idx);
      _mm_store_si128(reinterpret_cast<__m128i*>(src), from);
      _mm_store_si128(reinterpret_cast<__m128i*>(m), mask);
      for (int i = 0; i < VL; i++) {
        int tmp;
        if (m[2*i]) {
          tmp = *reinterpret_cast<const int*>(reinterpret_cast<const char*>(base) + scale * idxs[2*i]);
        } else {
          tmp = src[2*i];
        }
        result[2*i] = tmp;
        result[2*i+1] = tmp;
      }
      return _mm_load_si128(reinterpret_cast<__m128i*>(result));
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return ((a + b) & _mm_castsi128_pd(mask)) | _mm_andnot_pd(_mm_castsi128_pd(mask), src);
    }
    static void store(void *at, const fvec &a) {
      _mm_store_pd(reinterpret_cast<double*>(at), a);
    }
    static void int_store(int *at, const ivec &a) {
      _mm_store_si128(reinterpret_cast<__m128i*>(at), a);
      at[1] = at[2];
    }
    static void mask_store(int *at, const bvec &a) {
      _mm_store_si128(reinterpret_cast<__m128i*>(at), a);
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm_min_pd(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return reinterpret_cast<const int*>(&mask)[2*at];
    }
    static bool mask_testz(const bvec &mask) {
      return _mm_testz_si128(mask, mask);
    }
    static bvec mask_enable_lower(int n) {
      static const int base[4] __attribute__((aligned(32))) = {0,0,1,1};
      return int_cmplt(ivec(base), ivec(n));
    }
    static ivec int_load_vl(const int *a) {
      __m128i b = _mm_load_si128(reinterpret_cast<const __m128i*>(a));
      return _mm_unpacklo_epi32(b, b);
    }
    static void int_clear_arr(int *a) {
      _mm_store_si128(reinterpret_cast<__m128i *>(a), ivec(0));
    }
    static bvec full_mask() {
      return bvec(0xFFFFFFFF);
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm_store_si128(reinterpret_cast<__m128i *>(tmp), a);
      for (int i = 0; i < 4; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
};

template<>
struct vector_ops<float, SSE> {
    static const int VL = 4;
    static const int ALIGN = 16;
    typedef float fscal;
    typedef F32vec4 fvec;
    typedef ivec32x4 ivec;
    typedef sse_bvecx4 bvec;
    typedef float farr[4] __attribute__((aligned(16)));
    typedef int iarr[4] __attribute__((aligned(16)));
    static fvec recip(const fvec &a) {
      fvec b = _mm_rcp_ps(a);
      b = b + b - a * b * b;
      return  b + b - a * b * b;
    }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      farr result;
      farr src;
      iarr idxs, m;
      mask_store(m, mask);
      _mm_store_si128(reinterpret_cast<__m128i*>(idxs), idx);
      _mm_store_ps(reinterpret_cast<float*>(src), from);
      for (int i = 0; i < VL; i++) {
        result[i] = m[i]
            ? *reinterpret_cast<const float*>(reinterpret_cast<const char*>(base) + scale * idxs[i])
            : src[i];
      }
      return _mm_load_ps(reinterpret_cast<float*>(result));
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      *x = gather<1>(*x, mask, idxs, &base->x);
      *y = gather<1>(*y, mask, idxs, &base->y);
      *z = gather<1>(*z, mask, idxs, &base->z);
      *w = int_gather<1>(*w, mask, idxs, &base->w);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      fvec a = zero(), b = zero(), c = zero(), d = zero();
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 16, r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char*>(base) +  0);
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char*>(base) +  4);
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char*>(base) +  8);
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char*>(base) + 12);
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return (b & _mm_castsi128_ps(mask)) | _mm_andnot_ps(_mm_castsi128_ps(mask), a);
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return (b & mask) | _mm_andnot_si128(mask, a);
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return _mm_setzero_ps();
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return _mm_castps_si128(_mm_cmpeq_ps(a, b));
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return _mm_castps_si128(_mm_cmpnle_ps(a, b));
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return _mm_castps_si128(_mm_cmple_ps(a, b));
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return _mm_castps_si128(_mm_cmplt_ps(a, b));
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      __m128i we = _mm_undefined_si128();
      return _mm_andnot_si128(_mm_cmpeq_epi32(a, b), _mm_cmpeq_epi8(we, we));
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return _mm_cmplt_epi32(a, b);
    }
    static fvec invsqrt(const fvec &a) {
      return _mm_invsqrt_ps(a);
    }
    static fvec sincos(fvec *cos, const fvec &a) {
      return _mm_sincos_ps(reinterpret_cast<__m128 *>(cos), a);
    }
    static fscal reduce_add(const fvec &a) {
      __m128 t1 = _mm_hadd_ps(a, a);
      __m128 t2 = _mm_shuffle_ps(t1, t1, 0x1B); // reverse
      return _mm_cvtss_f32(_mm_add_ps(t1, t2));
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return _mm_mullo_epi32(a, b);
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return ((a + b) & mask) | _mm_andnot_si128(mask, src);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      iarr result;
      iarr src;
      iarr idxs, m;
      _mm_store_si128(reinterpret_cast<__m128i*>(idxs), idx);
      _mm_store_si128(reinterpret_cast<__m128i*>(src), from);
      _mm_store_si128(reinterpret_cast<__m128i*>(m), mask);
      for (int i = 0; i < VL; i++) {
        int tmp;
        if (m[i]) {
          tmp = *reinterpret_cast<const int*>(reinterpret_cast<const char*>(base) + scale * idxs[i]);
        } else {
          tmp = src[i];
        }
        result[i] = tmp;
      }
      return _mm_load_si128(reinterpret_cast<__m128i*>(result));
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return ((a + b) & _mm_castsi128_ps(mask)) | _mm_andnot_ps(_mm_castsi128_ps(mask), src);
    }
    static void store(void *at, const fvec &a) {
      _mm_store_ps(reinterpret_cast<float*>(at), a);
    }
    static void int_store(int *at, const ivec &a) {
      _mm_store_si128(reinterpret_cast<__m128i*>(at), a);
    }
    static void mask_store(int *at, const bvec &a) {
      _mm_store_si128(reinterpret_cast<__m128i*>(at), a);
    }
    static fvec min(const fvec &a, const fvec &b) {
      return _mm_min_ps(a, b);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return reinterpret_cast<const int*>(&mask)[at];
    }
    static bool mask_testz(const bvec &mask) {
      return _mm_testz_si128(mask, mask);
    }
    static bvec mask_enable_lower(int n) {
      static const int base[4] __attribute__((aligned(32))) = {0,1,2,3};
      return int_cmplt(ivec(base), ivec(n));
    }
    static ivec int_load_vl(const int *a) {
      return _mm_load_si128(reinterpret_cast<const __m128i*>(a));
    }
    static void int_clear_arr(int *a) {
      _mm_store_si128(reinterpret_cast<__m128i *>(a), ivec(0));
    }
    static bvec full_mask() {
      return bvec(0xFFFFFFFF);
    }
    static void int_print(const ivec &a) {
      iarr tmp;
      _mm_store_si128(reinterpret_cast<__m128i *>(tmp), a);
      for (int i = 0; i < 4; i++) printf("%d ", tmp[i]);
      printf("\n");
    }
    static fvec cvtdown(const vector_ops<double,SSE>::fvec &lo, const vector_ops<double,SSE>::fvec &hi) {
      __m128 t1 = _mm_cvtpd_ps(lo);
      __m128 t2 = _mm_cvtpd_ps(hi);
      return _mm_shuffle_ps(t1, t2, 0x44); // t2[1]t2[0]t1[1]t1[0]
    }
    static vector_ops<double,SSE>::fvec cvtup_lo(const fvec &a) {
      return _mm_cvtps_pd(a);
    }
    static vector_ops<double,SSE>::fvec cvtup_hi(const fvec &a) {
      return _mm_cvtps_pd(_mm_shuffle_ps(a, a, 0x4E)); // permute DCBA -> BADC
    }
    static void mask_cvtup(const bvec &a, vector_ops<double,SSE>::bvec *blo, vector_ops<double,SSE>::bvec *bhi) {
      *blo = _mm_unpacklo_epi32(a, a);
      *bhi = _mm_unpackhi_epi32(a, a);
    }
};


#endif

// Scalar implementation
template<class flt_t>
struct vector_ops<flt_t, NONE> {
    static const int VL = 1;
    static const int ALIGN = 4;
    typedef flt_t fscal;
    typedef flt_t fvec;
    typedef int ivec;
    typedef int bvec;
    typedef flt_t farr[1];
    typedef int iarr[1];
    static fvec recip(const fvec &a) {
      return ((flt_t) 1.) / a;
    }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      return mask ? *reinterpret_cast<const flt_t*>(reinterpret_cast<const char*>(base) + scale * idx) : from;
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      *x = gather<1>(*x, mask, idxs, &base->x);
      *y = gather<1>(*y, mask, idxs, &base->y);
      *z = gather<1>(*z, mask, idxs, &base->z);
      *w = int_gather<1>(*w, mask, idxs, &base->w);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      fvec a = zero(), b = zero(), c = zero(), d = zero();
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 4 * sizeof(fscal), r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base,
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char*>(base) +  0 * sizeof(fscal));
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char*>(base) +  1 * sizeof(fscal));
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char*>(base) +  2 * sizeof(fscal));
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char*>(base) +  3 * sizeof(fscal));
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return mask ? b : a;
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return mask ? b : a;
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return 0.;
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return a == b;
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return !(a <= b);
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return a <= b;
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return a < b;
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      return a != b;
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return a < b;
    }
    static fvec invsqrt(const fvec &a) {
      return 1. / std::sqrt(a);
    }
    static fvec sincos(fvec *c, const fvec &a) {
      *c = cos(a);
      return sin(a);
    }
    static fscal reduce_add(const fvec &a) {
      return a;
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return a * b;
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return mask ? a + b : src;
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      return mask ? *reinterpret_cast<const int*>(reinterpret_cast<const char*>(base) + scale * idx) : from;
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return mask ? a + b : src;
    }
    static void store(void *at, const fvec &a) {
      *reinterpret_cast<flt_t*>(at) = a;
    }
    static void int_store(int *at, const ivec &a) {
      *reinterpret_cast<int*>(at) = a;
    }
    static void mask_store(int *at, const bvec &a) {
      *at = a;
    }
    static fvec min(const fvec &a, const fvec &b) {
      return a < b ? a : b;
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return mask;
    }
    static bool mask_testz(const bvec &mask) {
      return ! mask;
    }
    static bvec mask_enable_lower(int n) {
      return n > 0 ? true : false;
    }
    static ivec int_load_vl(const int *a) {
      return *a;
    }
    static void int_clear_arr(int *a) {
      *a = 0;
    }
    static bvec full_mask() {
      return true;
    }
    static void int_print(const ivec &a) {
    }
};

// Mixins to implement mixed precision and single/single and double/double
// This one is for single/single and double/double
template<class BASE_flt_t, CalculationMode BASE_mic>
struct AccumulatorSameMixin {
  typedef vector_ops<BASE_flt_t, BASE_mic> BASE;
  typedef typename BASE::fvec avec;
  typedef typename BASE::farr aarr;

  static avec acc_mask_add(const avec &src, const typename BASE::bvec &m, const avec &a, const typename BASE::fvec &b) {
    return BASE::mask_add(src, m, a, b);
  }

  static typename BASE::fscal acc_reduce_add(const avec &a) {
    return BASE::reduce_add(a);
  }

  static avec acc_zero() {
    return BASE::zero();
  }

  static void acc_store(aarr mem, const avec &a) {
    BASE::store(mem, a);
  }

};

// Mixed precision for cases where double vectors contain fewer elements
template<class BASE_flt_t, class HIGH_flt_t, CalculationMode mic>
struct AccumulatorTwiceMixin {
  typedef vector_ops<BASE_flt_t, mic> BASE;
  typedef vector_ops<HIGH_flt_t, mic> HIGH;

  struct avec_t {
    typename HIGH::fvec lo, hi;
    avec_t(const typename HIGH::fvec &alo, const typename HIGH::fvec &ahi) : lo(alo), hi(ahi) {}
    avec_t(const typename BASE::fvec &a) {
      lo = BASE::cvtup_lo(a);
      hi = BASE::cvtup_hi(a);
    }
    friend avec_t operator +(const avec_t &a, const avec_t &b) {
      return avec_t(a.lo + b.lo, a.hi + b.hi);
    }
    friend avec_t operator -(const avec_t &a, const avec_t &b) {
      return avec_t(a.lo - b.lo, a.hi - b.hi);
    }
    friend avec_t operator *(const avec_t &a, const avec_t &b) {
      return avec_t(a.lo * b.lo, a.hi * b.hi);
    }
    operator typename BASE::fvec() const {
      return BASE::cvtdown(lo, hi);
    }
  };

  typedef avec_t avec;
  typedef typename HIGH::fscal aarr[BASE::VL] __attribute__((aligned(BASE::ALIGN)));

  static avec acc_mask_add(const avec &src, const typename BASE::bvec &m, const avec &a, const typename BASE::fvec &b) {
    typename HIGH::fvec blo = BASE::cvtup_lo(b);
    typename HIGH::fvec bhi = BASE::cvtup_hi(b);
    typename HIGH::bvec mlo, mhi;
    BASE::mask_cvtup(m, &mlo, &mhi);
    return avec(HIGH::mask_add(src.lo, mlo, a.lo, blo), HIGH::mask_add(src.hi, mhi, a.hi, bhi));
  }

  static typename HIGH::fscal acc_reduce_add(const avec &a) {
    return HIGH::reduce_add(a.lo + a.hi);
  }

  static avec acc_zero() {
    return avec(HIGH::zero(), HIGH::zero());
  }

  static void acc_store(aarr mem, const avec &a) {
    HIGH::store(mem, a.lo);
    HIGH::store(mem + BASE::VL / 2, a.hi);
  }

};

// For cases where vector_ops<float,x>::VL == vector_ops<double,x>::VL

template<class BASE_flt_t, class HIGH_flt_t, CalculationMode mic>
struct AccumulatorTwiceMixinNone {
  typedef vector_ops<BASE_flt_t, mic> BASE;
  typedef vector_ops<HIGH_flt_t, mic> HIGH;

  typedef typename HIGH::fvec avec;
  typedef typename HIGH::fscal aarr[BASE::VL];

  static avec acc_mask_add(const avec &src, const typename BASE::bvec &m, const avec &a, const typename BASE::fvec &b) {
     return HIGH::mask_add(src, m, a, static_cast<typename HIGH::fvec>(b));
  }
  static typename HIGH::fscal acc_reduce_add(const avec &a) {
    return HIGH::reduce_add(a);
  }

  static avec acc_zero() {
    return HIGH::zero();
  }

  static void acc_store(aarr mem, const avec &a) {
    HIGH::store(mem, a);
  }

};

// This is the interfact that the user will see in the end.
template<class flt_t, class acc_t, CalculationMode mic>
struct vector_routines {};

template<CalculationMode mic>
struct vector_routines<double,double,mic> : public vector_ops<double, mic>, public AccumulatorSameMixin<double, mic> {};

template<CalculationMode mic>
struct vector_routines<float,float,mic> : public vector_ops<float, mic>, public AccumulatorSameMixin<float, mic> {};

template<CalculationMode mic>
struct vector_routines<float,double,mic> : public vector_ops<float, mic>, public AccumulatorTwiceMixin<float,double, mic> {};

// Specialize for scalar
template<>
struct vector_routines<float,double,NONE> : public vector_ops<float, NONE>, public AccumulatorTwiceMixinNone<float,double, NONE> {};

} // namespace lmp_intel
