// clang-format off
/* -*- c++ -*- -------------------------------------------------------------
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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------

Vector intrinsics are temporarily being used for the Stillinger-Weber
potential to allow for advanced features in the AVX512 instruction set to
be exploited on early hardware. We hope to see compiler improvements for
AVX512 that will eliminate this requirement, so it is not recommended to
develop code based on the intrinsics implementation. Please e-mail the
authors for more details.

------------------------------------------------------------------------- */

#ifndef INTEL_SIMD_H
#define INTEL_SIMD_H

#include "intel_preprocess.h"
#include "immintrin.h"

#ifdef __AVX512F__

#ifndef _MM_SCALE_1
#define _MM_SCALE_1 1
#define _MM_SCALE_2 2
#define _MM_SCALE_4 4
#define _MM_SCALE_8 8
#endif

namespace ip_simd {

  typedef __mmask16 SIMD_mask;

  inline bool any(const SIMD_mask &m) { return m != 0; }

  struct SIMD_int {
    __m512i v;
    SIMD_int() {}
    SIMD_int(const __m512i in) : v(in) {}
    inline int & operator[](const int i) { return ((int *)&(v))[i]; }
    inline const int & operator[](const int i) const
      { return ((int *)&(v))[i]; }
    operator __m512i() const { return v;}
  };

  struct SIMD256_int {
    __m256i v;
    SIMD256_int() {}
    SIMD256_int(const __m256i in) : v(in) {}
    SIMD256_int(const int in) : v(_mm256_set1_epi32(in)) {}
    inline int & operator[](const int i) { return ((int *)&(v))[i]; }
    inline const int & operator[](const int i) const
      { return ((int *)&(v))[i]; }
#ifdef __INTEL_LLVM_COMPILER
    inline SIMD256_int operator&=(const int i)
      { v=_mm256_and_epi32(v, _mm256_set1_epi32(i)); return *this; };
#else
    inline SIMD256_int operator&=(const int i)
      { v=_mm256_and_si256(v, _mm256_set1_epi32(i)); return *this; };
#endif
    inline SIMD256_int operator+=(const int i)
      { v=_mm256_add_epi32(v, _mm256_set1_epi32(i)); return *this; };
    operator __m256i() const { return v;}
  };

  struct SIMD_float {
    __m512 v;
    SIMD_float() {}
    SIMD_float(const __m512 in) : v(in) {}
    operator __m512() const { return v;}
  };

  struct SIMD_double {
    __m512d v;
    SIMD_double() {}
    SIMD_double(const __m512d in) : v(in) {}
    SIMD_double(const double in) { v=_mm512_set1_pd(in); }
    inline double & operator[](const int i) { return ((double *)&(v))[i]; }
    inline const double & operator[](const int i) const
      { return ((double *)&(v))[i]; }
    operator __m512d() const { return v;}

    SIMD_double & operator=(const double i)
      { _mm512_set1_pd(i); return *this; }
    SIMD_double &operator=(const SIMD_double &i)
      { v = i.v; return *this; }

    SIMD_double operator-() { return _mm512_xor_pd(v, _mm512_set1_pd(-0.0)); }
    SIMD_double & operator+=(const SIMD_double & two)
      { v = _mm512_add_pd(v, two.v); return *this; }
    SIMD_double & operator-=(const SIMD_double & two)
      { v = _mm512_sub_pd(v, two.v); return *this; }
    SIMD_double & operator*=(const SIMD_double & two)
      { v = _mm512_mul_pd(v, two.v); return *this; }
  };

  template<class flt_t>
  class SIMD_type {
  };

  template<>
  class SIMD_type<float> {
   public:
    typedef SIMD_float SIMD_vec;
    static inline int width() { return 16; }
  };

  template<>
  class SIMD_type<double> {
   public:
    typedef SIMD_double SIMD_vec;
    static inline int width() { return 8; }
  };

  template<class flt_t, class acc_t>
  class is_same {
   public:
    static const int value = 1;
  };

  template<>
  class is_same<float,double> {
   public:
    static const int value = 0;
  };

  // ------- Set Operations

  inline SIMD256_int SIMD256_set(const int l0, const int l1, const int l2,
                                 const int l3, const int l4, const int l5,
                                 const int l6, const int l7) {
    return _mm256_setr_epi32(l0,l1,l2,l3,l4,l5,l6,l7);
  }

  inline SIMD_int SIMD_set(const int l0, const int l1, const int l2,
                           const int l3, const int l4, const int l5,
                           const int l6, const int l7, const int l8,
                           const int l9, const int l10, const int l11,
                           const int l12, const int l13, const int l14,
                           const int l15) {
    return _mm512_setr_epi32(l0,l1,l2,l3,l4,l5,l6,l7,
                             l8,l9,l10,l11,l12,l13,l14,l15);
  }

  inline SIMD256_int SIMD256_set(const int l) {
    return _mm256_set1_epi32(l);
  }

  inline SIMD_int SIMD_set(const int l) {
    return _mm512_set1_epi32(l);
  }

  inline SIMD_float SIMD_set(const float l) {
    return _mm512_set1_ps(l);
  }

  inline SIMD_double SIMD_set(const double l) {
    return _mm512_set1_pd(l);
  }

  inline SIMD256_int SIMD256_count() {
    return SIMD256_set(0,1,2,3,4,5,6,7);
  }

  inline SIMD_int SIMD_zero_masked(const SIMD_mask &m, const SIMD_int &one) {
    return _mm512_maskz_mov_epi32(m, one);
  }

  inline SIMD_float SIMD_zero_masked(const SIMD_mask &m,
                                     const SIMD_float &one) {
    return _mm512_maskz_mov_ps(m, one);
  }

  inline SIMD_double SIMD_zero_masked(const SIMD_mask &m,
                                     const SIMD_double &one) {
    return _mm512_maskz_mov_pd(m, one);
  }

  inline SIMD_float SIMD_set(const SIMD_float &src, const SIMD_mask &m,
                             const SIMD_float &one) {
    return _mm512_mask_mov_ps(src,m,one);
  }

  inline SIMD_double SIMD_set(const SIMD_double &src, const SIMD_mask &m,
                              const SIMD_double &one) {
    return _mm512_mask_mov_pd(src,m,one);
  }

  // -------- Load Operations

  inline SIMD256_int SIMD_load(const SIMD256_int *p) {
    return _mm256_load_epi32((int *)p);
  }

  inline SIMD_int SIMD_load(const int *p) {
    return _mm512_load_epi32(p);
  }

  inline SIMD_float SIMD_load(const float *p) {
    return _mm512_load_ps(p);
  }

  inline SIMD_double SIMD_load(const double *p) {
    return _mm512_load_pd(p);
  }

  inline SIMD_double SIMD_load(const SIMD_double *p) {
    return _mm512_load_pd((double *)p);
  }

  inline SIMD_int SIMD_loadz(const SIMD_mask &m, const int *p) {
    return _mm512_maskz_load_epi32(m, p);
  }

  inline SIMD_float SIMD_loadz(const SIMD_mask &m, const float *p) {
    return _mm512_maskz_load_ps(m, p);
  }

  inline SIMD_double SIMD_loadz(const SIMD_mask &m, const double *p) {
    return _mm512_maskz_load_pd(m, p);
  }

  inline SIMD256_int SIMD_gather(const int *p, const SIMD256_int &i) {
    return _mm256_i32gather_epi32(p, i, _MM_SCALE_4);
  }

  inline SIMD_int SIMD_gather(const int *p, const SIMD_int &i) {
    return _mm512_i32gather_epi32(i, p, _MM_SCALE_4);
  }

  inline SIMD_float SIMD_gather(const float *p, const SIMD_int &i) {
    return _mm512_i32gather_ps(i, p, _MM_SCALE_4);
  }

  inline SIMD_double SIMD_gather(const double *p, const SIMD256_int &i) {
    return _mm512_i32gather_pd(i, p, _MM_SCALE_8);
  }

  inline SIMD_double SIMD_gather(const double *p, const SIMD_int &i) {
    return _mm512_i32gather_pd(_mm512_castsi512_si256(i), p, _MM_SCALE_8);
  }

  inline SIMD_int SIMD_gather(const SIMD_mask &m, const int *p,
                              const SIMD_int &i) {
    return _mm512_mask_i32gather_epi32(_mm512_undefined_epi32(), m, i, p,
                                       _MM_SCALE_4);
  }

  inline SIMD_float SIMD_gather(const SIMD_mask &m, const float *p,
                                const SIMD_int &i) {
    return _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, p,
                                    _MM_SCALE_4);
  }

  inline SIMD_double SIMD_gather(const SIMD_mask &m, const double *p,
                                 const SIMD_int &i) {
    return _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), p, _MM_SCALE_8);
  }

  inline SIMD_double SIMD_gather(const SIMD_mask &m, const double *p,
                                 const SIMD256_int &i) {
    return _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    i, p, _MM_SCALE_8);
  }

  template <typename T>
  inline SIMD_int SIMD_gatherz_offset(const SIMD_mask &m, const int *p,
                                      const SIMD_int &i) {
  }

  template <>
  inline SIMD_int SIMD_gatherz_offset<float>(const SIMD_mask &m, const int *p,
                                             const SIMD_int &i) {
    return _mm512_mask_i32gather_epi32( _mm512_set1_epi32(0), m, i, p,
                                       _MM_SCALE_4);
  }

  template <>
  inline SIMD_int SIMD_gatherz_offset<double>(const SIMD_mask &m, const int *p,
                                              const SIMD_int &i) {
    return _mm512_mask_i32gather_epi32( _mm512_set1_epi32(0), m, i, p,
                                       _MM_SCALE_8);
  }

  inline SIMD_int SIMD_gatherz(const SIMD_mask &m, const int *p,
                               const SIMD_int &i) {
    return _mm512_mask_i32gather_epi32( _mm512_set1_epi32(0), m, i, p,
                                    _MM_SCALE_4);
  }

  inline SIMD_float SIMD_gatherz(const SIMD_mask &m, const float *p,
                                 const SIMD_int &i) {
    return _mm512_mask_i32gather_ps( _mm512_set1_ps((float)0), m, i, p,
                                    _MM_SCALE_4);
  }

  inline SIMD_double SIMD_gatherz(const SIMD_mask &m, const double *p,
                                  const SIMD_int &i) {
    return _mm512_mask_i32gather_pd( _mm512_set1_pd(0.0), m,
                                     _mm512_castsi512_si256(i),p, _MM_SCALE_8);
  }

  // ------- Store Operations

  inline void SIMD_store(int *p, const SIMD_int &one) {
    return _mm512_store_epi32(p,one);
  }

  inline void SIMD_store(float *p, const SIMD_float &one) {
    return _mm512_store_ps(p,one);
  }

  inline void SIMD_store(double *p, const SIMD_double &one) {
    return _mm512_store_pd(p,one);
  }

  inline void SIMD_store(SIMD_double *p, const SIMD_double &one) {
    return _mm512_store_pd((double *)p,one);
  }

  inline void SIMD_scatter(const SIMD_mask &m, int *p,
                           const SIMD256_int &i, const SIMD256_int &vec) {
    _mm256_mask_i32scatter_epi32(p, m, i, vec, _MM_SCALE_4);
  }

  inline void SIMD_scatter(const SIMD_mask &m, int *p,
                           const SIMD_int &i, const SIMD_int &vec) {
    _mm512_mask_i32scatter_epi32(p, m, i, vec, _MM_SCALE_4);
  }

  inline void SIMD_scatter(const SIMD_mask &m, float *p,
                           const SIMD_int &i, const SIMD_float &vec) {
    _mm512_mask_i32scatter_ps(p, m, i, vec, _MM_SCALE_4);
  }

  inline void SIMD_scatter(const SIMD_mask &m, double *p,
                           const SIMD_int &i, const SIMD_double &vec) {
    _mm512_mask_i32scatter_pd(p, m, _mm512_castsi512_si256(i), vec,
                              _MM_SCALE_8);
  }

  inline void SIMD_scatter(const SIMD_mask &m, double *p,
                           const SIMD256_int &i, const SIMD_double &vec) {
    _mm512_mask_i32scatter_pd(p, m, i, vec, _MM_SCALE_8);
  }

  inline void SIMD_scatter(double *p,
                           const SIMD256_int &i, const SIMD_double &vec) {
    _mm512_i32scatter_pd(p, i, vec, _MM_SCALE_8);
  }

  // ------- Arithmetic Operations

  inline SIMD256_int operator+(const SIMD256_int &one, const SIMD256_int &two) {
    return _mm256_add_epi32(one,two);
  }

  inline SIMD_int operator+(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_add_epi32(one,two);
  }

  inline SIMD_float operator+(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_add_ps(one,two);
  }

  inline SIMD_double operator+(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_add_pd(one,two);
  }

  inline SIMD_int operator+(const SIMD_int &one, const int two) {
    return _mm512_add_epi32(one,SIMD_set(two));
  }

  inline SIMD256_int operator+(const SIMD256_int &one, const int two) {
    return _mm256_add_epi32(one,SIMD256_set(two));
  }

  inline SIMD_float operator+(const SIMD_float &one, const float two) {
    return _mm512_add_ps(one,SIMD_set(two));
  }

  inline SIMD_double operator+(const SIMD_double &one, const double two) {
    return _mm512_add_pd(one,SIMD_set(two));
  }

  inline SIMD_int SIMD_add(const SIMD_mask &m,
                           const SIMD_int &one, const int two) {
    return _mm512_mask_add_epi32(one,m,one,SIMD_set(two));
  }

  inline SIMD256_int SIMD_add(const SIMD_mask &m,
                           const SIMD256_int &one, const int two) {
    return _mm256_mask_add_epi32(one,m,one,SIMD256_set(two));
  }

  inline SIMD_float SIMD_add(const SIMD_mask &m,
                             const SIMD_float &one, const float two) {
    return _mm512_mask_add_ps(one,m,one,SIMD_set(two));
  }

  inline SIMD_double SIMD_add(const SIMD_mask &m,
                              const SIMD_double &one, const double two) {
    return _mm512_mask_add_pd(one,m,one,SIMD_set(two));
  }

  inline SIMD_double SIMD_add(const SIMD_mask &m,
                              const SIMD_double &one, const SIMD_double &two) {
    return _mm512_mask_add_pd(one,m,one,two);
  }

  inline SIMD_int SIMD_add(const SIMD_int &s, const SIMD_mask &m,
                           const SIMD_int &one, const SIMD_int &two) {
    return _mm512_mask_add_epi32(s,m,one,two);
  }

  inline SIMD_float SIMD_add(const SIMD_float &s, const SIMD_mask &m,
                             const SIMD_float &one, const SIMD_float &two) {
    return _mm512_mask_add_ps(s,m,one,two);
  }

  inline SIMD_double SIMD_add(const SIMD_double &s, const SIMD_mask &m,
                              const SIMD_double &one, const SIMD_double &two) {
    return _mm512_mask_add_pd(s,m,one,two);
  }

  inline SIMD_int SIMD_sub(const SIMD_int &s, const SIMD_mask &m,
                           const SIMD_int &one, const SIMD_int &two) {
    return _mm512_mask_sub_epi32(s,m,one,two);
  }

  inline SIMD_float SIMD_sub(const SIMD_float &s, const SIMD_mask &m,
                             const SIMD_float &one, const SIMD_float &two) {
    return _mm512_mask_sub_ps(s,m,one,two);
  }

  inline SIMD_double SIMD_sub(const SIMD_double &s, const SIMD_mask &m,
                              const SIMD_double &one, const SIMD_double &two) {
    return _mm512_mask_sub_pd(s,m,one,two);
  }

  inline SIMD_int operator-(const SIMD_int &one) {
    return _mm512_sub_epi32(SIMD_set((int)0),one);
  }

  inline SIMD_float operator-(const SIMD_float &one) {
    return _mm512_sub_ps(SIMD_set((float)0),one);
  }

  inline SIMD_double operator-(const SIMD_double &one) {
    return _mm512_sub_pd(SIMD_set((double)0),one);
  }

  inline SIMD_int operator-(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_sub_epi32(one,two);
  }

  inline SIMD_float operator-(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_sub_ps(one,two);
  }

  inline SIMD_double operator-(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_sub_pd(one,two);
  }

  inline SIMD_int operator-(const SIMD_int &one, const int two) {
    return _mm512_sub_epi32(one,SIMD_set(two));
  }

  inline SIMD_float operator-(const SIMD_float &one, const float two) {
    return _mm512_sub_ps(one,SIMD_set(two));
  }

  inline SIMD_double operator-(const SIMD_double &one, const double two) {
    return _mm512_sub_pd(one,SIMD_set(two));
  }

  inline SIMD_int operator*(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_mullo_epi32(one,two);
  }

  inline SIMD_float operator*(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_mul_ps(one,two);
  }

  inline SIMD_double operator*(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_mul_pd(one,two);
  }

  inline SIMD256_int operator*(const SIMD256_int &one, const int two) {
    return _mm256_mullo_epi32(one,SIMD256_set(two));
  }

  inline SIMD_int operator*(const SIMD_int &one, const int two) {
    return _mm512_mullo_epi32(one,SIMD_set(two));
  }

  inline SIMD_float operator*(const SIMD_float &one, const float two) {
    return _mm512_mul_ps(one,SIMD_set(two));
  }

  inline SIMD_double operator*(const SIMD_double &one, const double two) {
    return _mm512_mul_pd(one,SIMD_set(two));
  }

  inline SIMD_float operator/(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_div_ps(one,two);
  }

  inline SIMD_double operator/(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_div_pd(one,two);
  }

  inline SIMD_float SIMD_fma(const SIMD_float &one, const SIMD_float &two,
                             const SIMD_float &three) {
    return _mm512_fmadd_ps(one,two,three);
  }

  inline SIMD_double SIMD_fma(const SIMD_double &one, const SIMD_double &two,
                              const SIMD_double &three) {
    return _mm512_fmadd_pd(one,two,three);
  }

  inline SIMD_double SIMD_fma(const SIMD_mask m, const SIMD_double &one,
                              const SIMD_double &two,
                              const SIMD_double &three) {
    return _mm512_mask3_fmadd_pd(one,two,three,m);
  }

  inline SIMD_float SIMD_fms(const SIMD_float &one, const SIMD_float &two,
                             const SIMD_float &three) {
    return _mm512_fmsub_ps(one,two,three);
  }

  inline SIMD_double SIMD_fms(const SIMD_double &one, const SIMD_double &two,
                              const SIMD_double &three) {
    return _mm512_fmsub_pd(one,two,three);
  }

  // ------- SVML operations

  inline SIMD_float SIMD_rcp(const SIMD_float &one) {
    #ifdef __AVX512ER__
    return _mm512_rcp28_ps(one);
    #else
    return _mm512_recip_ps(one);
    #endif
  }

  inline SIMD_double SIMD_rcp(const SIMD_double &one) {
    #ifdef __AVX512ER__
    return _mm512_rcp28_pd(one);
    #else
    return _mm512_recip_pd(one);
    #endif
  }

  inline SIMD_float SIMD_rcpz(const SIMD_mask &m, const SIMD_float &one) {
    #ifdef __AVX512ER__
    return _mm512_maskz_rcp28_ps(m, one);
    #else
    return _mm512_mask_recip_ps(_mm512_set1_ps(0), m, one);
    #endif
  }

  inline SIMD_double SIMD_rcpz(const SIMD_mask &m, const SIMD_double &one) {
    #ifdef __AVX512ER__
    return _mm512_maskz_rcp28_pd(m, one);
    #else
    return _mm512_mask_recip_pd(_mm512_set1_pd(0), m, one);
    #endif
  }

  inline SIMD_float SIMD_sqrt(const SIMD_float &one) {
    return _mm512_sqrt_ps(one);
  }

  inline SIMD_double SIMD_sqrt(const SIMD_double &one) {
    return _mm512_sqrt_pd(one);
  }

  inline SIMD_float SIMD_invsqrt(const SIMD_float &one) {
    #ifdef __AVX512ER__
    return _mm512_rsqrt28_ps(one);
    #else
    return _mm512_invsqrt_ps(one);
    #endif
  }

  inline SIMD_double SIMD_invsqrt(const SIMD_double &one) {
    #ifdef __AVX512ER__
    return _mm512_rsqrt28_pd(one);
    #else
    return _mm512_invsqrt_pd(one);
    #endif
  }

  inline SIMD_float SIMD_pow(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_pow_ps(one, two);
  }

  inline SIMD_double SIMD_pow(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_pow_pd(one, two);
  }

  inline SIMD_double SIMD_pow(const SIMD_double &one, const double two) {
    return _mm512_pow_pd(one, SIMD_set(two));
  }

  inline SIMD_float SIMD_exp(const SIMD_float &one) {
    return _mm512_exp_ps(one);
  }

  inline SIMD_double SIMD_exp(const SIMD_double &one) {
    return _mm512_exp_pd(one);
  }

  inline SIMD_double SIMD_cos(const SIMD_double &one) {
    return _mm512_cos_pd(one);
  }

  inline SIMD_double SIMD_sin(const SIMD_double &one) {
    return _mm512_sin_pd(one);
  }

  inline SIMD_double SIMD_tan(const SIMD_double &one) {
    return _mm512_tan_pd(one);
  }

  // ------- Comparison operations

  inline SIMD_mask SIMD_lt(SIMD_mask m, const SIMD_int &one,
                           const SIMD_int &two) {
    return _mm512_mask_cmplt_epi32_mask(m, one, two);
  }

  inline SIMD_mask SIMD_lt(SIMD_mask m, const SIMD_float &one,
                           const SIMD_float &two) {
    return _mm512_mask_cmplt_ps_mask(m, one, two);
  }

  inline SIMD_mask SIMD_lt(SIMD_mask m, const SIMD_double &one,
                           const SIMD_double &two) {
    return _mm512_mask_cmplt_pd_mask(m, one, two);
  }

  inline SIMD_mask SIMD_lt(SIMD_mask m, const int one,
                           const SIMD_int &two) {
    return _mm512_mask_cmplt_epi32_mask(m, SIMD_set(one), two);
  }

  inline SIMD_mask SIMD_lt(SIMD_mask m, const float one,
                           const SIMD_float &two) {
    return _mm512_mask_cmplt_ps_mask(m, SIMD_set(one), two);
  }

  inline SIMD_mask SIMD_lt(SIMD_mask m, const double one,
                           const SIMD_double &two) {
    return _mm512_mask_cmplt_pd_mask(m, SIMD_set(one), two);
  }

  inline SIMD_mask operator<(const SIMD256_int &one, const SIMD256_int &two) {
    return _mm256_cmplt_epi32_mask(one,two);
  }

  inline SIMD_mask operator<(const int one, const SIMD256_int &two) {
    return _mm256_cmplt_epi32_mask(SIMD256_set(one),two);
  }

  inline SIMD_mask operator<(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_cmplt_epi32_mask(one,two);
  }

  inline SIMD_mask operator<(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_cmplt_ps_mask(one,two);
  }

  inline SIMD_mask operator<(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_cmplt_pd_mask(one,two);
  }

  inline SIMD_mask operator<(const SIMD_int &one, const int two) {
    return _mm512_cmplt_epi32_mask(one,SIMD_set(two));
  }

  inline SIMD_mask operator<(const SIMD_float &one, const float two) {
    return _mm512_cmplt_ps_mask(one,SIMD_set(two));
  }

  inline SIMD_mask operator<(const SIMD_double &one, const double two) {
    return _mm512_cmplt_pd_mask(one,SIMD_set(two));
  }

  inline SIMD_mask operator<(const int one, const SIMD_int &two) {
    return _mm512_cmplt_epi32_mask(SIMD_set(one),two);
  }

  inline SIMD_mask operator<(const float one, const SIMD_float &two) {
    return _mm512_cmplt_ps_mask(SIMD_set(one),two);
  }

  inline SIMD_mask operator<(const double one, const SIMD_double &two) {
    return _mm512_cmplt_pd_mask(SIMD_set(one),two);
  }

  inline SIMD_mask operator<=(const int one, const SIMD_int &two) {
    return _mm512_cmple_epi32_mask(SIMD_set(one), two);
  }

  inline SIMD_mask operator<=(const float one, const SIMD_float &two) {
    return _mm512_cmple_ps_mask(SIMD_set(one), two);
  }

  inline SIMD_mask operator<=(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_cmple_pd_mask(one, two);
  }

  inline SIMD_mask operator<=(const double one, const SIMD_double &two) {
    return _mm512_cmple_pd_mask(SIMD_set(one), two);
  }

  inline SIMD_mask operator>(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_cmpgt_epi32_mask(one,two);
  }

  inline SIMD_mask operator>(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_cmplt_ps_mask(two,one);
  }

  inline SIMD_mask operator>(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_cmplt_pd_mask(two,one);
  }

  inline SIMD_mask operator>(const SIMD_double &one, const double two) {
    return _mm512_cmplt_pd_mask(SIMD_set(two),one);
  }

  inline SIMD_mask operator==(const SIMD256_int &one, const int two) {
    return _mm256_cmpeq_epi32_mask(one,_mm256_set1_epi32(two));
  }

  inline SIMD_mask operator==(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_cmpeq_epi32_mask(one,two);
  }

  inline SIMD_mask operator==(const SIMD_float &one, const SIMD_float &two) {
    return _mm512_cmpeq_ps_mask(one,two);
  }

  inline SIMD_mask operator==(const SIMD_double &one, const SIMD_double &two) {
    return _mm512_cmpeq_pd_mask(one,two);
  }

  // ------- Typecast operations

  inline void SIMD_cast(const SIMD_int &one, SIMD_float &two) {
    two = _mm512_cvtepi32_ps(one);
  }

  inline void SIMD_cast(const SIMD_int &one, SIMD_double &two) {
    two = _mm512_cvtepi32lo_pd(one);
  }

  // ------- Reduction operations

  inline int SIMD_max(const SIMD_int &i) {
    return _mm512_reduce_max_epi32(i);
  }

  inline float SIMD_max(const SIMD_float &i) {
    return _mm512_reduce_max_ps(i);
  }

  inline double SIMD_max(const SIMD_double &i) {
    return _mm512_reduce_max_pd(i);
  }

  inline int SIMD_sum(const SIMD_int &i) {
    return _mm512_reduce_add_epi32(i);
  }

  inline float SIMD_sum(const SIMD_float &i) {
    return _mm512_reduce_add_ps(i);
  }

  inline double SIMD_sum(const SIMD_double &i) {
    return _mm512_reduce_add_pd(i);
  }

  // i indices should be positive
  inline void SIMD_conflict_pi_reduce1(const SIMD_mask &m, const SIMD_int &i,
                                       SIMD_float &v1) {
    SIMD_int jc = _mm512_mask_mov_epi32(_mm512_set1_epi32(-1), m, i);
    SIMD_int cd = _mm512_maskz_conflict_epi32(m, jc);
    SIMD_mask todo_mask = _mm512_test_epi32_mask(cd, _mm512_set1_epi32(-1));
    if (todo_mask) {
      SIMD_int lz  = _mm512_lzcnt_epi32(cd);
      SIMD_int lid = _mm512_sub_epi32(_mm512_set1_epi32(31),
                                      _mm512_lzcnt_epi32(cd));

      while(todo_mask) {
        SIMD_int todo_bcast = _mm512_broadcastmw_epi32(todo_mask);
        SIMD_mask now_mask = _mm512_mask_testn_epi32_mask(todo_mask, cd,
                                                          todo_bcast);
        SIMD_float am_perm;
        am_perm = _mm512_mask_permutexvar_ps(_mm512_undefined_ps(),
                                             now_mask, lid, v1);
        v1 = _mm512_mask_add_ps(v1, now_mask, v1, am_perm);
        todo_mask = _mm512_kxor(todo_mask, now_mask);
      }
    }
  }

  // i indices should be positive
  inline void SIMD_conflict_pi_reduce1(const SIMD_mask &m, const SIMD_int &i,
                                       SIMD_double &v1) {
    SIMD_int jc = _mm512_mask_mov_epi32(_mm512_set1_epi32(-1), m, i);
    SIMD_int cd = _mm512_maskz_conflict_epi32(m, jc);
    SIMD_mask todo_mask = _mm512_test_epi32_mask(cd, _mm512_set1_epi32(-1));
    if (todo_mask) {
      SIMD_int lz  = _mm512_lzcnt_epi32(cd);
      SIMD_int lid = _mm512_sub_epi32(_mm512_set1_epi32(31),
                                      _mm512_lzcnt_epi32(cd));
      lid = _mm512_cvtepi32_epi64(_mm512_castsi512_si256(lid));

      while(todo_mask) {
        SIMD_int todo_bcast = _mm512_broadcastmw_epi32(todo_mask);
        SIMD_mask now_mask = _mm512_mask_testn_epi32_mask(todo_mask, cd,
                                                          todo_bcast);
        SIMD_double am_perm;
        am_perm = _mm512_mask_permutexvar_pd(_mm512_undefined_pd(),
                                             now_mask, lid, v1);
        v1 = _mm512_mask_add_pd(v1, now_mask, v1, am_perm);
        todo_mask = _mm512_kxor(todo_mask, now_mask);
      }
    }
  }

  // i indices should be positive
  inline void SIMD_conflict_pi_reduce3(const SIMD_mask &m, const SIMD_int &i,
                                       SIMD_float &v1, SIMD_float &v2,
                                       SIMD_float &v3) {
    SIMD_int jc = _mm512_mask_mov_epi32(_mm512_set1_epi32(-1), m, i);
    SIMD_int cd = _mm512_maskz_conflict_epi32(m, jc);
    SIMD_mask todo_mask = _mm512_test_epi32_mask(cd, _mm512_set1_epi32(-1));
    if (todo_mask) {
      SIMD_int lz  = _mm512_lzcnt_epi32(cd);
      SIMD_int lid = _mm512_sub_epi32(_mm512_set1_epi32(31),
                                      _mm512_lzcnt_epi32(cd));

      while(todo_mask) {
        SIMD_int todo_bcast = _mm512_broadcastmw_epi32(todo_mask);
        SIMD_mask now_mask = _mm512_mask_testn_epi32_mask(todo_mask, cd,
                                                          todo_bcast);
        SIMD_float am_perm;
        am_perm = _mm512_mask_permutexvar_ps(_mm512_undefined_ps(),
                                             now_mask, lid, v1);
        v1 = _mm512_mask_add_ps(v1, now_mask, v1, am_perm);
        am_perm = _mm512_mask_permutexvar_ps(_mm512_undefined_ps(),
                                             now_mask, lid, v2);
        v2 = _mm512_mask_add_ps(v2, now_mask, v2, am_perm);
        am_perm = _mm512_mask_permutexvar_ps(_mm512_undefined_ps(),
                                             now_mask, lid, v3);
        v3 = _mm512_mask_add_ps(v3, now_mask, v3, am_perm);
        todo_mask = _mm512_kxor(todo_mask, now_mask);
      }
    }
  }

  // i indices should be positive
  inline void SIMD_conflict_pi_reduce3(const SIMD_mask &m, const SIMD_int &i,
                                       SIMD_double &v1, SIMD_double &v2,
                                       SIMD_double &v3) {
    SIMD_int jc = _mm512_mask_mov_epi32(_mm512_set1_epi32(-1), m, i);
    SIMD_int cd = _mm512_maskz_conflict_epi32(m, jc);
    SIMD_mask todo_mask = _mm512_test_epi32_mask(cd, _mm512_set1_epi32(-1));
    if (todo_mask) {
      SIMD_int lz  = _mm512_lzcnt_epi32(cd);
      SIMD_int lid = _mm512_sub_epi32(_mm512_set1_epi32(31),
                                      _mm512_lzcnt_epi32(cd));
      lid = _mm512_cvtepi32_epi64(_mm512_castsi512_si256(lid));

      while(todo_mask) {
        SIMD_int todo_bcast = _mm512_broadcastmw_epi32(todo_mask);
        SIMD_mask now_mask = _mm512_mask_testn_epi32_mask(todo_mask, cd,
                                                          todo_bcast);
        SIMD_double am_perm;
        am_perm = _mm512_mask_permutexvar_pd(_mm512_undefined_pd(),
                                             now_mask, lid, v1);
        v1 = _mm512_mask_add_pd(v1, now_mask, v1, am_perm);
        am_perm = _mm512_mask_permutexvar_pd(_mm512_undefined_pd(),
                                             now_mask, lid, v2);
        v2 = _mm512_mask_add_pd(v2, now_mask, v2, am_perm);
        am_perm = _mm512_mask_permutexvar_pd(_mm512_undefined_pd(),
                                             now_mask, lid, v3);
        v3 = _mm512_mask_add_pd(v3, now_mask, v3, am_perm);
        todo_mask = _mm512_kxor(todo_mask, now_mask);
      }
    }
  }

  // ------- Bit shift operations

  inline SIMD_int operator&(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_and_epi32(one,two);
  }

  inline SIMD_int operator>>(const SIMD_int &one, const SIMD_int &two) {
    return _mm512_srlv_epi32(one,two);
  }

  inline SIMD_int operator<<(const SIMD_int &one, const unsigned two) {
    return _mm512_slli_epi32(one,two);
  }

  // -------- I/O operations

  inline void SIMD_print(const __m512i &vec) {
    for (int i = 0; i < 16; i++)
      printf("%d ",(*((int*)&(vec) + (i))));
  }

  inline void SIMD_print(const __m512 &vec) {
    for (int i = 0; i < 16; i++)
      printf("%f ",(*((float*)&(vec) + (i))));
  }

  inline void SIMD_print(const __m512d &vec) {
    for (int i = 0; i < 8; i++)
      printf("%f ",(*((double*)&(vec) + (i))));
  }

  inline void SIMD_print(const SIMD_mask &mask) {
    SIMD_print(_mm512_maskz_mov_epi32(mask,SIMD_set(1)));
  }

  inline void SIMD_print(const char *id, const SIMD_mask &mask) {
    printf("%s ",id);
    SIMD_print(mask);
    printf("\n");
  }

  inline void SIMD_print(const char *id, const SIMD_int &vec) {
    printf("%s ",id);
    SIMD_print(vec);
    printf("\n");
  }

  inline void SIMD_print(const char *id, const SIMD_float &vec) {
    printf("%s ",id);
    SIMD_print(vec);
    printf("\n");
  }

  inline void SIMD_print(const char *id, const SIMD_double &vec) {
    printf("%s ",id);
    SIMD_print(vec);
    printf("\n");
  }

  // ---------- LAMMPS operations
  #ifndef SW_GATHER_TEST
  inline void SIMD_atom_gather(const SIMD_mask &m, const float *atom,
                               const SIMD_int &i, SIMD_float &x, SIMD_float &y,
                               SIMD_float &z) {
    x = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, atom,
                                 _MM_SCALE_1);
    y = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, atom+1,
                                 _MM_SCALE_1);
    z = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, atom+2,
                                 _MM_SCALE_1);
  }

  inline void SIMD_atom_gather(const SIMD_mask &m, const float *atom,
                               const SIMD_int &i, SIMD_float &x, SIMD_float &y,
                               SIMD_float &z, SIMD_int &type) {
    x = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, atom,
                                 _MM_SCALE_1);
    y = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, atom+1,
                                 _MM_SCALE_1);
    z = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, atom+2,
                                 _MM_SCALE_1);
    type = _mm512_mask_i32gather_epi32(_mm512_undefined_epi32(), m, i, atom+3,
                                       _MM_SCALE_1);
  }
  #endif

  inline void SIMD_atom_gather(const SIMD_mask &m, const double *atom,
                               const SIMD_int &i, SIMD_double &x,
                               SIMD_double &y, SIMD_double &z) {
    x = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                 _mm512_castsi512_si256(i), atom,
                                 _MM_SCALE_2);
    y = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                 _mm512_castsi512_si256(i), atom+1,
                                 _MM_SCALE_2);
    z = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                 _mm512_castsi512_si256(i), atom+2,
                                 _MM_SCALE_2);
  }

  inline void SIMD_atom_gather(const SIMD_mask &m, const double *atom,
                               const SIMD_int &i, SIMD_double &x,
                               SIMD_double &y, SIMD_double &z, SIMD_int &type) {
    x = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                 _mm512_castsi512_si256(i), atom,
                                 _MM_SCALE_2);
    y = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                 _mm512_castsi512_si256(i), atom+1,
                                 _MM_SCALE_2);
    z = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                 _mm512_castsi512_si256(i), atom+2,
                                 _MM_SCALE_2);
    type = _mm512_mask_i32gather_epi32(_mm512_undefined_epi32(), m, i, atom+3,
                                       _MM_SCALE_2);
  }

  inline SIMD_float SIMD_ev_add(const SIMD_float &one,
                                const SIMD_float &two) {
    return _mm512_add_ps(one,two);
  }

  inline SIMD_double SIMD_ev_add(const SIMD_double &one,
                                 const SIMD_double &two) {
    return _mm512_add_pd(one,two);
  }

  inline SIMD_double SIMD_ev_add(const SIMD_double &one,
                                 const SIMD_float &two) {
    SIMD_double twod = _mm512_cvtps_pd(_mm512_castps512_ps256(two));
    SIMD_double ans = _mm512_add_pd(one,twod);
    twod = _mm512_cvtps_pd(_mm512_castps512_ps256(
                             _mm512_shuffle_f32x4(two,two,238)));
    return _mm512_add_pd(ans,twod);
  }

  inline void SIMD_jeng_update(const SIMD_mask &rmask, float *force,
                               const SIMD_int &joffset, SIMD_float &eng) {
    SIMD_float jeng;
    SIMD_conflict_pi_reduce1(rmask, joffset, eng);
    jeng = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), rmask, joffset,
                                    force, _MM_SCALE_1);
    jeng = jeng + eng;
    _mm512_mask_i32scatter_ps(force, rmask, joffset, jeng, _MM_SCALE_1);
  }

  inline void SIMD_jeng_update(const SIMD_mask &rmask, double *force,
                               const SIMD_int &joffset, SIMD_double &eng) {
    SIMD_double jeng;
    SIMD_conflict_pi_reduce1(rmask, joffset, eng);
    jeng = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask,
                                    _mm512_castsi512_si256(joffset),
                                    force, _MM_SCALE_2);
    jeng = jeng + eng;
    _mm512_mask_i32scatter_pd(force, rmask, _mm512_castsi512_si256(joffset),
                              jeng, _MM_SCALE_2);
  }

  inline void SIMD_jeng_update(const SIMD_mask &rmask, double *force,
                               const SIMD_int &joffset, SIMD_float &eng) {
    SIMD_double engd, jeng;
    engd = _mm512_cvtps_pd(_mm512_castps512_ps256(eng));
    SIMD_conflict_pi_reduce1(rmask, joffset, engd);
    jeng = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask,
                                    _mm512_castsi512_si256(joffset),
                                    force, _MM_SCALE_2);
    jeng = jeng + engd;
    _mm512_mask_i32scatter_pd(force, rmask, _mm512_castsi512_si256(joffset),
                              jeng, _MM_SCALE_2);

    SIMD_mask rmask2 = rmask >> 8;
    engd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                             _mm512_shuffle_f32x4(eng,eng,238)));
    SIMD_int joffset2 = _mm512_shuffle_i32x4(joffset, joffset, 238);
    SIMD_conflict_pi_reduce1(rmask2, joffset2, engd);
    jeng = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask2,
                                    _mm512_castsi512_si256(joffset2),
                                    force, _MM_SCALE_2);
    jeng = jeng + engd;
    _mm512_mask_i32scatter_pd(force, rmask2, _mm512_castsi512_si256(joffset2),
                              jeng, _MM_SCALE_2);
  }

  inline void SIMD_jeng_update_hi(const SIMD_mask &mask, float *force,
                                  const SIMD_int &joffset1, SIMD_float &eng) {
  }

  inline void SIMD_jeng_update_hi(const SIMD_mask &mask, double *force,
                                  const SIMD_int &joffset1, SIMD_double &eng) {
    SIMD_mask rmask = mask >> 8;
    SIMD_int joffset = _mm512_shuffle_i32x4(joffset1, joffset1, 238);

    SIMD_double jeng;
    SIMD_conflict_pi_reduce1(rmask, joffset, eng);
    jeng = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask,
                                    _mm512_castsi512_si256(joffset),
                                    force, _MM_SCALE_2);
    jeng = jeng + eng;
    _mm512_mask_i32scatter_pd(force, rmask, _mm512_castsi512_si256(joffset),
                              jeng, _MM_SCALE_2);
  }

  inline void SIMD_safe_jforce(const SIMD_mask &m, float *force,
                               const SIMD_int &i, SIMD_float &fx,
                               SIMD_float &fy, SIMD_float &fz) {
    SIMD_conflict_pi_reduce3(m, i, fx, fy, fz);
    SIMD_float jfrc;
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force,
                                    _MM_SCALE_1);
    jfrc = jfrc + fx;
    _mm512_mask_i32scatter_ps(force, m, i, jfrc, _MM_SCALE_1);
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force + 1,
                                    _MM_SCALE_1);
    jfrc = jfrc + fy;
    _mm512_mask_i32scatter_ps(force+1, m, i, jfrc, _MM_SCALE_1);
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force + 2,
                                    _MM_SCALE_1);
    jfrc = jfrc + fz;
    _mm512_mask_i32scatter_ps(force+2, m, i, jfrc, _MM_SCALE_1);
  }

  inline void SIMD_safe_jforce(const SIMD_mask &m, double *force,
                               const SIMD_int &i, SIMD_double &fx,
                               SIMD_double &fy, SIMD_double &fz) {
    SIMD_conflict_pi_reduce3(m, i, fx, fy, fz);
    SIMD_double jfrc;
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force,
                                    _MM_SCALE_2);
    jfrc = jfrc + fx;
    _mm512_mask_i32scatter_pd(force, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force + 1,
                                    _MM_SCALE_2);
    jfrc = jfrc + fy;
    _mm512_mask_i32scatter_pd(force+1, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force + 2,
                                    _MM_SCALE_2);
    jfrc = jfrc + fz;
    _mm512_mask_i32scatter_pd(force+2, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
  }

  inline void SIMD_safe_jforce(const SIMD_mask &rmask, double *force,
                               const SIMD_int &joffset, SIMD_float &amx,
                               SIMD_float &amy, SIMD_float &amz) {
    SIMD_double amxd, amyd, amzd;
    amxd = _mm512_cvtps_pd(_mm512_castps512_ps256(amx));
    amyd = _mm512_cvtps_pd(_mm512_castps512_ps256(amy));
    amzd = _mm512_cvtps_pd(_mm512_castps512_ps256(amz));
    SIMD_conflict_pi_reduce3(rmask, joffset, amxd, amyd, amzd);
    SIMD_double jfrc;
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask,
                                    _mm512_castsi512_si256(joffset),
                                    force, _MM_SCALE_2);
    jfrc = jfrc + amxd;
    _mm512_mask_i32scatter_pd(force, rmask, _mm512_castsi512_si256(joffset),
                              jfrc, _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask,
                                    _mm512_castsi512_si256(joffset),
                                    force + 1, _MM_SCALE_2);
    jfrc = jfrc + amyd;
    _mm512_mask_i32scatter_pd(force+1, rmask, _mm512_castsi512_si256(joffset),
                              jfrc, _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask,
                                    _mm512_castsi512_si256(joffset),
                                    force + 2, _MM_SCALE_2);
    jfrc = jfrc + amzd;
    _mm512_mask_i32scatter_pd(force+2, rmask, _mm512_castsi512_si256(joffset),
                              jfrc, _MM_SCALE_2);

    SIMD_mask rmask2 = rmask >> 8;
    amxd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                                  _mm512_shuffle_f32x4(amx,amx,238)));
    amyd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                                  _mm512_shuffle_f32x4(amy,amy,238)));
    amzd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                                  _mm512_shuffle_f32x4(amz,amz,238)));
    SIMD_int joffset2 = _mm512_shuffle_i32x4(joffset, joffset, 238);
    SIMD_conflict_pi_reduce3(rmask2, joffset2, amxd, amyd, amzd);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask2,
                                    _mm512_castsi512_si256(joffset2),
                                    force, _MM_SCALE_2);
    jfrc = jfrc + amxd;
    _mm512_mask_i32scatter_pd(force, rmask2, _mm512_castsi512_si256(joffset2),
                              jfrc, _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask2,
                                    _mm512_castsi512_si256(joffset2),
                                    force + 1, _MM_SCALE_2);
    jfrc = jfrc + amyd;
    _mm512_mask_i32scatter_pd(force+1, rmask2,
                              _mm512_castsi512_si256(joffset2), jfrc,
                              _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), rmask2,
                                    _mm512_castsi512_si256(joffset2),
                                    force + 2, _MM_SCALE_2);
    jfrc = jfrc + amzd;
    _mm512_mask_i32scatter_pd(force+2, rmask2,
                              _mm512_castsi512_si256(joffset2), jfrc,
                              _MM_SCALE_2);
  }

  inline void SIMD_jforce_update(const SIMD_mask &m, float *force,
                                 const SIMD_int &i, const SIMD_float &fx,
                                 const SIMD_float &fy, const SIMD_float &fz) {
    SIMD_float jfrc;
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force,
                                    _MM_SCALE_1);
    jfrc = jfrc - fx;
    _mm512_mask_i32scatter_ps(force, m, i, jfrc, _MM_SCALE_1);
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force + 1,
                                    _MM_SCALE_1);
    jfrc = jfrc - fy;
    _mm512_mask_i32scatter_ps(force+1, m, i, jfrc, _MM_SCALE_1);
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force + 2,
                                    _MM_SCALE_1);
    jfrc = jfrc - fz;
    _mm512_mask_i32scatter_ps(force+2, m, i, jfrc, _MM_SCALE_1);
  }

  template <class ft>
  inline void SIMD_scalar_update(const int jj, const int* ejnum, ft *force,
                                 const int* i, const double *fx,
                                 const double *fy, const double *fz,
                                 const double *fx2, const double *fy2,
                                 const double *fz2) {
    #pragma novector
    for (int k=0; k<8; k++) {
      if (jj < ejnum[k]) {
        const int j = i[k];
        force[j].x -= fx[k];
        force[j].y -= fy[k];
        force[j].z -= fz[k];
      }
    }

    #pragma novector
    for (int k=8; k<16; k++) {
      if (jj < ejnum[k]) {
        const int j = i[k];
        force[j].x -= fx2[k-8];
        force[j].y -= fy2[k-8];
        force[j].z -= fz2[k-8];
      }
    }
  }

  inline void SIMD_jforce_update(const SIMD_mask &m, double *force,
                                 const SIMD_int &i, const SIMD_double &fx,
                                 const SIMD_double &fy, const SIMD_double &fz)   {
    SIMD_double jfrc;
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force,
                                    _MM_SCALE_2);
    jfrc = jfrc - fx;
    _mm512_mask_i32scatter_pd(force, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force + 1,
                                    _MM_SCALE_2);
    jfrc = jfrc - fy;
    _mm512_mask_i32scatter_pd(force+1, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force + 2,
                                    _MM_SCALE_2);
    jfrc = jfrc - fz;
    _mm512_mask_i32scatter_pd(force+2, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
  }

  inline void SIMD_jforce_update(const SIMD_mask &rmask,
         double *force, const SIMD_int &joffset, SIMD_float &amx,
                                 SIMD_float &amy, SIMD_float &amz) {
    SIMD_double amxd, amyd, amzd;
    amxd = _mm512_cvtps_pd(_mm512_castps512_ps256(amx));
    amyd = _mm512_cvtps_pd(_mm512_castps512_ps256(amy));
    amzd = _mm512_cvtps_pd(_mm512_castps512_ps256(amz));
    SIMD_conflict_pi_reduce3(rmask, joffset, amxd, amyd, amzd);
    SIMD_jforce_update(rmask, force, joffset, amxd, amyd, amzd);

    SIMD_mask rmask2 = rmask >> 8;
    amxd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(amx,amx,238)));
    amyd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(amy,amy,238)));
    amzd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(amz,amz,238)));
    SIMD_int joffset2 = _mm512_shuffle_i32x4(joffset, joffset, 238);
    SIMD_conflict_pi_reduce3(rmask2, joffset2, amxd, amyd, amzd);
    SIMD_jforce_update(rmask2, force, joffset2, amxd, amyd, amzd);
  }

  inline void SIMD_cache3(float *pr, const int offset,
                          const SIMD_float &fx,
                          const SIMD_float &fy, const SIMD_float &fz) {
    float *p = pr;
    SIMD_float t;
    t = SIMD_load(p);
    t = t + fx;
    SIMD_store(p,t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fy;
    SIMD_store(p, t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fz;
    SIMD_store(p, t);
  }

  inline void SIMD_cache3(double *pr, const int offset,
                          const SIMD_double &fx,
                          const SIMD_double &fy, const SIMD_double &fz) {
    double *p = pr;
    SIMD_double t;
    t = SIMD_load(p);
    t = t + fx;
    SIMD_store(p,t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fy;
    SIMD_store(p, t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fz;
    SIMD_store(p, t);
  }

  inline void SIMD_cache3(double *pr, const int foffset,
                          const SIMD_float &fx,
                          const SIMD_float &fy, const SIMD_float &fz) {
    const int offset = foffset >> 1;
    double *p = pr;
    SIMD_double t, fd;

    fd = _mm512_cvtps_pd(_mm512_castps512_ps256(fx));
    t = SIMD_load(p);
    t = t + fd;
    SIMD_store(p,t);
    fd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fx,fx,238)));
    p = p + offset;
    t = SIMD_load(p);
    t = t + fd;
    SIMD_store(p,t);

    fd = _mm512_cvtps_pd(_mm512_castps512_ps256(fy));
    p = p + offset;
    t = SIMD_load(p);
    t = t + fd;
    SIMD_store(p,t);
    fd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fy,fy,238)));
    p = p + offset;
    t = SIMD_load(p);
    t = t + fd;
    SIMD_store(p,t);

    fd = _mm512_cvtps_pd(_mm512_castps512_ps256(fz));
    p = p + offset;
    t = SIMD_load(p);
    t = t + fd;
    SIMD_store(p,t);
    fd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fz,fz,238)));
    p = p + offset;
    t = SIMD_load(p);
    t = t + fd;
    SIMD_store(p,t);
  }

  inline void SIMD_cache3(float *pr, const int offset,
                          const SIMD_float &fx, const SIMD_float &fy,
                          const SIMD_float &fz, const SIMD_float &fx2,
                          const SIMD_float &fy2, const SIMD_float &fz2) {
  }

  inline void SIMD_cache3(double *pr, const int foffset,
                          const SIMD_double &fx, const SIMD_double &fy,
                          const SIMD_double &fz, const SIMD_double &fx2,
                          const SIMD_double &fy2, const SIMD_double &fz2) {
    const int offset = foffset >> 1;
    double *p = pr;
    SIMD_double t;

    t = SIMD_load(p);
    t = t + fx;
    SIMD_store(p,t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fx2;
    SIMD_store(p,t);

    p = p + offset;
    t = SIMD_load(p);
    t = t + fy;
    SIMD_store(p,t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fy2;
    SIMD_store(p,t);

    p = p + offset;
    t = SIMD_load(p);
    t = t + fz;
    SIMD_store(p,t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fz2;
    SIMD_store(p,t);
  }

  inline void SIMD_accumulate3(const SIMD_mask &kmask, const SIMD_float &fjx,
                               const SIMD_float &fjy, const SIMD_float &fjz,
                               SIMD_float &fxtmp, SIMD_float &fytmp,
                               SIMD_float &fztmp, SIMD_float &fjxtmp,
                               SIMD_float &fjytmp, SIMD_float &fjztmp,
                               SIMD_float &fxtmp2, SIMD_float &fytmp2,
                               SIMD_float &fztmp2, SIMD_float &fjxtmp2,
                               SIMD_float &fjytmp2, SIMD_float &fjztmp2) {
    fxtmp = SIMD_sub(fxtmp, kmask, fxtmp, fjx);
    fjxtmp = SIMD_sub(fjxtmp, kmask, fjxtmp, fjx);
    fytmp = SIMD_sub(fytmp, kmask, fytmp, fjy);
    fjytmp = SIMD_sub(fjytmp, kmask, fjytmp, fjy);
    fztmp = SIMD_sub(fztmp, kmask, fztmp, fjz);
    fjztmp = SIMD_sub(fjztmp, kmask, fjztmp, fjz);
  }

  inline void SIMD_accumulate3(const SIMD_mask &kmask, const SIMD_double &fjx,
                               const SIMD_double &fjy, const SIMD_double &fjz,
                               SIMD_double &fxtmp, SIMD_double &fytmp,
                               SIMD_double &fztmp, SIMD_double &fjxtmp,
                               SIMD_double &fjytmp, SIMD_double &fjztmp,
                               SIMD_double &fxtmp2, SIMD_double &fytmp2,
                               SIMD_double &fztmp2, SIMD_double &fjxtmp2,
                               SIMD_double &fjytmp2, SIMD_double &fjztmp2) {
    fxtmp = SIMD_sub(fxtmp, kmask, fxtmp, fjx);
    fjxtmp = SIMD_sub(fjxtmp, kmask, fjxtmp, fjx);
    fytmp = SIMD_sub(fytmp, kmask, fytmp, fjy);
    fjytmp = SIMD_sub(fjytmp, kmask, fjytmp, fjy);
    fztmp = SIMD_sub(fztmp, kmask, fztmp, fjz);
    fjztmp = SIMD_sub(fjztmp, kmask, fjztmp, fjz);
  }

  inline void SIMD_accumulate3(const SIMD_mask &kmask, const SIMD_float &fjx,
                               const SIMD_float &fjy, const SIMD_float &fjz,
                               SIMD_double &fxtmp, SIMD_double &fytmp,
                               SIMD_double &fztmp, SIMD_double &fjxtmp,
                               SIMD_double &fjytmp, SIMD_double &fjztmp,
                               SIMD_double &fxtmp2, SIMD_double &fytmp2,
                               SIMD_double &fztmp2, SIMD_double &fjxtmp2,
                               SIMD_double &fjytmp2, SIMD_double &fjztmp2) {
    SIMD_mask kmask2 = kmask >> 8;
    SIMD_double delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(fjx));
    fxtmp = SIMD_sub(fxtmp, kmask, fxtmp, delfd);
    fjxtmp = SIMD_sub(fjxtmp, kmask, fjxtmp, delfd);
    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fjx,fjx,238)));
    fxtmp2 = SIMD_sub(fxtmp2, kmask2, fxtmp2, delfd);
    fjxtmp2 = SIMD_sub(fjxtmp2, kmask2, fjxtmp2, delfd);

    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(fjy));
    fytmp = SIMD_sub(fytmp, kmask, fytmp, delfd);
    fjytmp = SIMD_sub(fjytmp, kmask, fjytmp, delfd);
    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fjy,fjy,238)));
    fytmp2 = SIMD_sub(fytmp2, kmask2, fytmp2, delfd);
    fjytmp2 = SIMD_sub(fjytmp2, kmask2, fjytmp2, delfd);

    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(fjz));
    fztmp = SIMD_sub(fztmp, kmask, fztmp, delfd);
    fjztmp = SIMD_sub(fjztmp, kmask, fjztmp, delfd);
    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fjz,fjz,238)));
    fztmp2 = SIMD_sub(fztmp2, kmask2, fztmp2, delfd);
    fjztmp2 = SIMD_sub(fjztmp2, kmask2, fjztmp2, delfd);
  }

  inline void SIMD_acc_cache3(const SIMD_mask &kmask, const SIMD_float &fjx,
                              const SIMD_float &fjy, const SIMD_float &fjz,
                              const SIMD_float &fkx, const SIMD_float &fky,
                              const SIMD_float &fkz,
                              SIMD_float &fxtmp, SIMD_float &fytmp,
                              SIMD_float &fztmp, SIMD_float &fjxtmp,
                              SIMD_float &fjytmp, SIMD_float &fjztmp,
                              SIMD_float &fxtmp2, SIMD_float &fytmp2,
                              SIMD_float &fztmp2, SIMD_float &fjxtmp2,
                              SIMD_float &fjytmp2, SIMD_float &fjztmp2,
                              float *pr, const int offset) {
    fxtmp = SIMD_sub(fxtmp, kmask, fxtmp, fjx - fkx);
    fjxtmp = SIMD_sub(fjxtmp, kmask, fjxtmp, fjx);
    fytmp = SIMD_sub(fytmp, kmask, fytmp, fjy - fky);
    fjytmp = SIMD_sub(fjytmp, kmask, fjytmp, fjy);
    fztmp = SIMD_sub(fztmp, kmask, fztmp, fjz - fkz);
    fjztmp = SIMD_sub(fjztmp, kmask, fjztmp, fjz);
    float *p = pr;
    SIMD_float t;
    t = SIMD_load(p);
    t = t + fkx;
    SIMD_store(p,t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fky;
    SIMD_store(p, t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fkz;
    SIMD_store(p, t);
  }

  inline void SIMD_acc_cache3(const SIMD_mask &kmask, const SIMD_double &fjx,
                              const SIMD_double &fjy, const SIMD_double &fjz,
                              const SIMD_double &fkx, const SIMD_double &fky,
                              const SIMD_double &fkz,
                              SIMD_double &fxtmp, SIMD_double &fytmp,
                              SIMD_double &fztmp, SIMD_double &fjxtmp,
                              SIMD_double &fjytmp, SIMD_double &fjztmp,
                              SIMD_double &fxtmp2, SIMD_double &fytmp2,
                              SIMD_double &fztmp2, SIMD_double &fjxtmp2,
                              SIMD_double &fjytmp2, SIMD_double &fjztmp2,
                              double *pr, const int offset) {
    fxtmp = SIMD_sub(fxtmp, kmask, fxtmp, fjx - fkx);
    fjxtmp = SIMD_sub(fjxtmp, kmask, fjxtmp, fjx);
    fytmp = SIMD_sub(fytmp, kmask, fytmp, fjy - fky);
    fjytmp = SIMD_sub(fjytmp, kmask, fjytmp, fjy);
    fztmp = SIMD_sub(fztmp, kmask, fztmp, fjz - fkz);
    fjztmp = SIMD_sub(fjztmp, kmask, fjztmp, fjz);
    double *p = pr;
    SIMD_double t;
    t = SIMD_load(p);
    t = t + fkx;
    SIMD_store(p,t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fky;
    SIMD_store(p, t);
    p = p + offset;
    t = SIMD_load(p);
    t = t + fkz;
    SIMD_store(p, t);
  }

  inline void SIMD_acc_cache3(const SIMD_mask &kmask, const SIMD_float &fjx,
                              const SIMD_float &fjy, const SIMD_float &fjz,
                              const SIMD_float &fkx, const SIMD_float &fky,
                              const SIMD_float &fkz,
                              SIMD_double &fxtmp, SIMD_double &fytmp,
                              SIMD_double &fztmp, SIMD_double &fjxtmp,
                              SIMD_double &fjytmp, SIMD_double &fjztmp,
                              SIMD_double &fxtmp2, SIMD_double &fytmp2,
                              SIMD_double &fztmp2, SIMD_double &fjxtmp2,
                              SIMD_double &fjytmp2, SIMD_double &fjztmp2,
                              double *pr, const int foffset) {
    SIMD_mask kmask2 = kmask >> 8;
    const int offset = foffset >> 1;
    double *p = pr;
    SIMD_double t;

    SIMD_double delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(fjx));
    SIMD_double delfdk = _mm512_cvtps_pd(_mm512_castps512_ps256(fkx));
    t = SIMD_load(p);
    t = t + delfdk;
    SIMD_store(p,t);
    fxtmp = SIMD_sub(fxtmp, kmask, fxtmp, delfd - delfdk);
    fjxtmp = SIMD_sub(fjxtmp, kmask, fjxtmp, delfd);
    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fjx,fjx,238)));
    delfdk = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fkx,fkx,238)));
    p = p + offset;
    t = SIMD_load(p);
    t = t + delfdk;
    SIMD_store(p,t);
    fxtmp2 = SIMD_sub(fxtmp2, kmask2, fxtmp2, delfd - delfdk);
    fjxtmp2 = SIMD_sub(fjxtmp2, kmask2, fjxtmp2, delfd);

    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(fjy));
    delfdk = _mm512_cvtps_pd(_mm512_castps512_ps256(fky));
    p = p + offset;
    t = SIMD_load(p);
    t = t + delfdk;
    SIMD_store(p,t);
    fytmp = SIMD_sub(fytmp, kmask, fytmp, delfd - delfdk);
    fjytmp = SIMD_sub(fjytmp, kmask, fjytmp, delfd);
    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fjy,fjy,238)));
    delfdk = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fky,fky,238)));
    p = p + offset;
    t = SIMD_load(p);
    t = t + delfdk;
    SIMD_store(p,t);
    fytmp2 = SIMD_sub(fytmp2, kmask2, fytmp2, delfd - delfdk);
    fjytmp2 = SIMD_sub(fjytmp2, kmask2, fjytmp2, delfd);

    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(fjz));
    delfdk = _mm512_cvtps_pd(_mm512_castps512_ps256(fkz));
    p = p + offset;
    t = SIMD_load(p);
    t = t + delfdk;
    SIMD_store(p,t);
    fztmp = SIMD_sub(fztmp, kmask, fztmp, delfd - delfdk);
    fjztmp = SIMD_sub(fjztmp, kmask, fjztmp, delfd);
    delfd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fjz,fjz,238)));
    delfdk = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(fkz,fkz,238)));
    p = p + offset;
    t = SIMD_load(p);
    t = t + delfdk;
    SIMD_store(p,t);
    fztmp2 = SIMD_sub(fztmp2, kmask2, fztmp2, delfd - delfdk);
    fjztmp2 = SIMD_sub(fjztmp2, kmask2, fjztmp2, delfd);
  }

  inline void SIMD_acc_energy3(const SIMD_mask &hmask,
                               const SIMD_float &evdwl, const int eatom,
                               SIMD_float &sevdwl, SIMD_float &fwtmp,
                               SIMD_float &fjtmp, SIMD_float &fwtmp2,
                               SIMD_float &fjtmp2) {
    sevdwl = SIMD_add(sevdwl, hmask, sevdwl, evdwl);
    if (eatom) {
      const SIMD_float hevdwl = evdwl * (float)0.5;
      fwtmp = SIMD_add(fwtmp, hmask, fwtmp, hevdwl);
      fjtmp = SIMD_add(fjtmp, hmask, fjtmp, hevdwl);
    }
  }

  inline void SIMD_acc_energy3(const SIMD_mask &hmask,
                               const SIMD_double &evdwl, const int eatom,
                               SIMD_double &sevdwl, SIMD_double &fwtmp,
                               SIMD_double &fjtmp, SIMD_double &fwtmp2,
                               SIMD_double &fjtmp2) {
    sevdwl = SIMD_add(sevdwl, hmask, sevdwl, evdwl);
    if (eatom) {
      const SIMD_double hevdwl = evdwl * (double)0.5;
      fwtmp = SIMD_add(fwtmp, hmask, fwtmp, hevdwl);
      fjtmp = SIMD_add(fjtmp, hmask, fjtmp, hevdwl);
    }
  }

  inline void SIMD_acc_energy3(const SIMD_mask &hmask,
                               const SIMD_float &evdwl, const int eatom,
                               SIMD_double &sevdwl, SIMD_double &fwtmp,
                               SIMD_double &fjtmp, SIMD_double &fwtmp2,
                               SIMD_double &fjtmp2) {
    SIMD_double evdwld;
    evdwld = _mm512_cvtps_pd(_mm512_castps512_ps256(evdwl));
    sevdwl = SIMD_add(sevdwl, hmask, sevdwl, evdwld);
    if (eatom) {
      const SIMD_double hevdwl = evdwld * (double)0.5;
      fwtmp = SIMD_add(fwtmp, hmask, fwtmp, hevdwl);
      fjtmp = SIMD_add(fjtmp, hmask, fjtmp, hevdwl);
    }
    SIMD_mask hmask2 = hmask >> 8;
    evdwld = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(evdwl,evdwl,238)));
    sevdwl = SIMD_add(sevdwl, hmask2, sevdwl, evdwld);
    if (eatom) {
      const SIMD_double hevdwl = evdwld * (double)0.5;
      fwtmp2 = SIMD_add(fwtmp2, hmask2, fwtmp2, hevdwl);
      fjtmp2 = SIMD_add(fjtmp2, hmask2, fjtmp2, hevdwl);
    }
  }

  inline void SIMD_acc_three(const SIMD_mask &hmask, const SIMD_float &facrad,
                             const int eatom, SIMD_float &sevdwl,
                             SIMD_float &fwtmp, SIMD_float &fjtmp,
                             SIMD_float &fwtmp2, SIMD_float &fjtmp2,
                             const SIMD_int &k, float *force) {
    sevdwl = SIMD_add(sevdwl, hmask, sevdwl, facrad);
    if (eatom) {
      SIMD_float hevdwl = facrad * SIMD_set((float)0.33333333);
      fwtmp = SIMD_add(fwtmp, hmask, fwtmp, hevdwl);
      fjtmp = SIMD_add(fjtmp, hmask, fjtmp, hevdwl);
      SIMD_conflict_pi_reduce1(hmask, k, hevdwl);
      SIMD_float keng = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), hmask,
                                                 k, force + 3, _MM_SCALE_1);
      keng = keng + hevdwl;
      _mm512_mask_i32scatter_ps(force + 3, hmask, k, keng, _MM_SCALE_1);
    }
  }

  inline void SIMD_acc_three(const SIMD_mask &hmask, const SIMD_double &facrad,
                             const int eatom, SIMD_double &sevdwl,
                             SIMD_double &fwtmp, SIMD_double &fjtmp,
                             SIMD_double &fwtmp2, SIMD_double &fjtmp2,
                             const SIMD_int &k, double *force) {
    sevdwl = SIMD_add(sevdwl, hmask, sevdwl, facrad);
    if (eatom) {
      SIMD_double hevdwl = facrad * SIMD_set((double)0.33333333);
      fwtmp = SIMD_add(fwtmp, hmask, fwtmp, hevdwl);
      fjtmp = SIMD_add(fjtmp, hmask, fjtmp, hevdwl);
      SIMD_conflict_pi_reduce1(hmask, k, hevdwl);
      SIMD_double keng = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), hmask,
                                                  _mm512_castsi512_si256(k),
                                                  force + 3, _MM_SCALE_2);
      keng = keng + hevdwl;
      _mm512_mask_i32scatter_pd(force + 3, hmask, _mm512_castsi512_si256(k),
                                keng, _MM_SCALE_2);
    }
  }

  inline void SIMD_acc_three(const SIMD_mask &hmask, const SIMD_float &facrad,
                             const int eatom, SIMD_double &sevdwl,
                             SIMD_double &fwtmp, SIMD_double &fjtmp,
                             SIMD_double &fwtmp2, SIMD_double &fjtmp2,
                             const SIMD_int &k, double *force) {
    SIMD_double facradd;
    facradd = _mm512_cvtps_pd(_mm512_castps512_ps256(facrad));
    sevdwl = SIMD_add(sevdwl, hmask, sevdwl, facradd);
    if (eatom) {
      SIMD_double hevdwl = facradd * SIMD_set((double)0.33333333);
      fwtmp = SIMD_add(fwtmp, hmask, fwtmp, hevdwl);
      fjtmp = SIMD_add(fjtmp, hmask, fjtmp, hevdwl);
      SIMD_conflict_pi_reduce1(hmask, k, hevdwl);
      SIMD_double keng = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), hmask,
                                                  _mm512_castsi512_si256(k),
                                                  force + 3, _MM_SCALE_2);
      keng = keng + hevdwl;
      _mm512_mask_i32scatter_pd(force + 3, hmask, _mm512_castsi512_si256(k),
                                keng, _MM_SCALE_2);
    }
    SIMD_mask hmask2 = hmask >> 8;
    facradd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(facrad,facrad,238)));
    sevdwl = SIMD_add(sevdwl, hmask2, sevdwl, facradd);
    if (eatom) {
      SIMD_double hevdwl = facradd * SIMD_set((double)0.33333333);
      fwtmp2 = SIMD_add(fwtmp2, hmask2, fwtmp2, hevdwl);
      fjtmp2 = SIMD_add(fjtmp2, hmask2, fjtmp2, hevdwl);
      SIMD_int k2 = _mm512_shuffle_i32x4(k, k, 238);
      SIMD_conflict_pi_reduce1(hmask2, k2, hevdwl);
      SIMD_double keng = _mm512_mask_i32gather_pd(_mm512_undefined_pd(),
                                                  hmask2,
                                                  _mm512_castsi512_si256(k2),
                                                  force + 3, _MM_SCALE_2);
      keng = keng + hevdwl;
      _mm512_mask_i32scatter_pd(force + 3, hmask2, _mm512_castsi512_si256(k2),
                                keng, _MM_SCALE_2);
    }
  }

  inline void SIMD_ev_tally_nbor(const SIMD_mask &m, const int vflag,
                                 const float ev_pre,
                 const SIMD_float &fpair, const SIMD_float &delx,
                 const SIMD_float &dely,  const SIMD_float &delz,
                 SIMD_float &sv0, SIMD_float &sv1, SIMD_float &sv2,
                 SIMD_float &sv3, SIMD_float &sv4, SIMD_float &sv5) {
    if (vflag == 1) {
      const SIMD_float prefpair = SIMD_set(ev_pre) * fpair;
      sv0 = SIMD_add(sv0, m, sv0, delx * delx * prefpair);
      sv1 = SIMD_add(sv1, m, sv1, dely * dely * prefpair);
      sv2 = SIMD_add(sv2, m, sv2, delz * delz * prefpair);
      sv3 = SIMD_add(sv3, m, sv3, delx * dely * prefpair);
      sv4 = SIMD_add(sv4, m, sv4, delx * delz * prefpair);
      sv5 = SIMD_add(sv5, m, sv5, dely * delz * prefpair);
    }
  }

  inline void SIMD_ev_tally_nbor(const SIMD_mask &m, const int vflag,
                                 const double ev_pre,
                 const SIMD_double &fpair, const SIMD_double &delx,
                 const SIMD_double &dely,  const SIMD_double &delz,
                 SIMD_double &sv0, SIMD_double &sv1, SIMD_double &sv2,
                 SIMD_double &sv3, SIMD_double &sv4, SIMD_double &sv5) {
    if (vflag == 1) {
      const SIMD_double prefpair = SIMD_set(ev_pre) * fpair;
      sv0 = SIMD_add(sv0, m, sv0, delx * delx * prefpair);
      sv1 = SIMD_add(sv1, m, sv1, dely * dely * prefpair);
      sv2 = SIMD_add(sv2, m, sv2, delz * delz * prefpair);
      sv3 = SIMD_add(sv3, m, sv3, delx * dely * prefpair);
      sv4 = SIMD_add(sv4, m, sv4, delx * delz * prefpair);
      sv5 = SIMD_add(sv5, m, sv5, dely * delz * prefpair);
    }
  }

  inline void SIMD_ev_tally_nbor(const SIMD_mask &m, const int vflag,
                                 const float ev_pre,
                 const SIMD_float &fpair, const SIMD_float &delx,
                 const SIMD_float &dely,  const SIMD_float &delz,
                 SIMD_double &sv0, SIMD_double &sv1, SIMD_double &sv2,
                 SIMD_double &sv3, SIMD_double &sv4, SIMD_double &sv5) {
    if (vflag == 1) {
      const SIMD_mask m2 = m >> 8;
      const SIMD_float prefpair = SIMD_set(ev_pre) * fpair;
      SIMD_float dpair = delx * delx * prefpair;
      SIMD_double dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv0 = SIMD_add(sv0, m, sv0, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv0 = SIMD_add(sv0, m2, sv0, dpaird);

      dpair = dely * dely * prefpair;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv1 = SIMD_add(sv1, m, sv1, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv1 = SIMD_add(sv1, m2, sv1, dpaird);

      dpair = delz * delz * prefpair;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv2 = SIMD_add(sv2, m, sv2, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv2 = SIMD_add(sv2, m2, sv2, dpaird);

      dpair = delx * dely * prefpair;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv3 = SIMD_add(sv3, m, sv3, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv3 = SIMD_add(sv3, m2, sv3, dpaird);

      dpair = delx * delz * prefpair;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv4 = SIMD_add(sv4, m, sv4, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv4 = SIMD_add(sv4, m2, sv4, dpaird);

      dpair = dely * delz * prefpair;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv5 = SIMD_add(sv5, m, sv5, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv5 = SIMD_add(sv5, m2, sv5, dpaird);
    }
  }

  inline void SIMD_ev_tally_nbor3v(const SIMD_mask &m, const int vflag,
                 const SIMD_float &fj0, const SIMD_float &fj1,
                 const SIMD_float &fj2, const SIMD_float &fk0,
                 const SIMD_float &fk1, const SIMD_float &fk2,
                 const SIMD_float &delx, const SIMD_float &dely,
                 const SIMD_float &delz, const SIMD_float &delr2x,
                 const SIMD_float &delr2y, const SIMD_float &delr2z,
                 SIMD_float &sv0, SIMD_float &sv1, SIMD_float &sv2,
                 SIMD_float &sv3, SIMD_float &sv4, SIMD_float &sv5) {
    if (vflag == 1) {
      sv0 = SIMD_add(sv0, m, sv0, delx * fj0 + delr2x * fk0);
      sv1 = SIMD_add(sv1, m, sv1, dely * fj1 + delr2y * fk1);
      sv2 = SIMD_add(sv2, m, sv2, delz * fj2 + delr2z * fk2);
      sv3 = SIMD_add(sv3, m, sv3, delx * fj1 + delr2x * fk1);
      sv4 = SIMD_add(sv4, m, sv4, delx * fj2 + delr2x * fk2);
      sv5 = SIMD_add(sv5, m, sv5, dely * fj2 + delr2y * fk2);
    }
  }

  inline void SIMD_ev_tally_nbor3v(const SIMD_mask &m, const int vflag,
                 const SIMD_double &fj0, const SIMD_double &fj1,
                 const SIMD_double &fj2, const SIMD_double &fk0,
                 const SIMD_double &fk1, const SIMD_double &fk2,
                 const SIMD_double &delx, const SIMD_double &dely,
                 const SIMD_double &delz, const SIMD_double &delr2x,
                 const SIMD_double &delr2y, const SIMD_double &delr2z,
                 SIMD_double &sv0, SIMD_double &sv1, SIMD_double &sv2,
                 SIMD_double &sv3, SIMD_double &sv4, SIMD_double &sv5) {
    if (vflag == 1) {
      sv0 = SIMD_add(sv0, m, sv0, delx * fj0 + delr2x * fk0);
      sv1 = SIMD_add(sv1, m, sv1, dely * fj1 + delr2y * fk1);
      sv2 = SIMD_add(sv2, m, sv2, delz * fj2 + delr2z * fk2);
      sv3 = SIMD_add(sv3, m, sv3, delx * fj1 + delr2x * fk1);
      sv4 = SIMD_add(sv4, m, sv4, delx * fj2 + delr2x * fk2);
      sv5 = SIMD_add(sv5, m, sv5, dely * fj2 + delr2y * fk2);
    }
  }

  inline void SIMD_ev_tally_nbor3v(const SIMD_mask &m, const int vflag,
                 const SIMD_float &fj0, const SIMD_float &fj1,
                 const SIMD_float &fj2, const SIMD_float &fk0,
                 const SIMD_float &fk1, const SIMD_float &fk2,
                 const SIMD_float &delx, const SIMD_float &dely,
                 const SIMD_float &delz, const SIMD_float &delr2x,
                 const SIMD_float &delr2y, const SIMD_float &delr2z,
                 SIMD_double &sv0, SIMD_double &sv1, SIMD_double &sv2,
                 SIMD_double &sv3, SIMD_double &sv4, SIMD_double &sv5) {
    if (vflag == 1) {
      const SIMD_mask m2 = m >> 8;
      SIMD_float dpair = delx * fj0 + delr2x * fk0;
      SIMD_double dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv0 = SIMD_add(sv0, m, sv0, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv0 = SIMD_add(sv0, m2, sv0, dpaird);

      dpair = dely * fj1 + delr2y * fk1;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv1 = SIMD_add(sv1, m, sv1, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                              _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv1 = SIMD_add(sv1, m2, sv1, dpaird);

      dpair = delz * fj2 + delr2z * fk2;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv2 = SIMD_add(sv2, m, sv2, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                              _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv2 = SIMD_add(sv2, m2, sv2, dpaird);

      dpair = delx * fj1 + delr2x * fk1;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv3 = SIMD_add(sv3, m, sv3, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                              _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv3 = SIMD_add(sv3, m2, sv3, dpaird);

      dpair = delx * fj2 + delr2x * fk2;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv4 = SIMD_add(sv4, m, sv4, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                              _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv4 = SIMD_add(sv4, m2, sv4, dpaird);

      dpair = dely * fj2 + delr2y * fk2;
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(dpair));
      sv5 = SIMD_add(sv5, m, sv5, dpaird);
      dpaird = _mm512_cvtps_pd(_mm512_castps512_ps256(
                              _mm512_shuffle_f32x4(dpair,dpair,238)));
      sv5 = SIMD_add(sv5, m2, sv5, dpaird);
    }
  }

  inline void SIMD_safe_force_accumulate(const SIMD_mask &rmask,
         float *force, const SIMD_int &joffset, SIMD_float &amx,
         SIMD_float &amy, SIMD_float &amz, SIMD_float &fxtmp,
         SIMD_float &fytmp, SIMD_float &fztmp, SIMD_float &fxtmp2,
         SIMD_float &fytmp2, SIMD_float &fztmp2) {
    fxtmp = SIMD_add(fxtmp, rmask, fxtmp, amx);
    fytmp = SIMD_add(fytmp, rmask, fytmp, amy);
    fztmp = SIMD_add(fztmp, rmask, fztmp, amz);
    SIMD_conflict_pi_reduce3(rmask, joffset, amx, amy, amz);
    SIMD_jforce_update(rmask, force, joffset, amx, amy, amz);
  }

  inline void SIMD_safe_force_accumulate(const SIMD_mask &rmask,
         double *force, const SIMD_int &joffset, SIMD_double &amx,
         SIMD_double &amy, SIMD_double &amz, SIMD_double &fxtmp,
         SIMD_double &fytmp, SIMD_double &fztmp, SIMD_double &fxtmp2,
         SIMD_double &fytmp2, SIMD_double &fztmp2) {
    fxtmp = SIMD_add(fxtmp, rmask, fxtmp, amx);
    fytmp = SIMD_add(fytmp, rmask, fytmp, amy);
    fztmp = SIMD_add(fztmp, rmask, fztmp, amz);
    SIMD_conflict_pi_reduce3(rmask, joffset, amx, amy, amz);
    SIMD_jforce_update(rmask, force, joffset, amx, amy, amz);
  }

  inline void SIMD_safe_force_accumulate(const SIMD_mask &rmask,
         double *force, const SIMD_int &joffset, SIMD_float &amx,
         SIMD_float &amy, SIMD_float &amz, SIMD_double &fxtmp,
         SIMD_double &fytmp, SIMD_double &fztmp, SIMD_double &fxtmp2,
         SIMD_double &fytmp2, SIMD_double &fztmp2) {
    SIMD_double amxd, amyd, amzd;
    amxd = _mm512_cvtps_pd(_mm512_castps512_ps256(amx));
    fxtmp = SIMD_add(fxtmp, rmask, fxtmp, amxd);
    amyd = _mm512_cvtps_pd(_mm512_castps512_ps256(amy));
    fytmp = SIMD_add(fytmp, rmask, fytmp, amyd);
    amzd = _mm512_cvtps_pd(_mm512_castps512_ps256(amz));
    fztmp = SIMD_add(fztmp, rmask, fztmp, amzd);
    SIMD_conflict_pi_reduce3(rmask, joffset, amxd, amyd, amzd);
    SIMD_jforce_update(rmask, force, joffset, amxd, amyd, amzd);

    SIMD_mask rmask2 = rmask >> 8;
    amxd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(amx,amx,238)));
    fxtmp2 = SIMD_add(fxtmp2, rmask2, fxtmp2, amxd);
    amyd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(amy,amy,238)));
    fytmp2 = SIMD_add(fytmp2, rmask2, fytmp2, amyd);
    amzd = _mm512_cvtps_pd(_mm512_castps512_ps256(
                                _mm512_shuffle_f32x4(amz,amz,238)));
    fztmp2 = SIMD_add(fztmp2, rmask2, fztmp2, amzd);
    SIMD_int joffset2 = _mm512_shuffle_i32x4(joffset, joffset, 238);
    SIMD_conflict_pi_reduce3(rmask2, joffset2, amxd, amyd, amzd);
    SIMD_jforce_update(rmask2, force, joffset2, amxd, amyd, amzd);
  }

  inline void SIMD_iforce_update(const SIMD_mask &m, float *force,
                                 const SIMD_int &i, const SIMD_float &fx,
                                 const SIMD_float &fy, const SIMD_float &fz,
                                 const int EFLAG, const int eatom,
                                 const SIMD_float &fwtmp) {
    SIMD_float jfrc;
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force,
                                    _MM_SCALE_1);
    jfrc = jfrc + fx;
    _mm512_mask_i32scatter_ps(force, m, i, jfrc, _MM_SCALE_1);
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force + 1,
                                    _MM_SCALE_1);
    jfrc = jfrc + fy;
    _mm512_mask_i32scatter_ps(force+1, m, i, jfrc, _MM_SCALE_1);
    jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force + 2,
                                    _MM_SCALE_1);
    jfrc = jfrc + fz;
    _mm512_mask_i32scatter_ps(force+2, m, i, jfrc, _MM_SCALE_1);
    if (EFLAG) {
      if (eatom) {
        jfrc = _mm512_mask_i32gather_ps(_mm512_undefined_ps(), m, i, force + 3,
                                        _MM_SCALE_1);
        jfrc = jfrc + fwtmp;
        _mm512_mask_i32scatter_ps(force+3, m, i, jfrc, _MM_SCALE_1);
      }
    }
  }

  inline void SIMD_iforce_update(const SIMD_mask &m, double *force,
                                 const SIMD_int &i, const SIMD_double &fx,
                                 const SIMD_double &fy, const SIMD_double &fz,
                                 const int EFLAG, const int eatom,
                                 const SIMD_double &fwtmp) {
    SIMD_double jfrc;
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force,
                                    _MM_SCALE_2);
    jfrc = jfrc + fx;
    _mm512_mask_i32scatter_pd(force, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force + 1,
                                    _MM_SCALE_2);
    jfrc = jfrc + fy;
    _mm512_mask_i32scatter_pd(force+1, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
    jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                    _mm512_castsi512_si256(i), force + 2,
                                    _MM_SCALE_2);
    jfrc = jfrc + fz;
    _mm512_mask_i32scatter_pd(force+2, m, _mm512_castsi512_si256(i), jfrc,
                              _MM_SCALE_2);
    if (EFLAG) {
      if (eatom) {
        jfrc = _mm512_mask_i32gather_pd(_mm512_undefined_pd(), m,
                                        _mm512_castsi512_si256(i),
                                        force + 3, _MM_SCALE_2);
        jfrc = jfrc + fwtmp;
        _mm512_mask_i32scatter_pd(force+3, m, _mm512_castsi512_si256(i), jfrc,
                                  _MM_SCALE_2);
      }
    }
  }

  #ifdef SW_GATHER_TEST
  template <class atom_t>
  inline void SIMD_atom_gather(const SIMD_mask &m, const atom_t *atom,
                               const SIMD_int &i, SIMD_float &x, SIMD_float &y,
                               SIMD_float &z, SIMD_int &type) {
    int jv_scalar[16] __attribute__((aligned(64)));
    int jm_scalar[16] __attribute__((aligned(64)));
    _mm512_store_epi32(jv_scalar, i);

    SIMD_float pl1, pl2, pl3, pl4;

    int js = jv_scalar[0];
    pl1 = _mm512_loadu_ps((float *)((char *)atom + js));
    js = jv_scalar[1];
    pl1 = _mm512_insertf32x4(pl1, _mm_load_ps((float *)((char *)atom +
                                                        js)), 1);
    js = jv_scalar[2];
    pl1 = _mm512_insertf32x4(pl1, _mm_load_ps((float *)((char *)atom +
                                                        js)), 2);
    js = jv_scalar[3];
    pl1 = _mm512_insertf32x4(pl1, _mm_load_ps((float *)((char *)atom +
                                                        js)), 3);

    js = jv_scalar[4];
    pl2 = _mm512_loadu_ps((float *)((char *)atom + js));
    js = jv_scalar[5];
    pl2 = _mm512_insertf32x4(pl2, _mm_load_ps((float *)((char *)atom +
                                                        js)), 1);
    js = jv_scalar[6];
    pl2 = _mm512_insertf32x4(pl2, _mm_load_ps((float *)((char *)atom +
                                                        js)), 2);
    js = jv_scalar[7];
    pl2 = _mm512_insertf32x4(pl2, _mm_load_ps((float *)((char *)atom +
                                                        js)), 3);

    js = jv_scalar[8];
    pl3 = _mm512_loadu_ps((float *)((char *)atom + js));
    js = jv_scalar[9];
    pl3 = _mm512_insertf32x4(pl3, _mm_load_ps((float *)((char *)atom +
                                                        js)), 1);
    js = jv_scalar[10];
    pl3 = _mm512_insertf32x4(pl3, _mm_load_ps((float *)((char *)atom +
                                                        js)), 2);
    js = jv_scalar[11];
    pl3 = _mm512_insertf32x4(pl3, _mm_load_ps((float *)((char *)atom +
                                                        js)), 3);

    js = jv_scalar[12];
    pl4 = _mm512_loadu_ps((float *)((char *)atom + js));
    js = jv_scalar[13];
    pl4 = _mm512_insertf32x4(pl4, _mm_load_ps((float *)((char *)atom +
                                                        js)), 1);
    js = jv_scalar[14];
    pl4 = _mm512_insertf32x4(pl4, _mm_load_ps((float *)((char *)atom +
                                                        js)), 2);
    js = jv_scalar[15];
    pl4 = _mm512_insertf32x4(pl4, _mm_load_ps((float *)((char *)atom +
                                                        js)), 3);

    SIMD_int c0 = _mm512_setr_epi32(0x0,0x4,0x8,0xc,0x10,0x14,0x18,0x1c,
                                    0x1,0x5,0x9,0xd,0x11,0x15,0x19,0x1d);
    SIMD_int c1 = _mm512_setr_epi32(0x1,0x5,0x9,0xd,0x11,0x15,0x19,0x1d,
                                    0x0,0x4,0x8,0xc,0x10,0x14,0x18,0x1c);
    SIMD_int c2 = _mm512_setr_epi32(0x2,0x6,0xa,0xe,0x12,0x16,0x1a,0x1e,
                                    0x3,0x7,0xb,0xf,0x13,0x17,0x1b,0x1f);
    SIMD_int c3 = _mm512_setr_epi32(0x3,0x7,0xb,0xf,0x13,0x17,0x1b,0x1f,
                                    0x2,0x6,0xa,0xe,0x12,0x16,0x1a,0x1e);
    SIMD_mask k_1 = _mm512_int2mask(65280);

    SIMD_float sl1 = _mm512_permutex2var_ps(pl3, c0, pl4);
    SIMD_float sl2 = _mm512_permutex2var_ps(pl1, c1, pl2);
    SIMD_float sl3 = _mm512_permutex2var_ps(pl3, c2, pl4);
    SIMD_float sl4 = _mm512_permutex2var_ps(pl1, c3, pl2);

    x = _mm512_shuffle_f32x4(sl2, sl1, 78);
    z = _mm512_shuffle_f32x4(sl4, sl3, 78);
    y = _mm512_mask_blend_ps(k_1, sl2, sl1);
    type = _mm512_castps_si512(_mm512_mask_blend_ps(k_1, sl4, sl3));
  }
  #endif
}

#endif

#endif
