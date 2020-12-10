#ifndef LMP_INTEL_AIREBO_SCALAR
# ifdef __INTEL_COMPILER
#  if defined(__MIC__) || defined(__AVX512F__)
#   define LMP_INTEL_AIREBO_512
#  elif defined(__AVX__)
#   define LMP_INTEL_AIREBO_256
#  else
#   define LMP_INTEL_AIREBO_SCALAR
#  endif
# else
#  define LMP_INTEL_AIREBO_SCALAR
# endif
#endif

#ifdef LMP_INTEL_AIREBO_512

#include <cassert>
#include <immintrin.h>

#define VEC_INLINE __attribute__((always_inline))


#ifndef FVEC_FIRST_PASS
#  define FVEC_LEN 8
#  define FVEC_SUFFIX(a) a##pd
#  define FVEC_SUFFIX_MASK(a) a##pd_mask
#  define FVEC_MASK_T __mmask8
#  define FVEC_VEC_T __m512d
#  define FVEC_SCAL_T double
#  define IVEC_NAME ivec8
#  define FVEC_NAME fvec8pd
#  define BVEC_NAME bvec8
#  define AVEC_NAME avec8pd
#else
#  undef FVEC_LEN
#  undef FVEC_SUFFIX
#  undef FVEC_SUFFIX_MASK
#  undef FVEC_MASK_T
#  undef FVEC_VEC_T
#  undef FVEC_SCAL_T
#  undef IVEC_NAME
#  undef FVEC_NAME
#  undef BVEC_NAME
#  undef AVEC_NAME

#  define FVEC_LEN 16
#  define FVEC_SUFFIX(a) a##ps
#  define FVEC_SUFFIX_MASK(a) a##ps_mask
#  define FVEC_MASK_T __mmask16
#  define FVEC_VEC_T __m512
#  define FVEC_SCAL_T float
#  define IVEC_NAME ivec16
#  define FVEC_NAME fvec16ps
#  define BVEC_NAME bvec16
#  define AVEC_NAME avec16ps
#endif

namespace mm512 {

#ifndef __AVX512F__

#ifndef FVEC_FIRST_PASS
VEC_INLINE static inline __m512i _mm512_mask_expand_epi32(__m512i src,
                                                          __mmask16 k,
                                                          __m512i a) {
  int buf[16] __attribute__((aligned(64)));
  _mm512_store_epi32(buf, a);
  return _mm512_mask_loadunpacklo_epi32(src, k, buf);
}
VEC_INLINE static inline __m512i _mm512_maskz_expand_epi32(__mmask16 k,
                                                           __m512i a) {
  int buf[16] __attribute__((aligned(64)));
  _mm512_store_epi32(buf, a);
  return _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), k, buf);
}
VEC_INLINE static inline __m512i _mm512_mask_compress_epi32(__m512i src,
                                                            __mmask16 k,
                                                            __m512i a) {
  int buf[16] __attribute__((aligned(64)));
  _mm512_store_epi32(buf, src);
  _mm512_mask_packstorelo_epi32(buf, k, a);
  return _mm512_load_epi32(buf);
}
VEC_INLINE static inline __m512i _mm512_maskz_compress_epi32(__mmask16 k,
                                                             __m512i a) {
  int buf[16] __attribute__((aligned(64))) = {0};
  _mm512_mask_packstorelo_epi32(buf, k, a);
  return _mm512_load_epi32(buf);
}

VEC_INLINE static inline void _mm512_mask_compressstoreu_epi32(int * dest,
                                                               __mmask16 mask,
                                                               __m512i src) {
  _mm512_mask_packstorelo_epi32(dest, mask, src);
  _mm512_mask_packstorehi_epi32(dest + 16, mask, src);
}

VEC_INLINE static inline __m512i _mm512_mask_loadu_epi32(__m512i src,
                                                         __mmask16 k,
                                                         const int * mem_addr) {
  assert((k & (k + 1)) == 0);
  __m512i ret = _mm512_mask_loadunpacklo_epi32(src, k, mem_addr);
  ret = _mm512_mask_loadunpackhi_epi32(ret, k, mem_addr + 16);
  return ret;
}
VEC_INLINE static inline __m512i _mm512_maskz_loadu_epi32(__mmask16 k,
                                                        const int * mem_addr) {
  assert((k & (k + 1)) == 0);
  __m512i ret = _mm512_mask_loadunpacklo_epi32(_mm512_setzero_epi32(), k,
                                               mem_addr);
  ret = _mm512_mask_loadunpackhi_epi32(ret, k, mem_addr + 16);
  return ret;
}
VEC_INLINE static inline void _mm512_mask_storeu_epi32(int * dest,
                                                       __mmask16 mask,
                                                       __m512i src) {
  assert((mask & (mask + 1)) == 0);
  _mm512_mask_packstorelo_epi32(dest, mask, src);
  _mm512_mask_packstorehi_epi32(dest + 16, mask, src);
}
#endif

VEC_INLINE static inline FVEC_VEC_T FVEC_SUFFIX(_mm512_mask_expand_)
  (FVEC_VEC_T src, __mmask16 k, FVEC_VEC_T a) {
  FVEC_SCAL_T buf[FVEC_LEN] __attribute__((aligned(64)));
  FVEC_SUFFIX(_mm512_store_)(buf, a);
  return FVEC_SUFFIX(_mm512_mask_loadunpacklo_)(src, k, buf);
}
VEC_INLINE static inline FVEC_VEC_T FVEC_SUFFIX(_mm512_maskz_expand_)
  (__mmask16 k, FVEC_VEC_T a) {
  FVEC_SCAL_T buf[FVEC_LEN] __attribute__((aligned(64)));
  FVEC_SUFFIX(_mm512_store_)(buf, a);
  return FVEC_SUFFIX(_mm512_mask_loadunpacklo_)(FVEC_SUFFIX(_mm512_setzero_)(),
                                                k, buf);
}
VEC_INLINE static inline FVEC_VEC_T FVEC_SUFFIX(_mm512_mask_compress_)
  (FVEC_VEC_T src, __mmask16 k, FVEC_VEC_T a) {
  FVEC_SCAL_T buf[FVEC_LEN] __attribute__((aligned(64)));
  FVEC_SUFFIX(_mm512_store_)(buf, src);
  FVEC_SUFFIX(_mm512_mask_packstorelo_)(buf, k, a);
  return FVEC_SUFFIX(_mm512_load_)(buf);
}
VEC_INLINE static inline FVEC_VEC_T FVEC_SUFFIX(_mm512_maskz_compress_)
  (__mmask16 k, FVEC_VEC_T a) {
  FVEC_SCAL_T buf[FVEC_LEN] __attribute__((aligned(64))) = {0};
  FVEC_SUFFIX(_mm512_mask_packstorelo_)(buf, k, a);
  return FVEC_SUFFIX(_mm512_load_)(buf);
}
VEC_INLINE static inline void FVEC_SUFFIX(_mm512_mask_storeu_)
  (FVEC_SCAL_T * dest, FVEC_MASK_T mask, FVEC_VEC_T src) {
  assert((mask & (mask + 1)) == 0);
  FVEC_SUFFIX(_mm512_mask_packstorelo_)(dest, mask, src);
  FVEC_SUFFIX(_mm512_mask_packstorehi_)(dest + FVEC_LEN, mask, src);
}
#endif


class FVEC_NAME;
class IVEC_NAME;
class AVEC_NAME;
class BVEC_NAME {
  friend class FVEC_NAME;
  friend class IVEC_NAME;
  friend class AVEC_NAME;
# if FVEC_LEN==16
  friend class avec16pd;
# endif
  FVEC_MASK_T val_;
  VEC_INLINE BVEC_NAME(const FVEC_MASK_T &v) : val_(v) {}
public:
  VEC_INLINE BVEC_NAME() {}
  VEC_INLINE static BVEC_NAME kand(const BVEC_NAME &a, const BVEC_NAME &b) {
    return _mm512_kand(a.val_, b.val_);
  }
  VEC_INLINE static BVEC_NAME kandn(const BVEC_NAME &a, const BVEC_NAME &b) {
    return _mm512_kandn(a.val_, b.val_);
  }
  VEC_INLINE static BVEC_NAME knot(const BVEC_NAME &a) {
    return _mm512_knot(a.val_);
  }
  VEC_INLINE static int kortestz(const BVEC_NAME &a, const BVEC_NAME &b) {
    return _mm512_kortestz(a.val_, b.val_);
  }
  VEC_INLINE static BVEC_NAME masku_compress(const BVEC_NAME &mask,
                                             const BVEC_NAME &a) {
    const __m512i c_i1 = _mm512_set1_epi32(1);
    __m512i a_int_vec = _mm512_mask_blend_epi32(a.val_, _mm512_setzero_epi32(),
                                                c_i1);
    __m512i compressed = _mm512_mask_compress_epi32(_mm512_undefined_epi32(),
                                                    mask.val_, a_int_vec);
    return _mm512_cmpeq_epi32_mask(compressed, c_i1);
  }
  VEC_INLINE static BVEC_NAME mask_expand(const BVEC_NAME &src,
                                          const BVEC_NAME &mask,
                                          const BVEC_NAME &a) {
    const __m512i c_i1 = _mm512_set1_epi32(1);
    __m512i a_int_vec = _mm512_mask_blend_epi32(a.val_, _mm512_setzero_epi32(),
                                                c_i1);
    __m512i src_int_vec = _mm512_mask_blend_epi32(src.val_,
                                                  _mm512_setzero_epi32(), c_i1);
    __m512i compressed = _mm512_mask_expand_epi32(src_int_vec, mask.val_,
                                                  a_int_vec);
    return _mm512_cmpeq_epi32_mask(compressed, c_i1);
  }
  VEC_INLINE static BVEC_NAME full() {
    return static_cast<FVEC_MASK_T>(0xFFFF);
  }
  VEC_INLINE static BVEC_NAME empty() {
    return 0;
  }
  VEC_INLINE static BVEC_NAME only(int n) {
    return full().val_ >> (FVEC_LEN - n);
  }
  VEC_INLINE static BVEC_NAME after(int n) {
    return full().val_ << n;
  }
  VEC_INLINE static BVEC_NAME onlyafter(int only, int after) {
    return (full().val_ >> (FVEC_LEN - only)) << after;
  }
  VEC_INLINE static int popcnt(const BVEC_NAME &a) {
    return _popcnt32(a.val_);
  }
  VEC_INLINE static bool test_all_unset(const BVEC_NAME &a) {
    return _mm512_kortestz(a.val_, a.val_);
  }
  VEC_INLINE static bool test_any_set(const BVEC_NAME &a) {
    return ! test_all_unset(a);
  }
  VEC_INLINE static bool test_at(const BVEC_NAME &a, int i) {
    assert(i < FVEC_LEN);
    return a.val_ & (1 << i);
  }
  VEC_INLINE BVEC_NAME operator &(const BVEC_NAME &b) const {
    return _mm512_kand(val_, b.val_);
  }
  VEC_INLINE BVEC_NAME operator |(const BVEC_NAME &b) const {
    return _mm512_kor(val_, b.val_);
  }
  VEC_INLINE BVEC_NAME operator ~() const {
    return _mm512_knot(val_);
  }
};

class IVEC_NAME {
  friend class FVEC_NAME;
  friend class AVEC_NAME;
# if FVEC_LEN==16
  friend class avec16pd;
# endif
  __m512i val_;
  VEC_INLINE IVEC_NAME(const __m512i &v) : val_(v) {}
public:
  static const int VL = 16;
  VEC_INLINE IVEC_NAME() {}

  #define IVEC_MASK_BINFN_B(the_name)                                \
    VEC_INLINE static BVEC_NAME the_name(const IVEC_NAME &a,         \
      const IVEC_NAME &b) {                                          \
      return _mm512_##the_name##_epi32_mask(a.val_, b.val_);         \
    }                                                                \
    VEC_INLINE static BVEC_NAME mask_##the_name(                        \
                                                const BVEC_NAME &mask,  \
                                                  const IVEC_NAME &a,   \
                                                  const IVEC_NAME &b    \
                                                  ) {                   \
      return _mm512_mask_##the_name##_epi32_mask(                       \
      mask.val_, a.val_, b.val_);                                       \
    }
  IVEC_MASK_BINFN_B(cmpeq)
  IVEC_MASK_BINFN_B(cmplt)
  IVEC_MASK_BINFN_B(cmpneq)
  IVEC_MASK_BINFN_B(cmpgt)

  #define IVEC_MASK_BINFN_I(the_name)                                   \
    VEC_INLINE static IVEC_NAME mask_##the_name(                        \
        const IVEC_NAME &src, const BVEC_NAME &mask,                    \
        const IVEC_NAME &a, const IVEC_NAME &b                          \
    ) {                                                                 \
       return _mm512_mask_##the_name##_epi32(                           \
        src.val_, mask.val_, a.val_, b.val_);                           \
    }
  IVEC_MASK_BINFN_I(add)
  VEC_INLINE static IVEC_NAME mask_blend(
      const BVEC_NAME &mask, const IVEC_NAME &a, const IVEC_NAME &b
  ) {
    return _mm512_mask_blend_epi32(mask.val_, a.val_, b.val_);
  }

  #define IVEC_BINFN_I(the_name)                                     \
    VEC_INLINE static IVEC_NAME the_name(const IVEC_NAME &a,         \
                                         const IVEC_NAME &b) {       \
      return _mm512_##the_name##_epi32(a.val_, b.val_);              \
    }
  IVEC_BINFN_I(mullo)
  IVEC_BINFN_I(srlv)
  VEC_INLINE static IVEC_NAME the_and(const IVEC_NAME &a, const IVEC_NAME &b) {
    return _mm512_and_epi32(a.val_, b.val_);
  }

  VEC_INLINE static IVEC_NAME mask_expand(
      const IVEC_NAME &src, const BVEC_NAME &a, const IVEC_NAME &b
  ) {
    return _mm512_mask_expand_epi32(src.val_,
      a.val_, b.val_);
  }
  VEC_INLINE static IVEC_NAME masku_compress(
      const BVEC_NAME &a, const IVEC_NAME &b
  ) {
    return _mm512_mask_compress_epi32(_mm512_undefined_epi32(), a.val_, b.val_);
  }

  VEC_INLINE static int at(const IVEC_NAME &a, int b) {
    int data[16] __attribute__((aligned(64)));
    _mm512_store_epi32(data, a.val_);
    return data[b];
  }

  VEC_INLINE static IVEC_NAME load(const int * src) {
    return _mm512_load_epi32(src);
  }
  VEC_INLINE static IVEC_NAME mask_loadu(const BVEC_NAME &mask,
                                         const int * src) {
    assert((mask.val_ & (mask.val_ + 1)) == 0);
    assert(mask.val_ <= BVEC_NAME::full().val_);
    return _mm512_mask_loadu_epi32(_mm512_undefined_epi32(), mask.val_, src);
  }
  VEC_INLINE static IVEC_NAME maskz_loadu(const BVEC_NAME &mask,
                                          const int * src) {
    assert((mask.val_ & (mask.val_ + 1)) == 0);
    assert(mask.val_ <= BVEC_NAME::full().val_);
    return _mm512_maskz_loadu_epi32(mask.val_, src);
  }
  VEC_INLINE static void mask_storeu(const BVEC_NAME &mask, int * dest,
    const IVEC_NAME &src) {
    assert((mask.val_ & (mask.val_ + 1)) == 0);
    assert(mask.val_ <= BVEC_NAME::full().val_);
    _mm512_mask_storeu_epi32(dest, mask.val_, src.val_);
  }
  VEC_INLINE static void store(int * dest, const IVEC_NAME &src) {
    _mm512_store_epi32(dest, src.val_);
  }

  VEC_INLINE static IVEC_NAME mask_gather(
      const IVEC_NAME &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const int * mem, const int scale
  ) {
    assert(mask.val_ <= BVEC_NAME::full().val_);
    assert(scale == sizeof(int));
    return _mm512_mask_i32gather_epi32(src.val_, mask.val_, idx.val_, mem,
      sizeof(int));
  }
  VEC_INLINE static void mask_i32scatter(
      int * mem, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const IVEC_NAME &a, const int scale
  ) {
    assert(mask.val_ <= BVEC_NAME::full().val_);
    assert(scale == sizeof(int));
    _mm512_mask_i32scatter_epi32(mem, mask.val_, idx.val_, a.val_, sizeof(int));
  }

  VEC_INLINE static void mask_compressstore(const BVEC_NAME &mask, int * dest,
    const IVEC_NAME &src) {
    _mm512_mask_compressstoreu_epi32(dest, mask.val_, src.val_);
  }

  VEC_INLINE static IVEC_NAME set1(int i) {
    return _mm512_set1_epi32(i);
  }
  VEC_INLINE static IVEC_NAME setzero() {
    return _mm512_setzero_epi32();
  }
  VEC_INLINE static IVEC_NAME undefined() {
    return _mm512_undefined_epi32();
  }

  VEC_INLINE IVEC_NAME operator +(const IVEC_NAME &b) const {
    return _mm512_add_epi32(this->val_, b.val_);
  }
  VEC_INLINE static void print(const char * str, const IVEC_NAME &a) {
    int data[8] __attribute__((aligned(32)));
    store(data, a);
    printf("%s:", str);
    for (int i = 0; i < FVEC_LEN; i++) {
      printf(" %d", data[i]);
    }
    printf("\n");
  }
};

class FVEC_NAME {
  friend class AVEC_NAME;
#if FVEC_LEN==16
  friend class avec16pd;
#endif
  FVEC_VEC_T val_;
  VEC_INLINE FVEC_NAME(const FVEC_VEC_T &v) : val_(v) {}
public:
  static const int VL = FVEC_LEN;
  VEC_INLINE FVEC_NAME() {}
  VEC_INLINE static FVEC_SCAL_T at(const FVEC_NAME &a, int i) {
    assert(i < FVEC_LEN);
    FVEC_SCAL_T data[FVEC_LEN] __attribute__((aligned(64)));
    FVEC_SUFFIX(_mm512_store_)(data, a.val_);
    return data[i];
  }
  VEC_INLINE static bool fast_compress() { return true; }

  #define FVEC_MASK_BINFN_B(the_name)                                \
    VEC_INLINE static BVEC_NAME the_name(const FVEC_NAME &a,         \
                                         const FVEC_NAME &b) {       \
      return FVEC_SUFFIX_MASK(_mm512_##the_name##_)(a.val_, b.val_); \
    }                                                                \
    VEC_INLINE static BVEC_NAME mask_##the_name(                     \
        const BVEC_NAME &mask,                                       \
        const FVEC_NAME &a, const FVEC_NAME &b                       \
    ) {                                                              \
      return FVEC_SUFFIX_MASK(_mm512_mask_##the_name##_)(            \
        mask.val_, a.val_, b.val_);                                  \
    }
  FVEC_MASK_BINFN_B(cmple)
  FVEC_MASK_BINFN_B(cmplt)
  FVEC_MASK_BINFN_B(cmpneq)
  FVEC_MASK_BINFN_B(cmpnle)
  FVEC_MASK_BINFN_B(cmpnlt)

  #define FVEC_UNFN_F(the_name)                                      \
    VEC_INLINE static FVEC_NAME the_name(const FVEC_NAME &a) {       \
      return FVEC_SUFFIX(_mm512_##the_name##_)(a.val_);              \
    }
  FVEC_UNFN_F(abs)
  FVEC_UNFN_F(exp)
  FVEC_UNFN_F(invsqrt)
  FVEC_UNFN_F(recip)
  FVEC_UNFN_F(sqrt)

  #define FVEC_MASK_UNFN_F(the_name)                                 \
    VEC_INLINE static FVEC_NAME mask_##the_name(                     \
        const FVEC_NAME &src, const BVEC_NAME &mask,                 \
        const FVEC_NAME &a                                           \
    ) {                                                              \
      return FVEC_SUFFIX(_mm512_mask_##the_name##_)(                 \
        src.val_, mask.val_, a.val_);                                \
    }
  FVEC_MASK_UNFN_F(cos)
  FVEC_MASK_UNFN_F(recip)
  FVEC_MASK_UNFN_F(sqrt)

  #define FVEC_BINFN_F(the_name)                                     \
    VEC_INLINE static FVEC_NAME the_name(const FVEC_NAME &a,         \
                                         const FVEC_NAME &b) {       \
      return FVEC_SUFFIX(_mm512_##the_name##_)(a.val_, b.val_);      \
    }
  FVEC_BINFN_F(max)
  FVEC_BINFN_F(min)

  #define FVEC_MASK_BINFN_F(the_name)                                \
    VEC_INLINE static FVEC_NAME mask_##the_name(                     \
        const FVEC_NAME &src, const BVEC_NAME &mask,                 \
        const FVEC_NAME &a, const FVEC_NAME &b                       \
    ) {                                                              \
      return FVEC_SUFFIX(_mm512_mask_##the_name##_)(                 \
        src.val_, mask.val_, a.val_, b.val_);                        \
    }
  FVEC_MASK_BINFN_F(add)
  FVEC_MASK_BINFN_F(div)
  FVEC_MASK_BINFN_F(mul)
  FVEC_MASK_BINFN_F(sub)
  VEC_INLINE static FVEC_NAME mask_blend(
      const BVEC_NAME &mask, const FVEC_NAME &a, const FVEC_NAME &b
  ) {
    return FVEC_SUFFIX(_mm512_mask_blend_)(mask.val_, a.val_, b.val_);
  }

  VEC_INLINE static FVEC_NAME mask_expand(
      const FVEC_NAME &src, const BVEC_NAME &a, const FVEC_NAME &b
  ) {
    return FVEC_SUFFIX(_mm512_mask_expand_)(src.val_,
      a.val_, b.val_);
  }
  VEC_INLINE static FVEC_NAME masku_compress(
      const BVEC_NAME &a, const FVEC_NAME &b
  ) {
    return FVEC_SUFFIX(_mm512_mask_compress_)(FVEC_SUFFIX(_mm512_undefined_)(),
                                                a.val_, b.val_);
  }

  VEC_INLINE static FVEC_NAME set1(const FVEC_SCAL_T &a) {
    return FVEC_SUFFIX(_mm512_set1_)(a);
  }
  VEC_INLINE static FVEC_NAME setzero() {
    return FVEC_SUFFIX(_mm512_setzero_)();
  }
  VEC_INLINE static FVEC_NAME undefined() {
    return FVEC_SUFFIX(_mm512_undefined_)();
  }

  VEC_INLINE static FVEC_NAME load(const FVEC_SCAL_T *mem) {
    return FVEC_SUFFIX(_mm512_load_)(mem);
  }
  VEC_INLINE static void mask_storeu(const BVEC_NAME &mask, FVEC_SCAL_T * dest,
                                       const FVEC_NAME &a) {
    FVEC_SUFFIX(_mm512_mask_storeu_)(dest, mask.val_, a.val_);
  }
  VEC_INLINE static void store(FVEC_SCAL_T * dest, const FVEC_NAME &a) {
    FVEC_SUFFIX(_mm512_store_)(dest, a.val_);
  }

  VEC_INLINE static FVEC_NAME gather(const IVEC_NAME &idx,
                                     const FVEC_SCAL_T * mem,
                                     const int scale) {
    assert(scale == sizeof(FVEC_SCAL_T));
#   if FVEC_LEN==8
    return FVEC_SUFFIX(_mm512_i32logather_)(idx.val_, mem, sizeof(FVEC_SCAL_T));
#   else
    return FVEC_SUFFIX(_mm512_i32gather_)(idx.val_, mem, sizeof(FVEC_SCAL_T));
#   endif
  }
  VEC_INLINE static FVEC_NAME mask_gather(
      const FVEC_NAME &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const FVEC_SCAL_T * mem, const int scale
  ) {
    assert(scale == sizeof(FVEC_SCAL_T));
#   if FVEC_LEN==8
    return FVEC_SUFFIX(_mm512_mask_i32logather_)(src.val_, mask.val_, idx.val_,
                       mem, sizeof(FVEC_SCAL_T));
#   else
    return FVEC_SUFFIX(_mm512_mask_i32gather_)(src.val_, mask.val_, idx.val_,
                       mem, sizeof(FVEC_SCAL_T));
#   endif
  }

  VEC_INLINE static void gather_3_adjacent(const IVEC_NAME &idx,
                                           const FVEC_SCAL_T * mem,
                                           const int scale,
                                           FVEC_NAME * out_0,
                                           FVEC_NAME * out_1,
                                           FVEC_NAME * out_2) {
    assert(scale == sizeof(FVEC_SCAL_T));
    *out_0 = FVEC_NAME::gather(idx, mem + 0, scale);
    *out_1 = FVEC_NAME::gather(idx, mem + 1, scale);
    *out_2 = FVEC_NAME::gather(idx, mem + 2, scale);
  }
  VEC_INLINE static void gather_4_adjacent(const IVEC_NAME &idx,
                                           const FVEC_SCAL_T * mem,
                                           const int scale, FVEC_NAME * out_0,
                                           FVEC_NAME * out_1,
                                           FVEC_NAME * out_2,
                                           FVEC_NAME * out_3) {
    assert(scale == sizeof(FVEC_SCAL_T));
    *out_0 = FVEC_NAME::gather(idx, mem + 0, scale);
    *out_1 = FVEC_NAME::gather(idx, mem + 1, scale);
    *out_2 = FVEC_NAME::gather(idx, mem + 2, scale);
    *out_3 = FVEC_NAME::gather(idx, mem + 3, scale);
  }

  VEC_INLINE static FVEC_SCAL_T mask_reduce_add(const BVEC_NAME &mask,
                                                const FVEC_NAME &a) {
    return FVEC_SUFFIX(_mm512_mask_reduce_add_)(mask.val_, a.val_);
  }
  VEC_INLINE static FVEC_SCAL_T reduce_add(const FVEC_NAME &a) {
    return FVEC_SUFFIX(_mm512_reduce_add_)(a.val_);
  }

  VEC_INLINE static IVEC_NAME unpackloepi32(const FVEC_NAME &a) {
#   if FVEC_LEN==8
    return _mm512_maskz_compress_epi32(0x5555, _mm512_castpd_si512(a.val_));
#   else
    return _mm512_castps_si512(a.val_);
#   endif
  }

  VEC_INLINE static FVEC_NAME mask_sincos(
      FVEC_NAME * cos, const FVEC_NAME &src_a, const FVEC_NAME &src_b,
      const BVEC_NAME &mask, const FVEC_NAME &arg
  ) {
    return FVEC_SUFFIX(_mm512_mask_sincos_)(&cos->val_, src_a.val_, src_b.val_,
      mask.val_, arg.val_);
  }

  #define FVEC_BINOP(the_sym, the_name)                                      \
    VEC_INLINE inline FVEC_NAME operator the_sym(const FVEC_NAME &b) const { \
    return FVEC_SUFFIX(_mm512_##the_name##_)(this->val_, b.val_);            \
  }
  FVEC_BINOP(+, add)
  FVEC_BINOP(-, sub)
  FVEC_BINOP(*, mul)
  FVEC_BINOP(/, div)

  VEC_INLINE static void gather_prefetch0(const IVEC_NAME &a, void * mem) {
    #ifdef __AVX512PF__
    _mm512_mask_prefetch_i32gather_ps(a.val_, BVEC_NAME::full().val_, mem,
      sizeof(FVEC_SCAL_T), _MM_HINT_T0);
    #endif
  }
};

class AVEC_NAME {
  FVEC_VEC_T val_;
  VEC_INLINE AVEC_NAME(const FVEC_VEC_T &a) : val_(a) {}
public:
  VEC_INLINE AVEC_NAME(const FVEC_NAME &a) : val_(a.val_) {}
  VEC_INLINE static AVEC_NAME undefined() {
    return FVEC_SUFFIX(_mm512_undefined_)();
  }
  VEC_INLINE static AVEC_NAME mask_gather(
      const AVEC_NAME &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const FVEC_SCAL_T * mem, const int scale
  ) {
    assert(scale == sizeof(FVEC_SCAL_T));
#   if FVEC_LEN==8
    return FVEC_SUFFIX(_mm512_mask_i32logather_)(src.val_, mask.val_, idx.val_,
                                                 mem, sizeof(FVEC_SCAL_T));
#   else
    return FVEC_SUFFIX(_mm512_mask_i32gather_)(src.val_, mask.val_, idx.val_,
                                               mem, sizeof(FVEC_SCAL_T));
#   endif
  }
  VEC_INLINE static void mask_i32loscatter(
      FVEC_SCAL_T * mem, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const AVEC_NAME &a, const int scale
  ) {
    assert(scale == sizeof(FVEC_SCAL_T));
#   if FVEC_LEN==8
    FVEC_SUFFIX(_mm512_mask_i32loscatter_)(mem, mask.val_, idx.val_, a.val_,
                                           sizeof(FVEC_SCAL_T));
#   else
    FVEC_SUFFIX(_mm512_mask_i32scatter_)(mem, mask.val_, idx.val_, a.val_,
                                         sizeof(FVEC_SCAL_T));
#   endif
  }

  #define AVEC_BINOP(the_sym, the_name)                                      \
    VEC_INLINE inline AVEC_NAME operator the_sym(const AVEC_NAME &b) const { \
    return FVEC_SUFFIX(_mm512_##the_name##_)(this->val_, b.val_);            \
  }
  AVEC_BINOP(-, sub)

  VEC_INLINE static void gather_prefetch0(const IVEC_NAME &a, void * mem) {
    _mm512_mask_prefetch_i32gather_ps(a.val_, BVEC_NAME::full().val_, mem,
                                      sizeof(FVEC_SCAL_T), _MM_HINT_T0);
  }
};

#if FVEC_LEN==16
class avec16pd {
  __m512d lo_, hi_;
  VEC_INLINE avec16pd(const __m512d &lo, const __m512d &hi) : lo_(lo), hi_(hi)
    {}
  VEC_INLINE static __mmask8 get_bvec_hi(__mmask16 a) {
    return a >> 8;
  }
  VEC_INLINE static __m512i get_ivec_hi(__m512i a) {
    return _mm512_permute4f128_epi32(a, _MM_PERM_BADC);
  }
public:
  VEC_INLINE avec16pd(const FVEC_NAME &a) {
    lo_ = _mm512_cvtpslo_pd(a.val_);
    hi_ = _mm512_cvtpslo_pd(_mm512_permute4f128_ps(a.val_, _MM_PERM_BADC));
  }
  VEC_INLINE static avec16pd undefined() {
    return avec16pd(_mm512_undefined_pd(), _mm512_undefined_pd());
  }
  VEC_INLINE static avec16pd mask_gather(
      const avec16pd &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const double * mem, const int scale
  ) {
    assert(scale == sizeof(double));
    __m512d lo = _mm512_mask_i32logather_pd(src.lo_, mask.val_, idx.val_, mem,
                                            sizeof(double));
    __m512d hi = _mm512_mask_i32logather_pd(src.hi_, get_bvec_hi(mask.val_),
                                            get_ivec_hi(idx.val_), mem,
                                            sizeof(double));
    return avec16pd(lo, hi);
  }
  VEC_INLINE static void mask_i32loscatter(
      double * mem, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const avec16pd &a, const int scale
  ) {
    assert(scale == sizeof(double));
    _mm512_mask_i32loscatter_pd(mem, mask.val_, idx.val_, a.lo_,
                                sizeof(double));
    _mm512_mask_i32loscatter_pd(mem, get_bvec_hi(mask.val_),
                                get_ivec_hi(idx.val_), a.hi_, sizeof(double));
  }

  #define AVEC2_BINOP(the_sym, the_name)                                    \
    VEC_INLINE inline avec16pd operator the_sym(const avec16pd &b) const {  \
    __m512d lo = _mm512_##the_name##_pd(this->lo_, b.lo_);                  \
    __m512d hi = _mm512_##the_name##_pd(this->hi_, b.hi_);                  \
    return avec16pd(lo, hi);                                                \
  }
  AVEC2_BINOP(-, sub)

  VEC_INLINE static void gather_prefetch0(const IVEC_NAME &a, void * mem) {
    _mm512_mask_prefetch_i32gather_ps(a.val_, BVEC_NAME::full().val_, mem,
                                      sizeof(double), _MM_HINT_T0);
  }
};
#endif

}


#ifdef FVEC_FIRST_PASS

template<typename flt_t, typename acc_t>
struct intr_types;

template<>
struct intr_types<double,double> {
  typedef mm512::fvec8pd fvec;
  typedef mm512::ivec8 ivec;
  typedef mm512::bvec8 bvec;
  typedef mm512::avec8pd avec;
};

template<>
struct intr_types<float,float> {
  typedef mm512::fvec16ps fvec;
  typedef mm512::ivec16 ivec;
  typedef mm512::bvec16 bvec;
  typedef mm512::avec16ps avec;
};

template<>
struct intr_types<float,double> {
  typedef mm512::fvec16ps fvec;
  typedef mm512::ivec16 ivec;
  typedef mm512::bvec16 bvec;
  typedef mm512::avec16pd avec;
};

#endif


#ifndef FVEC_FIRST_PASS
#  define FVEC_FIRST_PASS
#  include "intel_intrinsics_airebo.h"
#endif

#endif

#ifdef LMP_INTEL_AIREBO_256

#include <cassert>
#include <immintrin.h>
#include <stdint.h> // <cstdint> requires C++-11

#define VEC_INLINE __attribute__((always_inline))


#ifndef FVEC_FIRST_PASS
#  define FVEC_LEN 4
#  define FVEC_SUFFIX(a) a##pd
#  define FVEC_MASK_T __m256d
#  define FVEC_VEC_T __m256d
#  define FVEC_SCAL_T double
#  define IVEC_NAME ivec4
#  define FVEC_NAME fvec4pd
#  define BVEC_NAME bvec4
#  define AVEC_NAME avec4pd
#else
#  undef FVEC_LEN
#  undef FVEC_SUFFIX
#  undef FVEC_SUFFIX_MASK
#  undef FVEC_MASK_T
#  undef FVEC_VEC_T
#  undef FVEC_SCAL_T
#  undef IVEC_NAME
#  undef FVEC_NAME
#  undef BVEC_NAME
#  undef AVEC_NAME

#  define FVEC_LEN 8
#  define FVEC_SUFFIX(a) a##ps
#  define FVEC_MASK_T __m256
#  define FVEC_VEC_T __m256
#  define FVEC_SCAL_T float
#  define IVEC_NAME ivec8
#  define FVEC_NAME fvec8ps
#  define BVEC_NAME bvec8
#  define AVEC_NAME avec8ps
#endif



namespace mm256 {

//#define __AVX2__ __AVX2__

#if !defined(__AVX2__) && !defined(FVEC_FIRST_PASS)

#define IVEC_EM_BIN(op) \
  __m128i a_lo = _mm256_castsi256_si128(a);  \
  __m128i b_lo = _mm256_castsi256_si128(b);  \
  __m128i a_hi = _mm256_extractf128_si256(a, 1);  \
  __m128i b_hi = _mm256_extractf128_si256(b, 1);  \
  __m128i c_lo = op(a_lo, b_lo); \
  __m128i c_hi = op(a_hi, b_hi); \
  __m256i ret = _mm256_setr_m128i(c_lo, c_hi); \
  return ret;

VEC_INLINE inline __m256i _cm256_add_epi32(const __m256i &a, const __m256i &b) {
  IVEC_EM_BIN(_mm_add_epi32)
}

VEC_INLINE inline __m256i _cm256_and_si256(const __m256i &a, const __m256i &b) {
  IVEC_EM_BIN(_mm_and_si128)
}

VEC_INLINE inline __m256i _cm256_andnot_si256(const __m256i &a,
                                              const __m256i &b) {
  IVEC_EM_BIN(_mm_andnot_si128)
}

VEC_INLINE inline __m256i _cm256_cmpeq_epi32(const __m256i &a,
                                             const __m256i &b) {
  IVEC_EM_BIN(_mm_cmpeq_epi32)
}

VEC_INLINE inline __m256i _cm256_cmpgt_epi32(const __m256i &a,
                                             const __m256i &b) {
  IVEC_EM_BIN(_mm_cmpgt_epi32)
}

VEC_INLINE inline __m256i _cm256_cvtepu8_epi32(const __m128i &a) {
  __m128i a_hi = _mm_castps_si128(_mm_permute_ps(_mm_castsi128_ps(a), 1));
  __m128i c_lo = _mm_cvtepu8_epi32(a);
  __m128i c_hi = _mm_cvtepu8_epi32(a_hi);
  __m256i ret = _mm256_setr_m128i(c_lo, c_hi);
  return ret;

}

#define IVEC_EM_SCAL(op)                       \
  int buf_a[8] __attribute__((aligned(32)));   \
  int buf_b[8] __attribute__((aligned(32)));   \
  int dest[8] __attribute__((aligned(32)));    \
  _mm256_store_si256((__m256i*)buf_a, a);      \
  _mm256_store_si256((__m256i*)buf_b, b);      \
  for (int i = 0; i < 8; i++) {                \
    dest[i] = op;                              \
  }                                            \
  return _mm256_load_si256((__m256i*) dest);

VEC_INLINE inline __m256i _cm256_permutevar8x32_epi32(const __m256i &a,
                                                      const __m256i &b) {
  IVEC_EM_SCAL(buf_a[buf_b[i]])
}

VEC_INLINE inline __m256i _cm256_mullo_epi32(__m256i a, __m256i b) {
  IVEC_EM_BIN(_mm_mullo_epi32)
}

VEC_INLINE inline __m256i _cm256_srlv_epi32(__m256i a, __m256i b) {
  IVEC_EM_SCAL(buf_a[i] >> buf_b[i])
}


VEC_INLINE inline __m256 _cm256_permutevar8x32_ps(const __m256 &a,
                                                  const __m256i &b) {
  return _mm256_castsi256_ps(_cm256_permutevar8x32_epi32(_mm256_castps_si256(a),
                                                         b));
}

VEC_INLINE inline __m128i _cm_maskload_epi32(int const * mem, __m128i mask) {
  return _mm_castps_si128(_mm_maskload_ps((float const *) mem, mask));
}

VEC_INLINE inline __m256i _cm256_maskload_epi32(int const * mem, __m256i mask) {
  __m128i a_lo = _mm256_castsi256_si128(mask);
  __m128i a_hi = _mm256_extractf128_si256(mask, 1);
  __m128i c_lo = _cm_maskload_epi32(mem, a_lo);
  __m128i c_hi = _cm_maskload_epi32(mem + 4, a_hi);
  __m256i ret = _mm256_setr_m128i(c_lo, c_hi);
  return ret;
}


VEC_INLINE inline __m256i _cm256_mask_i32gather_epi32(__m256i src,
                                                      int const * base_addr,
                                                      __m256i index,
                                                      __m256i mask,
                                                      const int scale) {
  assert(scale == sizeof(int));
  int buf_index[8] __attribute__((aligned(32)));
  int buf_mask[8] __attribute__((aligned(32)));
  int dest[8] __attribute__((aligned(32)));
  _mm256_store_si256((__m256i*)dest, src);
  _mm256_store_si256((__m256i*)buf_index, index);
  _mm256_store_si256((__m256i*)buf_mask, mask);
  for (int i = 0; i < 8; i++) {
    if (buf_mask[i]) dest[i] = base_addr[buf_index[i]];
  }
  return _mm256_load_si256((__m256i*) dest);
}

VEC_INLINE inline __m256 _cm256_mask_i32gather_ps(__m256 src,
                                                  float const * base_addr,
                                                  __m256i index, __m256 mask,
                                                  const int scale) {
  return _mm256_castsi256_ps(_cm256_mask_i32gather_epi32(
    _mm256_castps_si256(src), (const int *) base_addr, index,
    _mm256_castps_si256(mask), scale));
}

VEC_INLINE inline __m256d _cm256_mask_i32gather_pd(__m256d src,
                                                   double const * base_addr,
                                                   __m128i index, __m256d mask,
                                                   const int scale) {
  assert(scale == sizeof(double));
  int buf_index[4] __attribute__((aligned(32)));
  int buf_mask[8] __attribute__((aligned(32)));
  double dest[4] __attribute__((aligned(32)));
  _mm256_store_pd(dest, src);
  _mm_store_si128((__m128i*)buf_index, index);
  _mm256_store_si256((__m256i*)buf_mask, _mm256_castpd_si256(mask));
  for (int i = 0; i < 4; i++) {
    if (buf_mask[2*i]) dest[i] = base_addr[buf_index[i]];
  }
  return _mm256_load_pd(dest);
}

VEC_INLINE inline __m256i _cm256_i32gather_epi32(int const * base_addr,
                                                 __m256i index,
                                                 const int scale) {
  assert(scale == sizeof(int));
  int buf_index[8] __attribute__((aligned(32)));
  int dest[8] __attribute__((aligned(32)));
  _mm256_store_si256((__m256i*)buf_index, index);
  for (int i = 0; i < 8; i++) {
    dest[i] = base_addr[buf_index[i]];
  }
  return _mm256_load_si256((__m256i*) dest);
}

VEC_INLINE inline __m256 _cm256_i32gather_ps(float const * base_addr,
                                             __m256i index, const int scale) {
  return _mm256_castsi256_ps(_cm256_i32gather_epi32((const int *) base_addr,
                                                    index, scale));
}

VEC_INLINE inline __m256d _cm256_i32gather_pd(double const * base_addr,
                                              __m128i index, const int scale) {
  assert(scale == sizeof(double));
  int buf_index[4] __attribute__((aligned(32)));
  double dest[4] __attribute__((aligned(32)));
  _mm_store_si128((__m128i*)buf_index, index);
  for (int i = 0; i < 4; i++) {
    dest[i] = base_addr[buf_index[i]];
  }
  return _mm256_load_pd(dest);
}

VEC_INLINE inline uint64_t _cdep_u64(uint64_t tmp, uint64_t mask) {
  uint64_t dst = 0;
  uint64_t k = 0;
  const uint64_t one = 1;
  const uint64_t zero = 0;
  for (uint64_t m = 0; m < 64; m++) {
    if (mask & (one << m)) {
      dst |= static_cast<uint64_t>((tmp & (one << k)) != zero) << m;
      k += 1;
    }
  }
  return dst;
}

VEC_INLINE inline uint64_t _cext_u64(uint64_t tmp, uint64_t mask) {
  uint64_t dst = 0;
  uint64_t k = 0;
  const uint64_t one = 1;
  const uint64_t zero = 0;
  for (uint64_t m = 0; m < 64; m++) {
    if (mask & (one << m)) {
      dst |= static_cast<uint64_t>((tmp & (one << m)) != zero) << k;
      k += 1;
    }
  }
  return dst;
}

#define _mm256_add_epi32 _cm256_add_epi32
#define _mm256_and_si256 _cm256_and_si256
#define _mm256_andnot_si256 _cm256_andnot_si256
#define _mm256_cmpeq_epi32 _cm256_cmpeq_epi32
#define _mm256_cmpgt_epi32 _cm256_cmpgt_epi32
#define _mm256_permutevar8x32_epi32 _cm256_permutevar8x32_epi32
#define _mm256_permutevar8x32_ps _cm256_permutevar8x32_ps
#define _mm_maskload_epi32 _cm_maskload_epi32
#define _mm256_maskload_epi32 _cm256_maskload_epi32
#define _mm256_mullo_epi32 _cm256_mullo_epi32
#define _mm256_srlv_epi32 _cm256_srlv_epi32
#define _mm256_mask_i32gather_epi32 _cm256_mask_i32gather_epi32
#define _mm256_mask_i32gather_pd _cm256_mask_i32gather_pd
#define _mm256_mask_i32gather_ps _cm256_mask_i32gather_ps
#define _mm256_i32gather_epi32 _cm256_i32gather_epi32
#define _mm256_i32gather_pd _cm256_i32gather_pd
#define _mm256_i32gather_ps _cm256_i32gather_ps
#define _pdep_u64 _cdep_u64
#define _pext_u64 _cext_u64
#define _mm256_cvtepu8_epi32 _cm256_cvtepu8_epi32

#endif

#ifndef FVEC_FIRST_PASS

VEC_INLINE inline __m256 _mm256_compress_ps(__m256 mask, __m256 a) {
# ifdef __AVX2__
  uint64_t expanded_mask = _pdep_u64(_mm256_movemask_ps(mask),
                                     0x0101010101010101);
  // unpack each bit to a byte
  expanded_mask *= 0xFF;   // mask |= mask<<1 | mask<<2 | ... | mask<<7;
  // the identity shuffle for vpermps, packed to one index per byte
  const uint64_t identity_indices = 0x0706050403020100;
  uint64_t wanted_indices = _pext_u64(identity_indices, expanded_mask);

  __m128i bytevec = _mm_cvtsi64_si128(wanted_indices);
  __m256i shufmask = _mm256_cvtepu8_epi32(bytevec);

  return _mm256_permutevar8x32_ps(a, shufmask);
# else
  int mask_buf[8] __attribute__((aligned(32)));
  float a_buf[8] __attribute__((aligned(32)));
  float dst_buf[8] __attribute__((aligned(32)));
  _mm256_store_si256((__m256i*) mask_buf, _mm256_castps_si256(mask));
  _mm256_store_ps(a_buf, a);
  int k = 0;
  for (int i = 0; i < 8; i++) {
    if (mask_buf[i]) {
      dst_buf[k++] = a_buf[i];
    }
  }
  return _mm256_load_ps(dst_buf);
# endif
}
VEC_INLINE inline __m256 _mm256_expand_ps(__m256 mask, __m256 a) {
# ifdef __AVX2__
  uint64_t expanded_mask = _pdep_u64(_mm256_movemask_ps(mask),
                                     0x0101010101010101);
  expanded_mask *= 0xFF;
  const uint64_t identity_indices = 0x0706050403020100;
  uint64_t wanted_indices = _pdep_u64(identity_indices, expanded_mask);
  __m128i bytevec = _mm_cvtsi64_si128(wanted_indices);
  __m256i shufmask = _mm256_cvtepu8_epi32(bytevec);
  return _mm256_permutevar8x32_ps(a, shufmask);
# else
  int mask_buf[8] __attribute__((aligned(32)));
  float a_buf[8] __attribute__((aligned(32)));
  float dst_buf[8] __attribute__((aligned(32))) = {0};
  _mm256_store_si256((__m256i*) mask_buf, _mm256_castps_si256(mask));
  _mm256_store_ps(a_buf, a);
  int k = 0;
  for (int i = 0; i < 8; i++) {
    if (mask_buf[i]) {
      dst_buf[i] = a_buf[k++];
    }
  }
  return _mm256_load_ps(dst_buf);
# endif
}

VEC_INLINE inline __m256d _mm256_compress_pd(__m256d mask, __m256d a) {
  return _mm256_castps_pd(_mm256_compress_ps(_mm256_castpd_ps(mask),
                                             _mm256_castpd_ps(a)));
}
VEC_INLINE inline __m256d _mm256_expand_pd(__m256d mask, __m256d a) {
  return _mm256_castps_pd(_mm256_expand_ps(_mm256_castpd_ps(mask),
                                           _mm256_castpd_ps(a)));
}
#endif


class FVEC_NAME;
class IVEC_NAME;
class AVEC_NAME;
class BVEC_NAME {
  friend class FVEC_NAME;
  friend class IVEC_NAME;
  friend class AVEC_NAME;
# if FVEC_LEN==8
  friend class avec8pd;
# endif
  FVEC_MASK_T val_;
  VEC_INLINE BVEC_NAME(const FVEC_MASK_T &v) : val_(v) {}
  VEC_INLINE BVEC_NAME(const __m256i &v) : val_(FVEC_SUFFIX(_mm256_castsi256_)
                                                (v)) {}
public:
  VEC_INLINE BVEC_NAME() {}
  VEC_INLINE static BVEC_NAME kand(const BVEC_NAME &a, const BVEC_NAME &b) {
    return FVEC_SUFFIX(_mm256_and_)(a.val_, b.val_);
  }
  VEC_INLINE static BVEC_NAME kandn(const BVEC_NAME &a, const BVEC_NAME &b) {
    return FVEC_SUFFIX(_mm256_andnot_)(a.val_, b.val_);
  }
  VEC_INLINE static BVEC_NAME masku_compress(const BVEC_NAME &mask,
                                             const BVEC_NAME &a) {
    return FVEC_SUFFIX(_mm256_compress_)(mask.val_, a.val_);
  }
  VEC_INLINE static BVEC_NAME mask_expand(const BVEC_NAME &src,
                                          const BVEC_NAME &mask,
                                          const BVEC_NAME &a) {
    FVEC_MASK_T ret = FVEC_SUFFIX(_mm256_expand_)(mask.val_, a.val_);
    ret = FVEC_SUFFIX(_mm256_and_)(mask.val_, ret);
    ret = FVEC_SUFFIX(_mm256_or_)(ret, FVEC_SUFFIX(_mm256_andnot_)
                                  (mask.val_, src.val_));
    return ret;
  }
  VEC_INLINE static BVEC_NAME full() {
    __m256i a = _mm256_undefined_si256();
    return FVEC_SUFFIX(_mm256_castsi256_)(_mm256_cmpeq_epi32(a, a));
  }
  VEC_INLINE static BVEC_NAME empty() {
    return FVEC_SUFFIX(_mm256_setzero_)();
  }
  VEC_INLINE static BVEC_NAME only(int n) {
    static const unsigned int FULL_ps = (unsigned int) -1;
    static const unsigned int LUT_ps[9][8] = {
      {0, 0, 0, 0, 0, 0, 0, 0},
      {FULL_ps, 0, 0, 0, 0, 0, 0, 0},
      {FULL_ps, FULL_ps, 0, 0, 0, 0, 0, 0},
      {FULL_ps, FULL_ps, FULL_ps, 0, 0, 0, 0, 0},
      {FULL_ps, FULL_ps, FULL_ps, FULL_ps, 0, 0, 0, 0},
      {FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, 0, 0, 0},
      {FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, 0, 0},
      {FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, 0},
      {FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps},
    };
    static const unsigned long long FULL_pd = (unsigned long long) -1;
    static const unsigned long long LUT_pd[5][4] = {
      {0, 0, 0, 0},
      {FULL_pd, 0, 0, 0},
      {FULL_pd, FULL_pd, 0, 0},
      {FULL_pd, FULL_pd, FULL_pd, 0},
      {FULL_pd, FULL_pd, FULL_pd, FULL_pd},
    };
    return FVEC_SUFFIX(_mm256_load_)((const FVEC_SCAL_T*) FVEC_SUFFIX(LUT_)[n]);
  }
  VEC_INLINE static BVEC_NAME after(int n) {
    static const unsigned int FULL_ps = (unsigned int) -1;
    static const unsigned int LUT_ps[9][8] = {
      {FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps},
      {0, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps},
      {0, 0, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps},
      {0, 0, 0, FULL_ps, FULL_ps, FULL_ps, FULL_ps, FULL_ps},
      {0, 0, 0, 0, FULL_ps, FULL_ps, FULL_ps, FULL_ps},
      {0, 0, 0, 0, 0, FULL_ps, FULL_ps, FULL_ps},
      {0, 0, 0, 0, 0, 0, FULL_ps, FULL_ps},
      {0, 0, 0, 0, 0, 0, 0, FULL_ps},
      {0, 0, 0, 0, 0, 0, 0, 0},
    };
    static const unsigned long long FULL_pd = (unsigned long long) -1;
    static const unsigned long long LUT_pd[5][4] = {
      {FULL_pd, FULL_pd, FULL_pd, FULL_pd},
      {0, FULL_pd, FULL_pd, FULL_pd},
      {0, 0, FULL_pd, FULL_pd},
      {0, 0, 0, FULL_pd},
      {0, 0, 0, 0},
    };
    return FVEC_SUFFIX(_mm256_load_)((const FVEC_SCAL_T*) FVEC_SUFFIX(LUT_)[n]);
  }
  VEC_INLINE static BVEC_NAME onlyafter(int only_, int after_) {
    return kand(after(after_), only(after_ + only_));
  }
  VEC_INLINE static int popcnt(const BVEC_NAME &a) {
    return _popcnt32(FVEC_SUFFIX(_mm256_movemask_)(a.val_));
  }
  VEC_INLINE static bool test_all_unset(const BVEC_NAME &a) {
    return FVEC_SUFFIX(_mm256_testz_)(a.val_, a.val_);
  }
  VEC_INLINE static bool test_any_set(const BVEC_NAME &a) {
    return ! test_all_unset(a);
  }
  VEC_INLINE static bool test_at(const BVEC_NAME &a, int i) {
    assert(i < FVEC_LEN);
    return FVEC_SUFFIX(_mm256_movemask_)(a.val_) & (1 << i);
  }
  VEC_INLINE BVEC_NAME operator &(const BVEC_NAME &b) const {
    return FVEC_SUFFIX(_mm256_and_)(val_, b.val_);
  }
  VEC_INLINE BVEC_NAME operator |(const BVEC_NAME &b) const {
    return FVEC_SUFFIX(_mm256_or_)(val_, b.val_);
  }
  VEC_INLINE BVEC_NAME operator ~() const {
    return FVEC_SUFFIX(_mm256_andnot_)(val_, full().val_);
  }
};

class IVEC_NAME {
  friend class FVEC_NAME;
  friend class AVEC_NAME;
# if FVEC_LEN==8
  friend class avec8pd;
# endif
  __m256i val_;
  VEC_INLINE IVEC_NAME(const __m256i &v) : val_(v) {}
  VEC_INLINE static __m256i to(const FVEC_VEC_T &a) {
#   if FVEC_LEN==4
    return _mm256_castpd_si256(a);
#   else
    return _mm256_castps_si256(a);
#   endif
  }
  VEC_INLINE static FVEC_VEC_T from(const __m256i &a) {
    return FVEC_SUFFIX(_mm256_castsi256_)(a);
  }
public:
  static const int VL = 8;
  VEC_INLINE IVEC_NAME() {}

  #define IVEC_MASK_BINFN_B(the_name)                                \
    VEC_INLINE static BVEC_NAME the_name(const IVEC_NAME &a,         \
                                         const IVEC_NAME &b) {       \
      return _mm256_##the_name##_epi32(a.val_, b.val_);              \
    }                                                                \
    VEC_INLINE static BVEC_NAME mask_##the_name(                     \
        const BVEC_NAME &mask,                                       \
        const IVEC_NAME &a, const IVEC_NAME &b                       \
    ) {                                                              \
      BVEC_NAME ret = _mm256_##the_name##_epi32(                     \
        a.val_, b.val_);                                             \
      return mask & ret;                                             \
    }
  IVEC_MASK_BINFN_B(cmpeq)
  IVEC_MASK_BINFN_B(cmpgt)

  VEC_INLINE static __m256i _mm256_cmplt_epi32(__m256i a, __m256i b) {
    __m256i le = _mm256_cmpgt_epi32(b, a);
    __m256i eq = _mm256_cmpeq_epi32(a, b);
    return _mm256_andnot_si256(eq, le);
  }

  VEC_INLINE static __m256i _mm256_cmpneq_epi32(__m256i a, __m256i b) {
    __m256i eq = _mm256_cmpeq_epi32(a, b);
    __m256i t = _mm256_undefined_si256();
    __m256i f = _mm256_cmpeq_epi32(t, t);
    return _mm256_andnot_si256(eq, f);
  }

  IVEC_MASK_BINFN_B(cmplt)
  IVEC_MASK_BINFN_B(cmpneq)
  #undef IVEC_MASK_BINFN_B

  VEC_INLINE static IVEC_NAME mask_blend(
      const BVEC_NAME &mask, const IVEC_NAME &a, const IVEC_NAME &b
  ) {
    return to(FVEC_SUFFIX(_mm256_blendv_)(from(a.val_), from(b.val_),
              mask.val_));
  }
  #define IVEC_MASK_BINFN_I(the_name)                                \
    VEC_INLINE static IVEC_NAME mask_##the_name(                     \
        const IVEC_NAME &src, const BVEC_NAME &mask,                 \
        const IVEC_NAME &a, const IVEC_NAME &b                       \
    ) {                                                              \
      IVEC_NAME ret = _mm256_##the_name##_epi32(                     \
                                                a.val_, b.val_);     \
        return mask_blend(mask, src, ret);                           \
    }
  IVEC_MASK_BINFN_I(add)
  #undef IVEC_MASK_BINFN_I

  #define IVEC_BINFN_I(the_name)                                     \
    VEC_INLINE static IVEC_NAME the_name(const IVEC_NAME &a,         \
                                         const IVEC_NAME &b) {       \
      return _mm256_##the_name##_epi32(a.val_, b.val_);              \
    }
  IVEC_BINFN_I(mullo)
  IVEC_BINFN_I(srlv)
  #undef IVEC_BINFN_I
  VEC_INLINE static IVEC_NAME the_and(const IVEC_NAME &a, const IVEC_NAME &b) {
    return _mm256_and_si256(a.val_, b.val_);
  }

  VEC_INLINE static IVEC_NAME masku_compress(const BVEC_NAME &mask,
                                             const IVEC_NAME &b) {
    return to(FVEC_SUFFIX(_mm256_compress_)(mask.val_, from(b.val_)));
  }
  VEC_INLINE static IVEC_NAME mask_expand(
      const IVEC_NAME &src, const BVEC_NAME &mask, const IVEC_NAME &b
  ) {
    FVEC_VEC_T ret = FVEC_SUFFIX(_mm256_expand_)(mask.val_, from(b.val_));
    ret = FVEC_SUFFIX(_mm256_and_)(mask.val_, ret);
    ret = FVEC_SUFFIX(_mm256_or_)(ret, FVEC_SUFFIX(_mm256_andnot_)
                                    (mask.val_, from(src.val_)));
    return to(ret);
  }

  VEC_INLINE static void store(int * dest, const IVEC_NAME &src) {
    _mm256_store_si256((__m256i*)dest, src.val_);
#   if FVEC_LEN==4
    dest[1] = dest[2];
    dest[2] = dest[4];
    dest[3] = dest[6];
#   endif
  }

  VEC_INLINE static int at(const IVEC_NAME &a, int b) {
    int data[8] __attribute__((aligned(32)));
    store(data, a);
    return data[b];
  }

  VEC_INLINE static void print(const char * str, const IVEC_NAME &a) {
    int data[8] __attribute__((aligned(32)));
    store(data, a);
    printf("%s:", str);
    for (int i = 0; i < FVEC_LEN; i++) {
      printf(" %d", data[i]);
    }
    printf("\n");
  }

  VEC_INLINE static IVEC_NAME maskz_loadu(const BVEC_NAME &mask,
                                          const int * src) {
    FVEC_VEC_T mask_val = mask.val_;
#   if FVEC_LEN==4
#    ifdef __AVX2__
    static const unsigned int mask_shuffle[8] __attribute__((aligned(32))) =
      {0, 2, 4, 6, 0, 0, 0, 0};
    __m256 m = _mm256_castpd_ps(mask_val);
    m = _mm256_permutevar8x32_ps(m, _mm256_load_si256((__m256i*)mask_shuffle));
    __m128i ret = _mm_maskload_epi32(src,
       _mm256_castsi256_si128(_mm256_castps_si256(m)));
    static const unsigned int load_shuffle[8] __attribute__((aligned(32))) =
      {0, 0, 1, 1, 2, 2, 3, 3};
    return _mm256_permutevar8x32_epi32(_mm256_castsi128_si256(ret),
      _mm256_load_si256((__m256i*)load_shuffle));
#    else
    int dest[8] __attribute__((aligned(32))) = {0};
    int mask_buf[8] __attribute__((aligned(32)));
    _mm256_store_pd((double*) mask_buf, mask.val_);
    for (int i = 0; i < 4; i++) {
      if (mask_buf[2*i]) {
        int val = src[i];
        dest[2*i+0] = val;
        dest[2*i+1] = val;
      }
    }
    return _mm256_load_si256((__m256i*) dest);
#    endif
#   else
    return _mm256_maskload_epi32(src, to(mask_val));
#   endif
  }

  VEC_INLINE static IVEC_NAME mask_gather(
      const IVEC_NAME &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const int * mem, const int scale
  ) {
    assert(scale == sizeof(int));
    return _mm256_mask_i32gather_epi32(src.val_, mem, idx.val_, to(mask.val_),
                                       sizeof(int));
  }

  VEC_INLINE static void mask_compressstore(const BVEC_NAME &mask, int * dest,
                                            const IVEC_NAME &src) {
    int buf[8] __attribute__((aligned(64)));
    const int stride = FVEC_LEN==4 ? 2 : 1;
    _mm256_store_si256((__m256i*)buf, src.val_);
    int mask_val = FVEC_SUFFIX(_mm256_movemask_)(mask.val_);
    int k = 0;
    #pragma unroll
    for (int i = 0; i < FVEC_LEN; i++) {
      if (mask_val & (1 << i))
        dest[k++] = buf[stride*i];
    }
  }

  VEC_INLINE static IVEC_NAME set1(int i) {
    return _mm256_set1_epi32(i);
  }
  VEC_INLINE static IVEC_NAME setzero() {
    return _mm256_setzero_si256();
  }
  VEC_INLINE static IVEC_NAME undefined() {
    return _mm256_undefined_si256();
  }

  VEC_INLINE IVEC_NAME operator +(const IVEC_NAME &b) const {
    return _mm256_add_epi32(this->val_, b.val_);
  }
};

class FVEC_NAME {
  friend class AVEC_NAME;
#if FVEC_LEN==8
  friend class avec8pd;
#endif
  FVEC_VEC_T val_;
  VEC_INLINE FVEC_NAME(const FVEC_VEC_T &v) : val_(v) {}
public:
  static const int VL = FVEC_LEN;
# if defined(__AVX2__) || defined(__MIC__) || defined(__AVX512F__)
  VEC_INLINE static bool fast_compress() { return true; }
# else
  VEC_INLINE static bool fast_compress() { return false; }
# endif
  VEC_INLINE FVEC_NAME() {}
  VEC_INLINE static FVEC_SCAL_T at(const FVEC_NAME &a, int i) {
    assert(i < FVEC_LEN);
    FVEC_SCAL_T data[FVEC_LEN] __attribute__((aligned(64)));
    FVEC_SUFFIX(_mm256_store_)(data, a.val_);
    return data[i];
  }

  #define FVEC_MASK_BINFN_B(the_name, the_imm)                          \
    VEC_INLINE static BVEC_NAME the_name(const FVEC_NAME &a,            \
                                         const FVEC_NAME &b) {          \
      return FVEC_SUFFIX(_mm256_cmp_)(a.val_, b.val_, the_imm);         \
    }                                                                   \
    VEC_INLINE static BVEC_NAME mask_##the_name(                        \
        const BVEC_NAME &mask,                                          \
        const FVEC_NAME &a, const FVEC_NAME &b                          \
    ) {                                                                 \
      BVEC_NAME ret = FVEC_SUFFIX(_mm256_cmp_)(                         \
        a.val_, b.val_, the_imm);                                       \
      return mask & ret;                                                \
    }
  FVEC_MASK_BINFN_B(cmple, _CMP_LE_OS)
  FVEC_MASK_BINFN_B(cmplt, _CMP_LT_OS)
  FVEC_MASK_BINFN_B(cmpneq, _CMP_NEQ_UQ)
  FVEC_MASK_BINFN_B(cmpnle, _CMP_NLE_US)
  FVEC_MASK_BINFN_B(cmpnlt, _CMP_NLT_US)
  #undef FVEC_MASK_BINFN_B

  VEC_INLINE static __m256d _mm256_recip_pd(__m256d a) {
    __m256d c_1 = _mm256_set1_pd(1);
    return _mm256_div_pd(c_1, a);
  }
  VEC_INLINE static __m256 _mm256_recip_ps(__m256 a) {
    return _mm256_rcp_ps(a);
  }
  VEC_INLINE static __m256d _mm256_abs_pd(__m256d a) {
    const unsigned long long abs_mask = 0x7FFFFFFFFFFFFFFF;
    const unsigned long long abs_full[8] =
        {abs_mask, abs_mask, abs_mask, abs_mask, abs_mask, abs_mask, abs_mask,
           abs_mask};
    return _mm256_and_pd(_mm256_load_pd((double*)abs_full), a);
  }
  VEC_INLINE static __m256 _mm256_abs_ps(__m256 a) {
    const unsigned long long abs_mask = 0x7FFFFFFF;
    const unsigned long long abs_full[16] =
        {abs_mask, abs_mask, abs_mask, abs_mask, abs_mask, abs_mask, abs_mask,
           abs_mask, abs_mask, abs_mask, abs_mask, abs_mask, abs_mask,
           abs_mask, abs_mask, abs_mask};
    return _mm256_and_ps(_mm256_load_ps((float*)abs_full), a);
  }

  #define FVEC_UNFN_F(the_name)                                      \
    VEC_INLINE static FVEC_NAME the_name(const FVEC_NAME &a) {       \
      return FVEC_SUFFIX(_mm256_##the_name##_)(a.val_);              \
    }
  FVEC_UNFN_F(abs)
  FVEC_UNFN_F(exp)
  FVEC_UNFN_F(invsqrt)
  FVEC_UNFN_F(recip)
  FVEC_UNFN_F(sqrt)
  #undef FVEC_UNFN_F

  VEC_INLINE static FVEC_NAME mask_blend(
      const BVEC_NAME &mask, const FVEC_NAME &a, const FVEC_NAME &b
  ) {
    return FVEC_SUFFIX(_mm256_blendv_)(a.val_, b.val_, mask.val_);
  }
  #define FVEC_MASK_UNFN_F(the_name)                                 \
    VEC_INLINE static FVEC_NAME mask_##the_name(                     \
        const FVEC_NAME &src, const BVEC_NAME &mask,                 \
        const FVEC_NAME &a                                           \
    ) {                                                              \
      FVEC_NAME ret = FVEC_SUFFIX(_mm256_##the_name##_)(             \
                                                        a.val_);     \
      return mask_blend(mask, src, ret);                             \
    }
  FVEC_MASK_UNFN_F(cos)
  FVEC_MASK_UNFN_F(recip)
  FVEC_MASK_UNFN_F(sqrt)
  #undef FVEC_MASK_UNFN_F

  #define FVEC_BINFN_F(the_name)                                     \
    VEC_INLINE static FVEC_NAME the_name(const FVEC_NAME &a,         \
                                         const FVEC_NAME &b) {       \
      return FVEC_SUFFIX(_mm256_##the_name##_)(a.val_, b.val_);      \
    }
  FVEC_BINFN_F(max)
  FVEC_BINFN_F(min)
  #undef FVEC_BINFN_F

  #define FVEC_MASK_BINFN_F(the_name)                                \
    VEC_INLINE static FVEC_NAME mask_##the_name(                     \
        const FVEC_NAME &src, const BVEC_NAME &mask,                 \
        const FVEC_NAME &a, const FVEC_NAME &b                       \
    ) {                                                              \
      FVEC_NAME ret = FVEC_SUFFIX(_mm256_##the_name##_)(             \
        a.val_, b.val_);                                             \
      return mask_blend(mask, src, ret);                             \
    }
  FVEC_MASK_BINFN_F(add)
  FVEC_MASK_BINFN_F(div)
  FVEC_MASK_BINFN_F(mul)
  FVEC_MASK_BINFN_F(sub)
  #undef FVEC_MASK_BINFN_F

  VEC_INLINE static FVEC_NAME mask_expand(
      const FVEC_NAME &src, const BVEC_NAME &mask, const FVEC_NAME &b
  ) {
    FVEC_VEC_T ret = FVEC_SUFFIX(_mm256_expand_)(mask.val_, b.val_);
    ret = FVEC_SUFFIX(_mm256_and_)(mask.val_, ret);
    ret = FVEC_SUFFIX(_mm256_or_)(ret, FVEC_SUFFIX(_mm256_andnot_)
      (mask.val_, src.val_));
    return ret;
  }
  VEC_INLINE static FVEC_NAME masku_compress(
      const BVEC_NAME &mask, const FVEC_NAME &b
  ) {
    return FVEC_SUFFIX(_mm256_compress_)(mask.val_, b.val_);
  }

  VEC_INLINE static FVEC_NAME set1(const FVEC_SCAL_T &a) {
    return FVEC_SUFFIX(_mm256_set1_)(a);
  }
  VEC_INLINE static FVEC_NAME setzero() {
    return FVEC_SUFFIX(_mm256_setzero_)();
  }
  VEC_INLINE static FVEC_NAME undefined() {
    return FVEC_SUFFIX(_mm256_undefined_)();
  }

  VEC_INLINE static FVEC_NAME load(const FVEC_SCAL_T *mem) {
    return FVEC_SUFFIX(_mm256_load_)(mem);
  }
  VEC_INLINE static void store(FVEC_SCAL_T * dest, const FVEC_NAME &a) {
    FVEC_SUFFIX(_mm256_store_)(dest, a.val_);
  }


  VEC_INLINE static FVEC_NAME gather(const IVEC_NAME &idx,
    const FVEC_SCAL_T * mem, const int scale) {
    assert(scale == sizeof(FVEC_SCAL_T));
#   if FVEC_LEN==4
#    ifdef __AVX2__
    static const unsigned int mask_shuffle[8] __attribute__((aligned(32))) =
      {0, 2, 4, 6, 0, 0, 0, 0};
    __m256i m = _mm256_permutevar8x32_epi32(idx.val_,
      _mm256_load_si256((__m256i*)mask_shuffle));
    __m128i idx_short = _mm256_castsi256_si128(m);
    return FVEC_SUFFIX(_mm256_i32gather_)(mem, idx_short, sizeof(FVEC_SCAL_T));
#    else
    int idx_buf[8] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i*) idx_buf, idx.val_);
    double dest[4] __attribute__((aligned(32)));
    for (int i = 0; i < 4; i++) {
      dest[i] = mem[idx_buf[2*i]];
    }
    return _mm256_load_pd(dest);
#    endif
#   else
    return FVEC_SUFFIX(_mm256_i32gather_)(mem, idx.val_, sizeof(FVEC_SCAL_T));
#   endif
  }
  VEC_INLINE static FVEC_NAME mask_gather(
      const FVEC_NAME &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const FVEC_SCAL_T * mem, const int scale
  ) {
    assert(scale == sizeof(FVEC_SCAL_T));
#   if FVEC_LEN==4
#    ifdef __AVX2__
    static const unsigned int mask_shuffle[8] __attribute__((aligned(32))) =
      {0, 2, 4, 6, 0, 0, 0, 0};
    __m256i m = _mm256_permutevar8x32_epi32(idx.val_,
      _mm256_load_si256((__m256i*)mask_shuffle));
    __m128i idx_short = _mm256_castsi256_si128(m);
    return FVEC_SUFFIX(_mm256_mask_i32gather_)(src.val_, mem, idx_short,
      mask.val_, sizeof(FVEC_SCAL_T));
#    else
    int idx_buf[8] __attribute__((aligned(32)));
    int mask_buf[8] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i*) idx_buf, idx.val_);
    _mm256_store_pd((double*) mask_buf, mask.val_);
    double dest[4] __attribute__((aligned(32)));
    _mm256_store_pd((double*) dest, src.val_);
    for (int i = 0; i < 4; i++) {
      if (mask_buf[2*i])
        dest[i] = mem[idx_buf[2*i]];
    }
    return _mm256_load_pd(dest);
#    endif
#   else
    return FVEC_SUFFIX(_mm256_mask_i32gather_)(src.val_, mem, idx.val_,
      mask.val_, sizeof(FVEC_SCAL_T));
#   endif
  }

  VEC_INLINE static void gather_4_adjacent(const IVEC_NAME &idx,
      const FVEC_SCAL_T * mem, const int scale, FVEC_NAME * out_0,
      FVEC_NAME * out_1, FVEC_NAME * out_2, FVEC_NAME * out_3) {
    assert(scale == sizeof(FVEC_SCAL_T));
    int idx_buf[8] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i*) idx_buf, idx.val_);
#   if FVEC_LEN==4
    __m256d a0 = _mm256_load_pd(&mem[idx_buf[0]]);
    __m256d a1 = _mm256_load_pd(&mem[idx_buf[2]]);
    __m256d a2 = _mm256_load_pd(&mem[idx_buf[4]]);
    __m256d a3 = _mm256_load_pd(&mem[idx_buf[6]]);
    __m256d b0 = _mm256_unpacklo_pd(a0, a1);
    __m256d b1 = _mm256_unpackhi_pd(a0, a1);
    __m256d b2 = _mm256_unpacklo_pd(a2, a3);
    __m256d b3 = _mm256_unpackhi_pd(a2, a3);
    *out_0 = _mm256_permute2f128_pd(b0, b2, 0x20);
    *out_1 = _mm256_permute2f128_pd(b1, b3, 0x20);
    *out_2 = _mm256_permute2f128_pd(b0, b2, 0x31);
    *out_3 = _mm256_permute2f128_pd(b1, b3, 0x31);
#   else
    const float *e0 = &mem[idx_buf[0]];
    const float *e1 = &mem[idx_buf[1]];
    const float *e2 = &mem[idx_buf[2]];
    const float *e3 = &mem[idx_buf[3]];
    const float *e4 = &mem[idx_buf[4]];
    const float *e5 = &mem[idx_buf[5]];
    const float *e6 = &mem[idx_buf[6]];
    const float *e7 = &mem[idx_buf[7]];
    __m256 a0 = _mm256_loadu2_m128(e4, e0);
    __m256 a1 = _mm256_loadu2_m128(e5, e1);
    __m256 b0 = _mm256_unpacklo_ps(a0, a1);
    __m256 b1 = _mm256_unpackhi_ps(a0, a1);
    __m256 a2 = _mm256_loadu2_m128(e6, e2);
    __m256 a3 = _mm256_loadu2_m128(e7, e3);
    __m256 b2 = _mm256_unpacklo_ps(a2, a3);
    __m256 b3 = _mm256_unpackhi_ps(a2, a3);
    *out_0 = _mm256_shuffle_ps(b0, b2, 0x44);
    *out_1 = _mm256_shuffle_ps(b0, b2, 0xEE);
    *out_2 = _mm256_shuffle_ps(b1, b3, 0x44);
    *out_3 = _mm256_shuffle_ps(b1, b3, 0xEE);
#   endif
  }
  VEC_INLINE static void gather_3_adjacent(const IVEC_NAME &idx,
                                           const FVEC_SCAL_T * mem,
                                           const int scale,
                                           FVEC_NAME * out_0,
                                           FVEC_NAME * out_1,
                                           FVEC_NAME * out_2) {
    assert(scale == sizeof(FVEC_SCAL_T));
    FVEC_NAME tmp_3;
    gather_4_adjacent(idx, mem, scale, out_0, out_1, out_2, &tmp_3);
  }

  VEC_INLINE static double _mm256_reduce_add_pd(__m256d a) {
    __m256d t1 = _mm256_hadd_pd(a, a);
    __m128d t2 = _mm256_extractf128_pd(t1, 1);
    __m128d t3 = _mm256_castpd256_pd128(t1);
    return _mm_cvtsd_f64(_mm_add_pd(t2, t3));
  }

  VEC_INLINE static float _mm256_reduce_add_ps(__m256 a) {
    __m256 t1 = _mm256_hadd_ps(a, a);
    __m128 t2 = _mm256_extractf128_ps(t1, 1);
    __m128 t3 = _mm256_castps256_ps128(t1);
    __m128 t4 = _mm_add_ps(t2, t3);
    __m128 t5 = _mm_permute_ps(t4, 0x1B); // 0x1B = reverse
    return _mm_cvtss_f32(_mm_add_ps(t4, t5));
  }

  VEC_INLINE static FVEC_SCAL_T reduce_add(const FVEC_NAME &a) {
    return FVEC_SUFFIX(_mm256_reduce_add_)(a.val_);
  }
  VEC_INLINE static FVEC_SCAL_T mask_reduce_add(const BVEC_NAME &mask,
                                                const FVEC_NAME &a) {
    return reduce_add(FVEC_SUFFIX(_mm256_and_)(mask.val_, a.val_));
  }

  VEC_INLINE static IVEC_NAME unpackloepi32(const FVEC_NAME &a) {
#   if FVEC_LEN==4
#    if __AVX2__
    static const unsigned int mask_shuffle[8] __attribute__((aligned(32))) =
      {0, 0, 2, 2, 4, 4, 6, 6};
    __m256 m = _mm256_permutevar8x32_ps(_mm256_castpd_ps(a.val_),
      _mm256_load_si256((__m256i*)mask_shuffle));
    return _mm256_castps_si256(m);
#    else
    __m128i a_lo = _mm256_castsi256_si128(_mm256_castpd_si256(a.val_));
    __m128i a_hi = _mm256_extractf128_si256(_mm256_castpd_si256(a.val_), 1);
    __m128i c_lo = _mm_shuffle_epi32(a_lo, 0xA0); /*1010 0000*/
    __m128i c_hi = _mm_shuffle_epi32(a_hi, 0xA0);
    __m256i ret = _mm256_setr_m128i(c_lo, c_hi);
    return ret;
#    endif
#   else
    return _mm256_castps_si256(a.val_);
#   endif
  }

  VEC_INLINE static FVEC_NAME mask_sincos(
      FVEC_NAME * cos, const FVEC_NAME &src_a, const FVEC_NAME &src_b,
      const BVEC_NAME &mask, const FVEC_NAME &arg
  ) {
    FVEC_VEC_T c, s = FVEC_SUFFIX(_mm256_sincos_)(&c, arg.val_);
    *cos = mask_blend(mask, src_b, c);
    return mask_blend(mask, src_a, s);
  }

  #define FVEC_BINOP(the_sym, the_name)                                      \
    VEC_INLINE inline FVEC_NAME operator the_sym(const FVEC_NAME &b) const { \
    return FVEC_SUFFIX(_mm256_##the_name##_)(this->val_, b.val_);            \
  }
  FVEC_BINOP(+, add)
  FVEC_BINOP(-, sub)
  FVEC_BINOP(*, mul)
  FVEC_BINOP(/, div)
  #undef FVEC_BINOP

  VEC_INLINE static void gather_prefetch0(const IVEC_NAME &a, void * mem) {
    /* NOP */
  }
};

class AVEC_NAME {
  friend class avec8pd;
  FVEC_VEC_T val_;
  VEC_INLINE AVEC_NAME(const FVEC_VEC_T &a) : val_(a) {}
public:
  VEC_INLINE AVEC_NAME(const FVEC_NAME &a) : val_(a.val_) {}
  VEC_INLINE static AVEC_NAME undefined() {
    return FVEC_SUFFIX(_mm256_undefined_)();
  }
  VEC_INLINE static AVEC_NAME mask_gather(
      const AVEC_NAME &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const FVEC_SCAL_T * mem, const int scale
  ) {
    assert(scale == sizeof(FVEC_SCAL_T));
    return FVEC_NAME::mask_gather(src.val_, mask, idx, mem, scale);
  }
  VEC_INLINE static void mask_i32loscatter(
      FVEC_SCAL_T * mem, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const AVEC_NAME &a, const int scale
  ) {
    assert(scale == sizeof(FVEC_SCAL_T));
    for (int l = 0; l < FVEC_NAME::VL; l++) {
      if (BVEC_NAME::test_at(mask, l))
        mem[IVEC_NAME::at(idx, l)] = FVEC_NAME::at(a.val_, l);
    }
  }

  #define AVEC_BINOP(the_sym, the_name)                                      \
    VEC_INLINE inline AVEC_NAME operator the_sym(const AVEC_NAME &b) const { \
    return FVEC_SUFFIX(_mm256_##the_name##_)(this->val_, b.val_);            \
  }
  AVEC_BINOP(-, sub)
  #undef AVEC_BINOP
};

#if FVEC_LEN==8
class avec8pd {
  __m256d lo_, hi_;
  VEC_INLINE avec8pd(const __m256d &lo, const __m256d &hi) : lo_(lo), hi_(hi) {}
  VEC_INLINE static __m128 get_ps_hi(__m256 a) {
    return _mm256_extractf128_ps(a, 1);
  }
  VEC_INLINE static __m128 get_ps_lo(__m256 a) {
    return _mm256_castps256_ps128(a);
  }
  VEC_INLINE static __m128i get_si_hi(__m256i a) {
    return _mm_castps_si128(get_ps_hi(_mm256_castsi256_ps(a)));
  }
  VEC_INLINE static __m128i get_si_lo(__m256i a) {
    return _mm_castps_si128(get_ps_lo(_mm256_castsi256_ps(a)));
  }
public:
  VEC_INLINE avec8pd(const FVEC_NAME &a) {
    lo_ = _mm256_cvtps_pd(get_ps_lo(a.val_));
    hi_ = _mm256_cvtps_pd(get_ps_hi(a.val_));
  }
  VEC_INLINE static avec8pd undefined() {
    return avec8pd(_mm256_undefined_pd(), _mm256_undefined_pd());
  }
  VEC_INLINE static avec8pd mask_gather(
      const avec8pd &src, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const double * mem, const int scale
  ) {
#   ifndef __AVX2__
    assert(scale == sizeof(double));
    int idx_buf[8] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i*) idx_buf, idx.val_);
    int mask_val = _mm256_movemask_ps(mask.val_);
    double ret_buf[8] __attribute__((aligned(32)));
    _mm256_store_pd(&ret_buf[0], src.lo_);
    _mm256_store_pd(&ret_buf[4], src.hi_);
    for (int i = 0; i < 8; i++) {
      if (mask_val & (1 << i)) {
        ret_buf[i] = mem[idx_buf[i]];
      }
    }
    __m256d lo = _mm256_load_pd(&ret_buf[0]);
    __m256d hi = _mm256_load_pd(&ret_buf[4]);
#   else
    static const unsigned int lo_shuffle[8] __attribute__((aligned(32))) =
      {0, 0, 1, 1, 2, 2, 3, 3};
    static const unsigned int hi_shuffle[8] __attribute__((aligned(32))) =
      {4, 4, 5, 5, 6, 6, 7, 7};
    __m256d lo_mask = _mm256_castps_pd(_mm256_permutevar8x32_ps(mask.val_,
      _mm256_load_si256((__m256i*) lo_shuffle)));
    __m256d hi_mask = _mm256_castps_pd(_mm256_permutevar8x32_ps(mask.val_,
      _mm256_load_si256((__m256i*) hi_shuffle)));
    __m256d lo = _mm256_mask_i32gather_pd(src.lo_, mem, get_si_lo(idx.val_),
                                          lo_mask, sizeof(double));
    __m256d hi = _mm256_mask_i32gather_pd(src.hi_, mem, get_si_hi(idx.val_),
                                          hi_mask, sizeof(double));
#   endif
    return avec8pd(lo, hi);
  }
  VEC_INLINE static void mask_i32loscatter(
      double * mem, const BVEC_NAME &mask, const IVEC_NAME &idx,
      const avec8pd &a, const int scale
  ) {
    assert(scale == sizeof(double));
    double a_buf[8] __attribute__((aligned(32)));
    _mm256_store_pd(a_buf, a.lo_);
    _mm256_store_pd(&a_buf[4], a.hi_);
    int idx_buf[8] __attribute__((aligned(32)));
    _mm256_store_si256((__m256i*)idx_buf, idx.val_);
    int mask_val = _mm256_movemask_ps(mask.val_);
    for (int i = 0; i < 8; i++) {
      if (mask_val & (1 << i))
        mem[idx_buf[i]] = a_buf[i];
    }
  }

  #define AVEC2_BINOP(the_sym, the_name)                                    \
    VEC_INLINE inline avec8pd operator the_sym(const avec8pd &b) const {    \
    __m256d lo = _mm256_##the_name##_pd(this->lo_, b.lo_);                  \
    __m256d hi = _mm256_##the_name##_pd(this->hi_, b.hi_);                  \
    return avec8pd(lo, hi);                                                 \
  }
  AVEC2_BINOP(-, sub)
};
#endif

}


#ifdef FVEC_FIRST_PASS

template<typename flt_t, typename acc_t>
struct intr_types;

template<>
struct intr_types<double,double> {
  typedef mm256::fvec4pd fvec;
  typedef mm256::ivec4 ivec;
  typedef mm256::bvec4 bvec;
  typedef mm256::avec4pd avec;
};

template<>
struct intr_types<float,float> {
  typedef mm256::fvec8ps fvec;
  typedef mm256::ivec8 ivec;
  typedef mm256::bvec8 bvec;
  typedef mm256::avec8ps avec;
};

template<>
struct intr_types<float,double> {
  typedef mm256::fvec8ps fvec;
  typedef mm256::ivec8 ivec;
  typedef mm256::bvec8 bvec;
  typedef mm256::avec8pd avec;
};

#endif

#ifndef FVEC_FIRST_PASS
#  define FVEC_FIRST_PASS
#  include "intel_intrinsics_airebo.h"
#endif

#endif

#ifdef LMP_INTEL_AIREBO_SCALAR

#include <cassert>
#include <cmath>
#include <immintrin.h>

#define VEC_INLINE __attribute__((always_inline))

template<typename flt_t, typename acc_t>
struct intr_types {

class fvec;
class ivec;
class avec;
class bvec {
  friend class fvec;
  friend class ivec;
  friend class avec;
  bool val_;
  VEC_INLINE bvec(const bool &v) : val_(v) {}
public:
  VEC_INLINE bvec() {}
  VEC_INLINE static bvec kand(const bvec &a, const bvec &b) {
    return a.val_ && b.val_;
  }
  VEC_INLINE static bvec kandn(const bvec &a, const bvec &b) {
    return (! a.val_) && b.val_;
  }
  VEC_INLINE static bvec knot(const bvec &a) {
    return ! a.val_;
  }
  VEC_INLINE static int kortestz(const bvec &a, const bvec &b) {
    return (! a.val_) && (! b.val_) ? true : false;
  }
  VEC_INLINE static bvec masku_compress(const bvec &mask, const bvec &a) {
    return mask.val_ ? a.val_ : false;
  }
  VEC_INLINE static bvec mask_expand(const bvec &src, const bvec &mask,
                                     const bvec &a) {
    return mask.val_ ? a.val_ : src.val_;
  }
  VEC_INLINE static bvec full() {
    return true;
  }
  VEC_INLINE static bvec empty() {
    return false;
  }
  VEC_INLINE static bvec only(int n) {
    return n == 1 ? true : false;
  }
  VEC_INLINE static bvec after(int n) {
    return n == 0 ? true : false;
  }
  VEC_INLINE static bvec onlyafter(int only, int after) {
    return after == 0 && only == 1 ? true : false;
  }
  VEC_INLINE static int popcnt(const bvec &a) {
    return static_cast<int>(a.val_);
  }
  VEC_INLINE static bool test_all_unset(const bvec &a) {
    return kortestz(a, a);
  }
  VEC_INLINE static bool test_any_set(const bvec &a) {
    return ! test_all_unset(a);
  }
  VEC_INLINE static bool test_at(const bvec &a, int i) {
    assert(i < 1);
    return a.val_;
  }
  VEC_INLINE bvec operator &(const bvec &b) const {
    return val_ && b.val_;
  }
  VEC_INLINE bvec operator |(const bvec &b) const {
    return val_ || b.val_;
  }
  VEC_INLINE bvec operator ~() const {
    return ! val_;
  }
};

class ivec {
  friend class fvec;
  friend class avec;
  int val_;
  VEC_INLINE ivec(const int &v) : val_(v) {}
public:
  static const int VL = 1;
  VEC_INLINE ivec() {}

  #define IVEC_MASK_BINFN_B(the_name, the_op)                        \
    VEC_INLINE static bvec the_name(const ivec &a, const ivec &b) {  \
      return a.val_ the_op b.val_;                                   \
    }                                                                \
    VEC_INLINE static bvec mask_##the_name(                          \
        const bvec &mask,                                            \
        const ivec &a, const ivec &b                                 \
    ) {                                                              \
      return mask.val_ && (a.val_ the_op b.val_);                    \
                                                                     \
    }
  IVEC_MASK_BINFN_B(cmpeq, ==)
  IVEC_MASK_BINFN_B(cmplt, <)
  IVEC_MASK_BINFN_B(cmpneq, !=)
  IVEC_MASK_BINFN_B(cmpgt, >)

  #define IVEC_MASK_BINFN_I(the_name, the_op)                        \
    VEC_INLINE static ivec mask_##the_name(                          \
        const ivec &src, const bvec &mask,                           \
        const ivec &a, const ivec &b                                 \
    ) {                                                              \
      return mask.val_ ? a.val_ the_op b.val_ : src.val_;            \
    }
  IVEC_MASK_BINFN_I(add, +)
  VEC_INLINE static ivec mask_blend(
      const bvec &mask, const ivec &a, const ivec &b
  ) {
    return mask.val_ ? b.val_ : a.val_;
  }

  #define IVEC_BINFN_I(the_name, the_op)                             \
    VEC_INLINE static ivec the_name(const ivec &a, const ivec &b) {  \
      return a.val_ the_op b.val_;                                   \
    }
  IVEC_BINFN_I(mullo, *)
  IVEC_BINFN_I(srlv, >>)
  VEC_INLINE static ivec the_and(const ivec &a, const ivec &b) {
    return a.val_ & b.val_;
  }

  VEC_INLINE static ivec mask_expand(
      const ivec &src, const bvec &a, const ivec &b
  ) {
    return a.val_ ? b.val_ : src.val_;
  }
  VEC_INLINE static ivec masku_compress(
      const bvec &a, const ivec &b
  ) {
    return a.val_ ? b.val_ : 0;
  }

  VEC_INLINE static int at(const ivec &a, int b) {
    assert(b == 0);
    return a.val_;
  }

  VEC_INLINE static ivec load(const int * src) {
    return *src;
  }
  VEC_INLINE static ivec mask_loadu(const bvec &mask, const int * src) {
    return mask.val_ ? *src : 0xDEAD;
  }
  VEC_INLINE static ivec maskz_loadu(const bvec &mask, const int * src) {
    return mask.val_ ? *src : 0;
  }
  VEC_INLINE static void mask_storeu(const bvec &mask, int * dest,
    const ivec &src) {
    if (mask.val_) *dest = src.val_;
  }
  VEC_INLINE static void store(int * dest, const ivec &src) {
    *dest = src.val_;
  }

  VEC_INLINE static ivec mask_gather(
      const ivec &src, const bvec &mask, const ivec &idx, const int * mem,
        const int scale
  ) {
    return mask.val_ ? *reinterpret_cast<const int *>
      (reinterpret_cast<const char*>(mem) + scale * idx.val_) : src.val_;
  }
  VEC_INLINE static void mask_i32scatter(
      int * mem, const bvec &mask, const ivec &idx, const ivec &a,
        const int scale
  ) {
    if (mask.val_) *reinterpret_cast<int *>(reinterpret_cast<char*>(mem) +
      scale * idx.val_) = a.val_;
  }

  VEC_INLINE static void mask_compressstore(const bvec &mask, int * dest,
      const ivec &src) {
    if (mask.val_) *dest = src.val_;
  }

  VEC_INLINE static ivec set(
      int /*i15*/, int /*i14*/, int /*i13*/, int /*i12*/, int /*i11*/, int /*i10*/, int /*i9*/, int /*i8*/,
      int /*i7*/, int /*i6*/, int /*i5*/, int /*i4*/, int /*i3*/, int /*i2*/, int /*i1*/, int i0
  ) {
    return i0;
  }
  VEC_INLINE static ivec set1(int i) {
    return i;
  }
  VEC_INLINE static ivec setzero() {
    return 0;
  }
  VEC_INLINE static ivec undefined() {
    return 0xDEAD;
  }

  VEC_INLINE ivec operator +(const ivec &b) const {
    return val_ + b.val_;
  }
};

class fvec {
  friend class avec;
  flt_t val_;
  VEC_INLINE fvec(const flt_t &v) : val_(v) {}
public:
  static const int VL = 1;
  VEC_INLINE fvec() {}
  VEC_INLINE static flt_t at(const fvec &a, int i) {
    assert(i < 1);
    return a.val_;
  }
  VEC_INLINE static bool fast_compress() { return false; }

  #define FVEC_MASK_BINFN_B(the_name, the_op)                        \
    VEC_INLINE static bvec the_name(const fvec &a, const fvec &b) {  \
      return a.val_ the_op b.val_;                                   \
    }                                                                \
    VEC_INLINE static bvec mask_##the_name(                          \
        const bvec &mask,                                            \
        const fvec &a, const fvec &b                                 \
    ) {                                                              \
      return mask.val_ && (a.val_ the_op b.val_);                    \
    }
  FVEC_MASK_BINFN_B(cmple, <=)
  FVEC_MASK_BINFN_B(cmplt, <)
  FVEC_MASK_BINFN_B(cmpneq, !=)
  FVEC_MASK_BINFN_B(cmpnle, >)
  FVEC_MASK_BINFN_B(cmpnlt, >=)

  #define FVEC_UNFN_F(the_name, the_fn)                              \
    VEC_INLINE static fvec the_name(const fvec &a) {                 \
      return the_fn(a.val_);                                         \
    }
  FVEC_UNFN_F(abs, fabs)
  FVEC_UNFN_F(exp, ::exp)
  FVEC_UNFN_F(invsqrt, 1/std::sqrt)
  FVEC_UNFN_F(recip, 1/)
  FVEC_UNFN_F(sqrt, std::sqrt)

  #define FVEC_MASK_UNFN_F(the_name, the_fn)                         \
    VEC_INLINE static fvec mask_##the_name(                          \
        const fvec &src, const bvec &mask,                           \
        const fvec &a                                                \
    ) {                                                              \
      return mask.val_ ? the_fn(a.val_) : src.val_;                  \
    }
  FVEC_MASK_UNFN_F(cos, std::cos)
  FVEC_MASK_UNFN_F(recip, 1/)
  FVEC_MASK_UNFN_F(sqrt, std::sqrt)

  #define FVEC_BINFN_F(the_name, the_fn)                             \
    VEC_INLINE static fvec the_name(const fvec &a, const fvec &b) {  \
      return the_fn(a.val_, b.val_);                                 \
    }
  FVEC_BINFN_F(max, ::fmax)
  FVEC_BINFN_F(min, ::fmin)

  #define FVEC_MASK_BINFN_F(the_name, the_op)                        \
    VEC_INLINE static fvec mask_##the_name(                          \
        const fvec &src, const bvec &mask,                           \
        const fvec &a, const fvec &b                                 \
    ) {                                                              \
      return mask.val_ ? a.val_ the_op b.val_ : src.val_;            \
    }
  FVEC_MASK_BINFN_F(add, +)
  FVEC_MASK_BINFN_F(div, /)
  FVEC_MASK_BINFN_F(mul, *)
  FVEC_MASK_BINFN_F(sub, -)
  VEC_INLINE static fvec mask_blend(
      const bvec &mask, const fvec &a, const fvec &b
  ) {
    return mask.val_ ? b.val_ : a.val_;
  }

  VEC_INLINE static fvec mask_expand(
      const fvec &src, const bvec &a, const fvec &b
  ) {
    return a.val_ ? b.val_ : src.val_;
  }
  VEC_INLINE static fvec masku_compress(
      const bvec &a, const fvec &b
  ) {
    return a.val_ ? b.val_ : 0;
  }

  VEC_INLINE static fvec set1(const flt_t &a) {
    return a;
  }
  VEC_INLINE static fvec setzero() {
    return 0;
  }
  VEC_INLINE static fvec undefined() {
    return 1337.1337;
  }

  VEC_INLINE static fvec load(const flt_t *mem) {
    return *mem;
  }
  VEC_INLINE static void mask_storeu(const bvec &mask, flt_t * dest,
                                     const fvec &a) {
    if (mask.val_) *dest = a.val_;
  }
  VEC_INLINE static void store(flt_t * dest, const fvec &a) {
    *dest = a.val_;
  }

  VEC_INLINE static fvec gather(const ivec &idx, const flt_t * mem,
                                const int scale) {
    return *reinterpret_cast<const flt_t*>(reinterpret_cast<const char*>(mem) +
      scale * idx.val_);
  }
  VEC_INLINE static fvec mask_gather(
      const fvec &src, const bvec &mask, const ivec &idx,
      const flt_t * mem, const int scale
  ) {
    return mask.val_ ? *reinterpret_cast<const flt_t*>
      (reinterpret_cast<const char*>(mem) + scale * idx.val_) : src.val_;
  }

  VEC_INLINE static void gather_3_adjacent(const ivec &idx, const flt_t * mem,
                                           const int scale, fvec * out_0,
                                           fvec * out_1, fvec * out_2) {
    assert(scale == sizeof(flt_t));
    *out_0 = gather(idx, mem + 0, scale);
    *out_1 = gather(idx, mem + 1, scale);
    *out_2 = gather(idx, mem + 2, scale);
  }
  VEC_INLINE static void gather_4_adjacent(const ivec &idx, const flt_t * mem,
                                           const int scale, fvec * out_0,
                                           fvec * out_1, fvec * out_2,
                                           fvec * out_3) {
    assert(scale == sizeof(flt_t));
    *out_0 = gather(idx, mem + 0, scale);
    *out_1 = gather(idx, mem + 1, scale);
    *out_2 = gather(idx, mem + 2, scale);
    *out_3 = gather(idx, mem + 3, scale);
  }

  VEC_INLINE static flt_t mask_reduce_add(const bvec &mask, const fvec &a) {
    return mask.val_ ? a.val_ : 0;
  }
  VEC_INLINE static flt_t reduce_add(const fvec &a) {
    return a.val_;
  }

  VEC_INLINE static ivec unpackloepi32(const fvec &a) {
    union { int i; flt_t f; } atype;
    atype.f = a.val_;
    return ivec(atype.i);
  }

  VEC_INLINE static fvec mask_sincos(
      fvec * cos_out, const fvec &src_a, const fvec &src_b,
      const bvec &mask, const fvec &arg
  ) {
    cos_out->val_ = mask.val_ ? ::cos(arg.val_) : src_b.val_;
    return mask.val_ ? ::sin(arg.val_) : src_a.val_;
  }

  #define FVEC_BINOP(the_sym, the_name)                              \
    VEC_INLINE inline fvec operator the_sym(const fvec &b) const {   \
    return this->val_ the_sym b.val_;                                \
  }
  FVEC_BINOP(+, add)
  FVEC_BINOP(-, sub)
  FVEC_BINOP(*, mul)
  FVEC_BINOP(/, div)

  VEC_INLINE static void gather_prefetch0(const ivec & /*idx*/, const void * /*mem*/) {}
};

class avec {
  acc_t val_;
  VEC_INLINE avec(const acc_t &a) : val_(a) {}
public:
  VEC_INLINE avec(const fvec &a) : val_(a.val_) {}
  VEC_INLINE static avec undefined() {
    return 1337.1337;
  }
  VEC_INLINE static avec mask_gather(const avec &src, const bvec &mask,
                                     const ivec &idx, const acc_t * mem,
                                     const int scale) {
    return mask.val_ ? *reinterpret_cast<const acc_t*>
      (reinterpret_cast<const char*>(mem) + scale * idx.val_) : src.val_;
  }
  VEC_INLINE static void mask_i32loscatter(acc_t * mem, const bvec &mask,
                                           const ivec &idx, const avec &a,
                                           const int scale) {
    if (mask.val_) *reinterpret_cast<acc_t*>(reinterpret_cast<char*>(mem) +
                                             idx.val_ * scale) = a.val_;
  }

  #define AVEC_BINOP(the_sym, the_name)                              \
    VEC_INLINE inline avec operator the_sym(const avec &b) const {   \
    return this->val_ the_sym b.val_;                                \
  }
  AVEC_BINOP(-, sub)
};

};

#endif
