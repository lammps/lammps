/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
   we use a stripped down KISS FFT as default FFT for LAMMPS
   this code is adapted from kiss_fft_v1_2_9
   homepage: http://kissfft.sf.net/

   changes 2008-2011 by Axel Kohlmeyer <akohlmey@gmail.com>
*/

#ifndef LMP_FFT_KISSFFT
#define LMP_FFT_KISSFFT

#include <cmath>
#include <cstdlib>
#include <cstring>

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944
#endif

/*
  Copyright (c) 2003-2010, Mark Borgerding

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    * Neither the author nor the names of any contributors may be used
      to endorse or promote products derived from this software without
      specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#define KISS_FFT_MALLOC malloc
#define KISS_FFT_FREE free
#define MAXFACTORS 32
/* e.g. an fft of length 128 has 4 factors
 as far as kissfft is concerned: 4*4*4*2  */
struct kiss_fft_state {
  int nfft;
  int inverse;
  int factors[2 * MAXFACTORS];
  FFT_DATA twiddles[1];
};

#ifdef KISS_FFT_USE_ALLOCA
// define this to allow use of alloca instead of malloc for temporary buffers
// Temporary buffers are used in two case:
// 1. FFT sizes that have "bad" factors. i.e. not 2,3 and 5
// 2. "in-place" FFTs.  Notice the quotes, since kissfft does not really do an in-place transform.
#include <alloca.h>
#define KISS_FFT_TMP_ALLOC(nbytes) alloca(nbytes)
#define KISS_FFT_TMP_FREE(ptr)
#else
#define KISS_FFT_TMP_ALLOC(nbytes) KISS_FFT_MALLOC(nbytes)
#define KISS_FFT_TMP_FREE(ptr) KISS_FFT_FREE(ptr)
#endif

static kiss_fft_cfg kiss_fft_alloc(int, int, void *, size_t *);
static void kiss_fft(kiss_fft_cfg, const FFT_DATA *, FFT_DATA *);

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */

#define S_MUL(a, b) ((a) * (b))
#define C_MUL(m, a, b)                          \
  do {                                          \
    (m).re = (a).re * (b).re - (a).im * (b).im; \
    (m).im = (a).re * (b).im + (a).im * (b).re; \
  } while (0)
#define C_FIXDIV(c, div) /* NOOP */
#define C_MULBYSCALAR(c, s) \
  do {                      \
    (c).re *= (s);          \
    (c).im *= (s);          \
  } while (0)

#ifndef CHECK_OVERFLOW_OP
#define CHECK_OVERFLOW_OP(a, op, b) /* noop */
#endif

#define C_ADD(res, a, b)                 \
  do {                                   \
    CHECK_OVERFLOW_OP((a).re, +, (b).re) \
    CHECK_OVERFLOW_OP((a).im, +, (b).im) \
    (res).re = (a).re + (b).re;          \
    (res).im = (a).im + (b).im;          \
  } while (0)
#define C_SUB(res, a, b)                 \
  do {                                   \
    CHECK_OVERFLOW_OP((a).re, -, (b).re) \
    CHECK_OVERFLOW_OP((a).im, -, (b).im) \
    (res).re = (a).re - (b).re;          \
    (res).im = (a).im - (b).im;          \
  } while (0)
#define C_ADDTO(res, a)                    \
  do {                                     \
    CHECK_OVERFLOW_OP((res).re, +, (a).re) \
    CHECK_OVERFLOW_OP((res).im, +, (a).im) \
    (res).re += (a).re;                    \
    (res).im += (a).im;                    \
  } while (0)

#define C_SUBFROM(res, a)                  \
  do {                                     \
    CHECK_OVERFLOW_OP((res).re, -, (a).re) \
    CHECK_OVERFLOW_OP((res).im, -, (a).im) \
    (res).re -= (a).re;                    \
    (res).im -= (a).im;                    \
  } while (0)

#define KISS_FFT_COS(phase) (kiss_fft_scalar) cos(phase)
#define KISS_FFT_SIN(phase) (kiss_fft_scalar) sin(phase)
#define HALF_OF(x) ((x) *.5)

#define kf_cexp(x, phase)          \
  do {                             \
    (x)->re = KISS_FFT_COS(phase); \
    (x)->im = KISS_FFT_SIN(phase); \
  } while (0)

static void kf_bfly2(FFT_DATA *Fout, const size_t fstride, const kiss_fft_cfg st, int m)
{
  FFT_DATA *Fout2;
  FFT_DATA *tw1 = st->twiddles;
  FFT_DATA t;

  Fout2 = Fout + m;
  do {
    C_FIXDIV(*Fout, 2);
    C_FIXDIV(*Fout2, 2);

    C_MUL(t, *Fout2, *tw1);
    tw1 += fstride;
    C_SUB(*Fout2, *Fout, t);
    C_ADDTO(*Fout, t);
    ++Fout2;
    ++Fout;
  } while (--m);
}

static void kf_bfly4(FFT_DATA *Fout, const size_t fstride, const kiss_fft_cfg st, const size_t m)
{
  FFT_DATA *tw1, *tw2, *tw3;
  FFT_DATA scratch[6];
  size_t k = m;
  const size_t m2 = 2 * m;
  const size_t m3 = 3 * m;

  tw3 = tw2 = tw1 = st->twiddles;

  do {
    C_FIXDIV(*Fout, 4);
    C_FIXDIV(Fout[m], 4);
    C_FIXDIV(Fout[m2], 4);
    C_FIXDIV(Fout[m3], 4);

    C_MUL(scratch[0], Fout[m], *tw1);
    C_MUL(scratch[1], Fout[m2], *tw2);
    C_MUL(scratch[2], Fout[m3], *tw3);

    C_SUB(scratch[5], *Fout, scratch[1]);
    C_ADDTO(*Fout, scratch[1]);
    C_ADD(scratch[3], scratch[0], scratch[2]);
    C_SUB(scratch[4], scratch[0], scratch[2]);
    C_SUB(Fout[m2], *Fout, scratch[3]);
    tw1 += fstride;
    tw2 += fstride * 2;
    tw3 += fstride * 3;
    C_ADDTO(*Fout, scratch[3]);

    if (st->inverse) {
      Fout[m].re = scratch[5].re - scratch[4].im;
      Fout[m].im = scratch[5].im + scratch[4].re;
      Fout[m3].re = scratch[5].re + scratch[4].im;
      Fout[m3].im = scratch[5].im - scratch[4].re;
    } else {
      Fout[m].re = scratch[5].re + scratch[4].im;
      Fout[m].im = scratch[5].im - scratch[4].re;
      Fout[m3].re = scratch[5].re - scratch[4].im;
      Fout[m3].im = scratch[5].im + scratch[4].re;
    }
    ++Fout;
  } while (--k);
}

static void kf_bfly3(FFT_DATA *Fout, const size_t fstride, const kiss_fft_cfg st, size_t m)
{
  size_t k = m;
  const size_t m2 = 2 * m;
  FFT_DATA *tw1, *tw2;
  FFT_DATA scratch[5];
  FFT_DATA epi3;
  epi3 = st->twiddles[fstride * m];

  tw1 = tw2 = st->twiddles;

  do {
    C_FIXDIV(*Fout, 3);
    C_FIXDIV(Fout[m], 3);
    C_FIXDIV(Fout[m2], 3);

    C_MUL(scratch[1], Fout[m], *tw1);
    C_MUL(scratch[2], Fout[m2], *tw2);

    C_ADD(scratch[3], scratch[1], scratch[2]);
    C_SUB(scratch[0], scratch[1], scratch[2]);
    tw1 += fstride;
    tw2 += fstride * 2;

    Fout[m].re = Fout->re - HALF_OF(scratch[3].re);
    Fout[m].im = Fout->im - HALF_OF(scratch[3].im);

    C_MULBYSCALAR(scratch[0], epi3.im);

    C_ADDTO(*Fout, scratch[3]);

    Fout[m2].re = Fout[m].re + scratch[0].im;
    Fout[m2].im = Fout[m].im - scratch[0].re;

    Fout[m].re -= scratch[0].im;
    Fout[m].im += scratch[0].re;

    ++Fout;
  } while (--k);
}

static void kf_bfly5(FFT_DATA *Fout, const size_t fstride, const kiss_fft_cfg st, int m)
{
  FFT_DATA *Fout0, *Fout1, *Fout2, *Fout3, *Fout4;
  int u;
  FFT_DATA scratch[13];
  FFT_DATA *twiddles = st->twiddles;
  FFT_DATA *tw;
  FFT_DATA ya, yb;
  ya = twiddles[fstride * m];
  yb = twiddles[fstride * 2 * m];

  Fout0 = Fout;
  Fout1 = Fout0 + m;
  Fout2 = Fout0 + 2 * m;
  Fout3 = Fout0 + 3 * m;
  Fout4 = Fout0 + 4 * m;

  tw = st->twiddles;
  for (u = 0; u < m; ++u) {
    C_FIXDIV(*Fout0, 5);
    C_FIXDIV(*Fout1, 5);
    C_FIXDIV(*Fout2, 5);
    C_FIXDIV(*Fout3, 5);
    C_FIXDIV(*Fout4, 5);
    scratch[0] = *Fout0;

    C_MUL(scratch[1], *Fout1, tw[u * fstride]);
    C_MUL(scratch[2], *Fout2, tw[2 * u * fstride]);
    C_MUL(scratch[3], *Fout3, tw[3 * u * fstride]);
    C_MUL(scratch[4], *Fout4, tw[4 * u * fstride]);

    C_ADD(scratch[7], scratch[1], scratch[4]);
    C_SUB(scratch[10], scratch[1], scratch[4]);
    C_ADD(scratch[8], scratch[2], scratch[3]);
    C_SUB(scratch[9], scratch[2], scratch[3]);

    Fout0->re += scratch[7].re + scratch[8].re;
    Fout0->im += scratch[7].im + scratch[8].im;

    scratch[5].re = scratch[0].re + S_MUL(scratch[7].re, ya.re) + S_MUL(scratch[8].re, yb.re);
    scratch[5].im = scratch[0].im + S_MUL(scratch[7].im, ya.re) + S_MUL(scratch[8].im, yb.re);

    scratch[6].re = S_MUL(scratch[10].im, ya.im) + S_MUL(scratch[9].im, yb.im);
    scratch[6].im = -S_MUL(scratch[10].re, ya.im) - S_MUL(scratch[9].re, yb.im);

    C_SUB(*Fout1, scratch[5], scratch[6]);
    C_ADD(*Fout4, scratch[5], scratch[6]);

    scratch[11].re = scratch[0].re + S_MUL(scratch[7].re, yb.re) + S_MUL(scratch[8].re, ya.re);
    scratch[11].im = scratch[0].im + S_MUL(scratch[7].im, yb.re) + S_MUL(scratch[8].im, ya.re);
    scratch[12].re = -S_MUL(scratch[10].im, yb.im) + S_MUL(scratch[9].im, ya.im);
    scratch[12].im = S_MUL(scratch[10].re, yb.im) - S_MUL(scratch[9].re, ya.im);

    C_ADD(*Fout2, scratch[11], scratch[12]);
    C_SUB(*Fout3, scratch[11], scratch[12]);

    ++Fout0;
    ++Fout1;
    ++Fout2;
    ++Fout3;
    ++Fout4;
  }
}

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(FFT_DATA *Fout, const size_t fstride, const kiss_fft_cfg st, int m,
                            int p)
{
  int u, k, q1, q;
  FFT_DATA *twiddles = st->twiddles;
  FFT_DATA t;
  int Norig = st->nfft;

  FFT_DATA *scratch = (FFT_DATA *) KISS_FFT_TMP_ALLOC(sizeof(FFT_DATA) * p);
  for (u = 0; u < m; ++u) {
    k = u;
    for (q1 = 0; q1 < p; ++q1) {
      scratch[q1] = Fout[k];
      C_FIXDIV(scratch[q1], p);
      k += m;
    }

    k = u;
    for (q1 = 0; q1 < p; ++q1) {
      int twidx = 0;
      Fout[k] = scratch[0];
      for (q = 1; q < p; ++q) {
        twidx += fstride * k;
        if (twidx >= Norig) twidx -= Norig;
        C_MUL(t, scratch[q], twiddles[twidx]);
        C_ADDTO(Fout[k], t);
      }
      k += m;
    }
  }
  KISS_FFT_TMP_FREE(scratch);
}

static void kf_work(FFT_DATA *Fout, const FFT_DATA *f, const size_t fstride, int in_stride,
                    int *factors, const kiss_fft_cfg st)
{
  FFT_DATA *Fout_beg = Fout;
  const int p = *factors++; /* the radix  */
  const int m = *factors++; /* stage's fft length/p */
  const FFT_DATA *Fout_end = Fout + p * m;

  if (m == 1) {
    do {
      *Fout = *f;
      f += fstride * in_stride;
    } while (++Fout != Fout_end);
  } else {
    do {
      /* recursive call:
               DFT of size m*p performed by doing
               p instances of smaller DFTs of size m,
               each one takes a decimated version of the input */
      kf_work(Fout, f, fstride * p, in_stride, factors, st);
      f += fstride * in_stride;
    } while ((Fout += m) != Fout_end);
  }

  Fout = Fout_beg;

  /* recombine the p smaller DFTs */
  switch (p) {
    case 2:
      kf_bfly2(Fout, fstride, st, m);
      break;
    case 3:
      kf_bfly3(Fout, fstride, st, m);
      break;
    case 4:
      kf_bfly4(Fout, fstride, st, m);
      break;
    case 5:
      kf_bfly5(Fout, fstride, st, m);
      break;
    default:
      kf_bfly_generic(Fout, fstride, st, m, p);
      break;
  }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static void kf_factor(int n, int *facbuf)
{
  int p = 4, nf = 0;
  double floor_sqrt;
  floor_sqrt = floor(sqrt((double) n));

  /* factor out the remaining powers of 4, powers of 2,
       and then any other remaining primes */
  do {
    if (nf == MAXFACTORS) p = n; /* make certain that we don't run out of space */
    while (n % p) {
      switch (p) {
        case 4:
          p = 2;
          break;
        case 2:
          p = 3;
          break;
        default:
          p += 2;
          break;
      }
      if (p > floor_sqrt) p = n; /* no more factors, skip to end */
    }
    n /= p;
    *facbuf++ = p;
    *facbuf++ = n;
    ++nf;
  } while (n > 1);
}

/*
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 */
static kiss_fft_cfg kiss_fft_alloc(int nfft, int inverse_fft, void *mem, size_t *lenmem)
{
  kiss_fft_cfg st = nullptr;
  size_t memneeded =
      sizeof(struct kiss_fft_state) + sizeof(FFT_DATA) * (nfft - 1); /* twiddle factors */

  if (lenmem == nullptr) {
    st = (kiss_fft_cfg) KISS_FFT_MALLOC(memneeded);
  } else {
    if (mem != nullptr && *lenmem >= memneeded) st = (kiss_fft_cfg) mem;
    *lenmem = memneeded;
  }

  if (st) {
    int i;
    st->nfft = nfft;
    st->inverse = inverse_fft;

    for (i = 0; i < nfft; ++i) {
      const double phase = (st->inverse ? 2.0 * M_PI : -2.0 * M_PI) * i / nfft;
      kf_cexp(st->twiddles + i, phase);
    }

    kf_factor(nfft, st->factors);
  }
  return st;
}

static void kiss_fft_stride(kiss_fft_cfg st, const FFT_DATA *fin, FFT_DATA *fout, int in_stride)
{
  if (fin == fout) {
    // NOTE: this is not really an in-place FFT algorithm.
    // It just performs an out-of-place FFT into a temp buffer
    FFT_DATA *tmpbuf = (FFT_DATA *) KISS_FFT_TMP_ALLOC(sizeof(FFT_DATA) * st->nfft);
    kf_work(tmpbuf, fin, 1, in_stride, st->factors, st);
    memcpy(fout, tmpbuf, sizeof(FFT_DATA) * st->nfft);
    KISS_FFT_TMP_FREE(tmpbuf);
  } else {
    kf_work(fout, fin, 1, in_stride, st->factors, st);
  }
}

static void kiss_fft(kiss_fft_cfg cfg, const FFT_DATA *fin, FFT_DATA *fout)
{
  kiss_fft_stride(cfg, fin, fout, 1);
}

#endif
