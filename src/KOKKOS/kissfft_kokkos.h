/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
   we use a stripped down KISS FFT as default FFT for LAMMPS
   this code is adapted from kiss_fft_v1_2_9
   homepage: http://kissfft.sf.net/

   changes 2008-2011 by Axel Kohlmeyer <akohlmey@gmail.com>
*/

#ifndef LMP_KISSFFT_KOKKOS_H
#define LMP_KISSFFT_KOKKOS_H

#include "fftdata_kokkos.h"
#include "kokkos_type.h"

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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kissfft_kokkos.h"
#include "fftdata_kokkos.h"

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

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
   C_EQ( res , a)       : res = a
 * */

#define S_MUL(a,b) ( (a)*(b) )

#define C_MUL(m,a,a_index,b,b_index) \
    do{ (m)[0] = (a)(a_index,0)*(b)(b_index,0) - (a)(a_index,1)*(b)(b_index,1);\
        (m)[1] = (a)(a_index,0)*(b)(b_index,1) + (a)(a_index,1)*(b)(b_index,0); }while(0)

#define C_FIXDIV(c,div) /* NOOP */



#define C_MULBYSCALAR( c, s ) \
    do{ (c)[0] *= (s);\
        (c)[1] *= (s); }while(0)

//#define  C_ADD( res, a,b)\
//    do { \
//            (res)[0]=(a)[0]+(b)[0];  (res)[1]=(a)[1]+(b)[1]; \
//    }while(0)

#define  C_SUB( res, a,b)\
    do { \
            (res)[0]=(a)[0]-(b)[0];  (res)[1]=(a)[1]-(b)[1]; \
    }while(0)

#define C_ADDTO( res , a)\
    do { \
            (res)[0] += (a)[0];  (res)[1] += (a)[1];\
    }while(0)

#define C_SUBFROM( res , a)\
    do {\
            (res)[0] -= (a)[0];  (res)[1] -= (a)[1]; \
    }while(0)

#define C_EQ(res, a)\
    do {\
            (res)[0] = (a)[0];  (res)[1] = (a)[1]; \
    }while(0)


#define KISS_FFT_COS(phase) (FFT_SCALAR) cos(phase)
#define KISS_FFT_SIN(phase) (FFT_SCALAR) sin(phase)
#define HALF_OF(x) ((x)*.5)

#define  kf_cexp(x,x_index,phase) \
        do{ \
                (x)(x_index,0) = KISS_FFT_COS(phase);\
                (x)(x_index,1) = KISS_FFT_SIN(phase);\
        }while(0)

namespace LAMMPS_NS {

#define MAXFACTORS 32
/* e.g. an fft of length 128 has 4 factors
 as far as kissfft is concerned: 4*4*4*2  */
template<class DeviceType>
struct kiss_fft_state_kokkos {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  int nfft;
  int inverse;
  typename AT::t_int_64 d_factors;
  typename AT::t_FFT_DATA_1d d_twiddles;
  typename AT::t_FFT_DATA_1d d_scratch;
};

template<class DeviceType>
class KissFFTKokkos {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  KOKKOS_INLINE_FUNCTION
  static void kf_bfly2(typename AT::t_FFT_DATA_1d_um &d_Fout, const size_t fstride,
                       const kiss_fft_state_kokkos<DeviceType> &st, int m, int Fout_count)
  {
      typename AT::t_FFT_DATA_1d_um d_twiddles = st.d_twiddles;
      FFT_SCALAR t[2];
      int Fout2_count;
      int tw1_count = 0;

      Fout2_count = Fout_count + m;
      do {
          //C_FIXDIV(d_Fout[Fout_count],2); C_FIXDIV(d_Fout[Fout2_count],2);

          C_MUL(t,d_Fout,Fout2_count,d_twiddles,tw1_count);
          tw1_count += fstride;
          //C_SUB(*Fout2,*Fout,t);
          d_Fout(Fout2_count,0) = d_Fout(Fout_count,0) - t[0];
          d_Fout(Fout2_count,1) = d_Fout(Fout_count,1) - t[1];
          //C_ADDTO(d_Fout[Fout_count],t);
          d_Fout(Fout_count,0) += t[0];
          d_Fout(Fout_count,1) += t[1];
          ++Fout2_count;
          ++Fout_count;
      } while(--m);
  }

  KOKKOS_INLINE_FUNCTION
  static void kf_bfly4(typename AT::t_FFT_DATA_1d_um &d_Fout, const size_t fstride,
                       const kiss_fft_state_kokkos<DeviceType> &st, const size_t m, int Fout_count)
  {
      typename AT::t_FFT_DATA_1d_um d_twiddles = st.d_twiddles;
      FFT_SCALAR scratch[6][2];
      size_t k=m;
      const size_t m2=2*m;
      const size_t m3=3*m;

      int tw3_count,tw2_count,tw1_count;
      tw3_count = tw2_count = tw1_count = 0;

      do {
          //C_FIXDIV(d_Fout[Fout_count],4); C_FIXDIV(d_Fout[m],4); C_FIXDIV(d_Fout[m2],4); C_FIXDIV(d_Fout[m3],4);

          C_MUL(scratch[0],d_Fout,Fout_count + m,d_twiddles,tw1_count);
          C_MUL(scratch[1],d_Fout,Fout_count + m2,d_twiddles,tw2_count);
          C_MUL(scratch[2],d_Fout,Fout_count + m3,d_twiddles,tw3_count);

          //C_SUB(scratch[5],d_Fout[Fout_count],scratch[1] );
          scratch[5][0] = d_Fout(Fout_count,0) - scratch[1][0];
          scratch[5][1] = d_Fout(Fout_count,1) - scratch[1][1];
          //C_ADDTO(d_Fout[Fout_count], scratch[1]);
          d_Fout(Fout_count,0) += scratch[1][0];
          d_Fout(Fout_count,1) += scratch[1][1];
          //C_ADD(scratch[3],scratch[0],scratch[2]);
          scratch[3][0] = scratch[0][0] + scratch[2][0];
          scratch[3][1] = scratch[0][1] + scratch[2][1];
          //C_SUB( scratch[4] , scratch[0] , scratch[2] );
          scratch[4][0] = scratch[0][0] - scratch[2][0];
          scratch[4][1] = scratch[0][1] - scratch[2][1];
          //C_SUB(d_Fout[m2],d_Fout[Fout_count],scratch[3]);
          d_Fout(Fout_count + m2,0) = d_Fout(Fout_count,0) - scratch[3][0];
          d_Fout(Fout_count + m2,1) = d_Fout(Fout_count,1) - scratch[3][1];

          tw1_count += fstride;
          tw2_count += fstride*2;
          tw3_count += fstride*3;
          //C_ADDTO(d_Fout[Fout_count],scratch[3]);
          d_Fout(Fout_count,0) += scratch[3][0];
          d_Fout(Fout_count,1) += scratch[3][1];

          if (st.inverse) {
              d_Fout(Fout_count + m,0) = scratch[5][0] - scratch[4][1];
              d_Fout(Fout_count + m,1) = scratch[5][1] + scratch[4][0];
              d_Fout(Fout_count + m3,0) = scratch[5][0] + scratch[4][1];
              d_Fout(Fout_count + m3,1) = scratch[5][1] - scratch[4][0];
          } else{
              d_Fout(Fout_count + m,0) = scratch[5][0] + scratch[4][1];
              d_Fout(Fout_count + m,1) = scratch[5][1] - scratch[4][0];
              d_Fout(Fout_count + m3,0) = scratch[5][0] - scratch[4][1];
              d_Fout(Fout_count + m3,1) = scratch[5][1] + scratch[4][0];
          }
          ++Fout_count;
      } while(--k);
  }

  KOKKOS_INLINE_FUNCTION
  static void kf_bfly3(typename AT::t_FFT_DATA_1d_um &d_Fout, const size_t fstride,
                       const kiss_fft_state_kokkos<DeviceType> &st, size_t m, int Fout_count)
  {
      size_t k=m;
      const size_t m2 = 2*m;
      typename AT::t_FFT_DATA_1d_um d_twiddles = st.d_twiddles;
      FFT_SCALAR scratch[5][2];
      FFT_SCALAR epi3[2];
      //C_EQ(epi3,d_twiddles[fstride*m]);
      epi3[0] = d_twiddles(fstride*m,0);
      epi3[1] = d_twiddles(fstride*m,1);

      int tw1_count,tw2_count;
      tw1_count = tw2_count = 0;

      do {
          //C_FIXDIV(d_Fout[Fout_count],3); C_FIXDIV(d_Fout[m],3); C_FIXDIV(d_Fout[m2],3);

          C_MUL(scratch[1],d_Fout,Fout_count + m,d_twiddles,tw1_count);
          C_MUL(scratch[2],d_Fout,Fout_count + m2,d_twiddles,tw2_count);

          //C_ADD(scratch[3],scratch[1],scratch[2]);
          scratch[3][0] = scratch[1][0] + scratch[2][0];
          scratch[3][1] = scratch[1][1] + scratch[2][1];
          //C_SUB(scratch[0],scratch[1],scratch[2]);
          scratch[0][0] = scratch[1][0] - scratch[2][0];
          scratch[0][1] = scratch[1][1] - scratch[2][1];
          tw1_count += fstride;
          tw2_count += fstride*2;

          d_Fout(Fout_count + m,0) = d_Fout(Fout_count,0) - HALF_OF(scratch[3][0]);
          d_Fout(Fout_count + m,1) = d_Fout(Fout_count,1) - HALF_OF(scratch[3][1]);

          //C_MULBYSCALAR(scratch[0],epi3[1]);
          scratch[0][0] *= epi3[1];
          scratch[0][1] *= epi3[1];

          //C_ADDTO(d_Fout[Fout_count],scratch[3]);
          d_Fout(Fout_count,0) += scratch[3][0];
          d_Fout(Fout_count,1) += scratch[3][1];

          d_Fout(Fout_count + m2,0) = d_Fout(Fout_count + m,0) + scratch[0][1];
          d_Fout(Fout_count + m2,1) = d_Fout(Fout_count + m,1) - scratch[0][0];

          d_Fout(Fout_count + m,0) -= scratch[0][1];
          d_Fout(Fout_count + m,1) += scratch[0][0];

          ++Fout_count;
      } while(--k);
  }

  KOKKOS_INLINE_FUNCTION
  static void kf_bfly5(typename AT::t_FFT_DATA_1d_um &d_Fout, const size_t fstride,
                       const kiss_fft_state_kokkos<DeviceType> &st, int m, int Fout_count)
  {
      int u;
      FFT_SCALAR scratch[13][2];
      typename AT::t_FFT_DATA_1d_um d_twiddles = st.d_twiddles;
      FFT_SCALAR ya[2],yb[2];
      //C_EQ(ya,d_twiddles[fstride*m]);
      ya[1] = d_twiddles(fstride*m,1);
      ya[0] = d_twiddles(fstride*m,0);
      //C_EQ(yb,d_twiddles[fstride*2*m]);
      yb[1] = d_twiddles(fstride*2*m,1);
      yb[0] = d_twiddles(fstride*2*m,0);

      int Fout0_count=Fout_count;
      int Fout1_count=Fout0_count+m;
      int Fout2_count=Fout0_count+2*m;
      int Fout3_count=Fout0_count+3*m;
      int Fout4_count=Fout0_count+4*m;

      for ( u=0; u<m; ++u ) {
          //C_FIXDIV( d_Fout[Fout0_count],5); C_FIXDIV( d_Fout[Fout1_count],5); C_FIXDIV( d_Fout[Fout2_count],5);
          //C_FIXDIV( d_Fout[Fout3_count],5); C_FIXDIV( d_Fout[Fout4_count],5);
          //C_EQ(scratch[0],d_Fout[Fout0_count]);
          scratch[0][0] = d_Fout(Fout0_count,0);
          scratch[0][1] = d_Fout(Fout0_count,1);

          C_MUL(scratch[1],d_Fout,Fout1_count,d_twiddles,u*fstride  );
          C_MUL(scratch[2],d_Fout,Fout2_count,d_twiddles,2*u*fstride);
          C_MUL(scratch[3],d_Fout,Fout3_count,d_twiddles,3*u*fstride);
          C_MUL(scratch[4],d_Fout,Fout4_count,d_twiddles,4*u*fstride);

          //C_ADD(scratch[7],scratch[1],scratch[4]);
          scratch[7][0] = scratch[1][0] + scratch[4][0];
          scratch[7][1] = scratch[1][1] + scratch[4][1];
          //C_SUB(scratch[10],scratch[1],scratch[4]);
          scratch[10][0] = scratch[1][0] - scratch[4][0];
          scratch[10][1] = scratch[1][1] - scratch[4][1];
          //C_ADD(scratch[8],scratch[2],scratch[3]);
          scratch[8][0] = scratch[2][0] + scratch[3][0];
          scratch[8][1] = scratch[2][1] + scratch[3][1];
          //C_SUB(scratch[9],scratch[2],scratch[3]);
          scratch[9][0] = scratch[2][0] - scratch[3][0];
          scratch[9][1] = scratch[2][1] - scratch[3][1];

          d_Fout(Fout0_count,0) += scratch[7][0] + scratch[8][0];
          d_Fout(Fout0_count,1) += scratch[7][1] + scratch[8][1];

          scratch[5][0] = scratch[0][0] + S_MUL(scratch[7][0],ya[0]) + S_MUL(scratch[8][0],yb[0]);
          scratch[5][1] = scratch[0][1] + S_MUL(scratch[7][1],ya[0]) + S_MUL(scratch[8][1],yb[0]);

          scratch[6][0] =  S_MUL(scratch[10][1],ya[1]) + S_MUL(scratch[9][1],yb[1]);
          scratch[6][1] = -S_MUL(scratch[10][0],ya[1]) - S_MUL(scratch[9][0],yb[1]);

          //C_SUB(d_Fout[Fout1_count],scratch[5],scratch[6]);
          d_Fout(Fout1_count,0) = scratch[5][0] - scratch[6][0];
          d_Fout(Fout1_count,1) = scratch[5][1] - scratch[6][1];
          //C_ADD(d_Fout[Fout4_count],scratch[5],scratch[6]);
          d_Fout(Fout4_count,0) = scratch[5][0] + scratch[6][0];
          d_Fout(Fout4_count,1) = scratch[5][1] + scratch[6][1];

          scratch[11][0] = scratch[0][0] + S_MUL(scratch[7][0],yb[0]) + S_MUL(scratch[8][0],ya[0]);
          scratch[11][1] = scratch[0][1] + S_MUL(scratch[7][1],yb[0]) + S_MUL(scratch[8][1],ya[0]);
          scratch[12][0] = - S_MUL(scratch[10][1],yb[1]) + S_MUL(scratch[9][1],ya[1]);
          scratch[12][1] = S_MUL(scratch[10][0],yb[1]) - S_MUL(scratch[9][0],ya[1]);

          //C_ADD(d_Fout[Fout2_count],scratch[11],scratch[12]);
          d_Fout(Fout2_count,0) = scratch[11][0] + scratch[12][0];
          d_Fout(Fout2_count,1) = scratch[11][1] + scratch[12][1];
          //C_SUB(d_Fout3[Fout3_count],scratch[11],scratch[12]);
          d_Fout(Fout3_count,0) = scratch[11][0] - scratch[12][0];
          d_Fout(Fout3_count,1) = scratch[11][1] - scratch[12][1];

          ++Fout0_count;++Fout1_count;++Fout2_count;++Fout3_count;++Fout4_count;
      }
  }

  /* perform the butterfly for one stage of a mixed radix FFT */

  KOKKOS_INLINE_FUNCTION
  static void kf_bfly_generic(typename AT::t_FFT_DATA_1d_um &d_Fout, const size_t fstride,
                              const kiss_fft_state_kokkos<DeviceType> &st, int m, int p, int Fout_count)
  {
      int u,k,q1,q;
      typename AT::t_FFT_DATA_1d_um d_twiddles = st.d_twiddles;
      FFT_SCALAR t[2];
      int Norig = st.nfft;

      typename AT::t_FFT_DATA_1d_um d_scratch = st.d_scratch;
      for ( u=0; u<m; ++u ) {
          k=u;
          for ( q1=0 ; q1<p ; ++q1 ) {
              //C_EQ(d_scratch[q1],d_Fout[k]);
              d_scratch(q1,0) = d_Fout(Fout_count + k,0);
              d_scratch(q1,1) = d_Fout(Fout_count + k,1);
              //C_FIXDIV(d_scratch[q1],p);
              k += m;
          }

          k=u;
          for ( q1=0 ; q1<p ; ++q1 ) {
              int twidx=0;
              //C_EQ(d_Fout[k],d_scratch[0]);
              d_Fout(Fout_count + k,0) = d_scratch(0,0);
              d_Fout(Fout_count + k,1) = d_scratch(0,1);
              for (q=1;q<p;++q ) {
                  twidx += fstride * k;
                  if (twidx>=Norig) twidx-=Norig;
                  C_MUL(t,d_scratch,q,d_twiddles,twidx);
                  //C_ADDTO(d_Fout[k],t);
                  d_Fout(Fout_count + k,0) += t[0];
                  d_Fout(Fout_count + k,1) += t[1];
              }
              k += m;
          }
      }
      //KISS_FFT_TMP_FREE(d_scratch);
  }

  KOKKOS_INLINE_FUNCTION
  static void kf_work(typename AT::t_FFT_DATA_1d_um &d_Fout, const typename AT::t_FFT_DATA_1d_um &d_f,
                      const size_t fstride, int in_stride,
                      const typename AT::t_int_64_um &d_factors, const kiss_fft_state_kokkos<DeviceType> &st, int Fout_count, int f_count, int factors_count)
  {
      const int beg = Fout_count;
      const int p = d_factors[factors_count++]; /* the radix  */
      const int m = d_factors[factors_count++]; /* stage's fft length/p */
      const int end = Fout_count + p*m;

      if (m == 1) {
          do {
              //C_EQ(d_Fout[Fout_count],d_f[f_count]);
              d_Fout(Fout_count,0) = d_f(f_count,0);
              d_Fout(Fout_count,1) = d_f(f_count,1);
              f_count += fstride*in_stride;
          } while (++Fout_count != end);
      } else {
          do {
              /* recursive call:
                 DFT of size m*p performed by doing
                 p instances of smaller DFTs of size m,
                 each one takes a decimated version of the input */
              kf_work(d_Fout, d_f, fstride*p, in_stride, d_factors, st, Fout_count, f_count, factors_count);
              f_count += fstride*in_stride;
          } while( (Fout_count += m) != end);
      }

      Fout_count=beg;

      /* recombine the p smaller DFTs */
      switch (p) {
        case 2: kf_bfly2(d_Fout,fstride,st,m,Fout_count); break;
        case 3: kf_bfly3(d_Fout,fstride,st,m,Fout_count); break;
        case 4: kf_bfly4(d_Fout,fstride,st,m,Fout_count); break;
        case 5: kf_bfly5(d_Fout,fstride,st,m,Fout_count); break;
        default: kf_bfly_generic(d_Fout,fstride,st,m,p,Fout_count); break;
      }
  }

  /*  facbuf is populated by p1,m1,p2,m2, ...
      where
      p[i] * m[i] = m[i-1]
      m0 = n                  */

  static int kf_factor(int n, HAT::t_int_64 h_facbuf)
  {
      int p=4, nf=0;
      double floor_sqrt;
      floor_sqrt = floor( sqrt((double)n) );
      int facbuf_count = 0;
      int p_max = 0;

      /* factor out the remaining powers of 4, powers of 2,
         and then any other remaining primes */
      do {
          if (nf == MAXFACTORS) p = n; /* make certain that we don't run out of space */
          while (n % p) {
              switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
              }
              if (p > floor_sqrt)
                  p = n;          /* no more factors, skip to end */
          }
          n /= p;
          h_facbuf[facbuf_count++] = p;
          h_facbuf[facbuf_count++] = n;
          p_max = MAX(p,p_max);
          ++nf;
      } while (n > 1);
      return p_max;
  }

  /*
   * User-callable function to allocate all necessary storage space for the fft.
   *
   * The return value is a contiguous block of memory, allocated with malloc.  As such,
   * It can be freed with free(), rather than a kiss_fft-specific function.
   */

  static kiss_fft_state_kokkos<DeviceType> kiss_fft_alloc_kokkos(int nfft, int inverse_fft, void *mem, size_t *lenmem)
  {
      kiss_fft_state_kokkos<DeviceType> st;
      int i;
      st.nfft = nfft;
      st.inverse = inverse_fft;

      typename AT::tdual_int_64 k_factors = typename AT::tdual_int_64();
      typename AT::tdual_FFT_DATA_1d k_twiddles = typename AT::tdual_FFT_DATA_1d();

      if (nfft > 0) {
          k_factors = typename AT::tdual_int_64("kissfft:factors",MAXFACTORS*2);
          k_twiddles = typename AT::tdual_FFT_DATA_1d("kissfft:twiddles",nfft);

          for (i=0;i<nfft;++i) {
              const double phase = (st.inverse ? 2.0*M_PI:-2.0*M_PI)*i / nfft;
              kf_cexp(k_twiddles.h_view,i,phase );
          }

          int p_max = kf_factor(nfft,k_factors.h_view);
          st.d_scratch = typename AT::t_FFT_DATA_1d("kissfft:scratch",p_max);
      }

      k_factors.template modify<LMPHostType>();
      k_factors.template sync<LMPDeviceType>();
      st.d_factors = k_factors.template view<DeviceType>();

      k_twiddles.template modify<LMPHostType>();
      k_twiddles.template sync<LMPDeviceType>();
      st.d_twiddles = k_twiddles.template view<DeviceType>();

      return st;
  }

  KOKKOS_INLINE_FUNCTION
  static void kiss_fft_stride(const kiss_fft_state_kokkos<DeviceType> &st, const typename AT::t_FFT_DATA_1d_um &d_fin, typename AT::t_FFT_DATA_1d_um &d_fout, int in_stride, int offset)
  {
      //if (d_fin == d_fout) {
      //    // NOTE: this is not really an in-place FFT algorithm.
      //    // It just performs an out-of-place FFT into a temp buffer
      //    typename AT::t_FFT_DATA_1d_um d_tmpbuf = typename AT::t_FFT_DATA_1d("kissfft:tmpbuf",st.nfft);
      //    kf_work(d_tmpbuf,d_fin,1,in_stride,st.d_factors,st,offset,offset,0);
      //    memcpy(d_fout,d_tmpbuf,sizeof(FFT_DATA)*st.nfft);
      //} else {
          kf_work(d_fout,d_fin,1,in_stride,st.d_factors,st,offset,offset,0);
      //}
  }


  KOKKOS_INLINE_FUNCTION
  static void kiss_fft_kokkos(const kiss_fft_state_kokkos<DeviceType> &cfg, const typename AT::t_FFT_DATA_1d_um d_fin, typename AT::t_FFT_DATA_1d_um d_fout, int offset)
  {
      kiss_fft_stride(cfg,d_fin,d_fout,1,offset);
  }

};

}

#endif