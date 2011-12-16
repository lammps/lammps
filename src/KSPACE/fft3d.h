/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// User-settable FFT precision 

// FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag) 
// FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag) 

#ifdef FFT_SINGLE
#define FFT_PRECISION 1
typedef float FFT_SCALAR;
#else
#define FFT_PRECISION 2
typedef double FFT_SCALAR;
#endif


// set default fftw library. switch to FFT_FFTW3 when convenient.
#ifdef FFT_FFTW
#define FFT_FFTW2
#endif

// ------------------------------------------------------------------------- 

// Data types for single-precision complex 

#if FFT_PRECISION == 1

#if defined(FFT_SGI)
#include "fft.h"
typedef complex FFT_DATA;
#define FFT_1D cfft1d
#define FFT_1D_INIT cfft1di
extern "C" {
  int cfft1d(int, int, FFT_DATA *, int, FFT_DATA *);
  FFT_DATA *cfft1di(int, FFT_DATA *);
}

#elif defined(FFT_SCSL)
#include <scsl_fft.h>
typedef scsl_complex FFT_DATA;
typedef float FFT_PREC;
#define FFT_1D ccfft
#define FFT_1D_INIT ccfft
extern "C" {
  int ccfft(int, int, FFT_PREC, FFT_DATA *, FFT_DATA *,
                      FFT_PREC *, FFT_PREC *, int *);
}

#elif defined(FFT_ACML)
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft1m_
extern "C" {
  void cfft1m_(int *, int *, int *, FFT_DATA *, FFT_DATA *, int *);
}

#elif defined(FFT_INTEL)
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft1d_
#define FFT_1D_INIT cfft1d_
extern "C" {
  void cfft1d_(FFT_DATA *, int *, int *, FFT_DATA *);
}

#elif defined(FFT_MKL)
#include "mkl_dfti.h"
typedef float _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_SINGLE

#elif defined(FFT_DEC)
typedef struct {
  float re;
  float im;
} FFT_DATA;
#define FFT_1D cfft_
extern "C" {
  void cfft_(char *, char *, char *, FFT_DATA *, FFT_DATA *, int *, int *);
}

#elif defined(FFT_T3E)
#include <complex.h>
typedef complex single FFT_DATA;
#define FFT_1D GGFFT
#define FFT_1D_INIT GGFFT
extern "C" {
  void GGFFT(int *, int *, double *, FFT_DATA *, FFT_DATA *,
	     double *, double *, int *);
}

#elif defined(FFT_FFTW2)
#if defined(FFTW_SIZE)
#include "sfftw.h"
#else
#include "fftw.h"
#endif
typedef FFTW_COMPLEX FFT_DATA;

#elif defined(FFT_FFTW3)
#include "fftw3.h"
typedef fftwf_complex FFT_DATA;
#define FFTW_API(function)  fftwf_ ## function

#else

/* use a stripped down version of kiss fft as default fft */
#ifndef FFT_KISSFFT
#define FFT_KISSFFT
#endif
#define kiss_fft_scalar float
typedef struct {
    kiss_fft_scalar re;
    kiss_fft_scalar im;
} FFT_DATA;

struct kiss_fft_state;
typedef struct kiss_fft_state* kiss_fft_cfg;
#endif

// ------------------------------------------------------------------------- 

// Data types for double-precision complex 

#elif FFT_PRECISION == 2

#if defined(FFT_SGI)
#include "fft.h"
typedef zomplex FFT_DATA;
#define FFT_1D zfft1d
#define FFT_1D_INIT zfft1di
extern "C" {
  int zfft1d(int, int, FFT_DATA *, int, FFT_DATA *);
  FFT_DATA *zfft1di(int, FFT_DATA *);
}

#elif defined(FFT_SCSL)
#include <scsl_fft.h>
typedef scsl_zomplex FFT_DATA;
typedef double FFT_PREC;
#define FFT_1D zzfft
#define FFT_1D_INIT zzfft
extern "C" {
  int zzfft(int, int, FFT_PREC, FFT_DATA *, FFT_DATA *,
                      FFT_PREC *, FFT_PREC *, int *);
}

#elif defined(FFT_ACML)
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft1m_
extern "C" {
  void zfft1m_(int *, int *, int *, FFT_DATA *, FFT_DATA *, int *);
}

#elif defined(FFT_INTEL)
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft1d_
#define FFT_1D_INIT zfft1d_
extern "C" {
  void zfft1d_(FFT_DATA *, int *, int *, FFT_DATA *);
}

#elif defined(FFT_MKL)
#include "mkl_dfti.h"
typedef double _Complex FFT_DATA;
#define FFT_MKL_PREC DFTI_DOUBLE

#elif defined(FFT_DEC)
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define FFT_1D zfft_
extern "C" {
  void zfft_(char *, char *, char *, FFT_DATA *, FFT_DATA *, int *, int *);
}

#elif defined(FFT_T3E)
#include <complex.h>
typedef complex double FFT_DATA;
#define FFT_1D CCFFT
#define FFT_1D_INIT CCFFT
extern "C" {
  void CCFFT(int *, int *, double *, FFT_DATA *, FFT_DATA *,
	     double *, double *, int *);
}

#elif defined(FFT_FFTW2)
#if defined(FFTW_SIZE)
#include "dfftw.h"
#else
#include "fftw.h"
#endif
typedef FFTW_COMPLEX FFT_DATA;

#elif defined(FFT_FFTW3)
#include "fftw3.h"
typedef fftw_complex FFT_DATA;
#define FFTW_API(function)  fftw_ ## function

#else

/* use a stripped down version of kiss fft as default fft */
#ifndef FFT_KISSFFT
#define FFT_KISSFFT
#endif
#define kiss_fft_scalar double
typedef struct {
    kiss_fft_scalar re;
    kiss_fft_scalar im;
} FFT_DATA;

struct kiss_fft_state;
typedef struct kiss_fft_state* kiss_fft_cfg;
#endif

#else
#error "FFT_PRECISION needs to be either 1 (=single) or 2 (=double)"
#endif

// ------------------------------------------------------------------------- 

// details of how to do a 3d FFT 

struct fft_plan_3d {
  struct remap_plan_3d *pre_plan;       // remap from input -> 1st FFTs 
  struct remap_plan_3d *mid1_plan;      // remap from 1st -> 2nd FFTs 
  struct remap_plan_3d *mid2_plan;      // remap from 2nd -> 3rd FFTs 
  struct remap_plan_3d *post_plan;      // remap from 3rd FFTs -> output 
  FFT_DATA *copy;                   // memory for remap results (if needed) 
  FFT_DATA *scratch;                // scratch space for remaps 
  int total1,total2,total3;         // # of 1st,2nd,3rd FFTs (times length) 
  int length1,length2,length3;      // length of 1st,2nd,3rd FFTs 
  int pre_target;                   // where to put remap results 
  int mid1_target,mid2_target;
  int scaled;                       // whether to scale FFT results 
  int normnum;                      // # of values to rescale 
  double norm;                      // normalization factor for rescaling 

                                    // system specific 1d FFT info 
#if defined(FFT_SGI)
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#elif defined(FFT_SCSL)
  FFT_PREC *coeff1;
  FFT_PREC *coeff2;
  FFT_PREC *coeff3;
  FFT_PREC *work1;
  FFT_PREC *work2;
  FFT_PREC *work3;
#elif defined(FFT_ACML)
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#elif defined(FFT_INTEL)
  FFT_DATA *coeff1;
  FFT_DATA *coeff2;
  FFT_DATA *coeff3;
#elif defined(FFT_MKL)
  DFTI_DESCRIPTOR *handle_fast;
  DFTI_DESCRIPTOR *handle_mid;
  DFTI_DESCRIPTOR *handle_slow;
#elif defined(FFT_T3E)
  double *coeff1;
  double *coeff2;
  double *coeff3;
  double *work1;
  double *work2;
  double *work3;
#elif defined(FFT_FFTW2)
  fftw_plan plan_fast_forward;
  fftw_plan plan_fast_backward;
  fftw_plan plan_mid_forward;
  fftw_plan plan_mid_backward;
  fftw_plan plan_slow_forward;
  fftw_plan plan_slow_backward;
#elif defined(FFT_FFTW3)
  FFTW_API(plan) plan_fast_forward;
  FFTW_API(plan) plan_fast_backward;
  FFTW_API(plan) plan_mid_forward;
  FFTW_API(plan) plan_mid_backward;
  FFTW_API(plan) plan_slow_forward;
  FFTW_API(plan) plan_slow_backward;
#elif defined(FFT_KISSFFT)
  kiss_fft_cfg cfg_fast_forward;
  kiss_fft_cfg cfg_fast_backward;
  kiss_fft_cfg cfg_mid_forward;
  kiss_fft_cfg cfg_mid_backward;
  kiss_fft_cfg cfg_slow_forward;
  kiss_fft_cfg cfg_slow_backward;
#endif
};

// function prototypes 

void fft_3d(FFT_DATA *, FFT_DATA *, int, struct fft_plan_3d *);
struct fft_plan_3d *fft_3d_create_plan(MPI_Comm, int, int, int,
  int, int, int, int, int, int, int, int, int, int, int, int,
  int, int, int *);
void fft_3d_destroy_plan(struct fft_plan_3d *);
void factor(int, int *, int *);
void bifactor(int, int *, int *);
void fft_1d_only(FFT_DATA *, int, int, struct fft_plan_3d *);
/* ERROR/WARNING messages:

*/
