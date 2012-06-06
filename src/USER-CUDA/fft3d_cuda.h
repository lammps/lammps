/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

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

// User-settable FFT precision

// FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag)
#include "cuda_precision.h"
//#define FFT_PRECISION 2

// -------------------------------------------------------------------------

// Data types for single-precision complex

#if FFT_PRECISION_CU == 1

#ifdef FFT_CUFFT
#include "cuda_runtime.h"
#include "cufft.h"
typedef struct {
  float re;
  float im;
} FFT_DATA;
typedef cufftComplex cufftData;
typedef cufftReal cufftDataInit;
#define cufft cufftExecC2C
#define cufftinit cufftExecR2C
#define CUFFT_PLAN CUFFT_C2C
#define CUFFT_PLAN_INIT CUFFT_R2C
#else
typedef struct {
  float re;
  float im;
} FFT_DATA;
#endif

#endif

// -------------------------------------------------------------------------

// Data types for double-precision complex

#if FFT_PRECISION_CU == 2


#ifdef FFT_CUFFT
#include "cuda_runtime.h"
#include "cufft.h"
typedef cufftDoubleComplex cufftData;
typedef cufftDoubleReal cufftDataInit;
typedef struct {
  double re;
  double im;
} FFT_DATA;
#define cufft cufftExecZ2Z
#define cufftinit cufftExecD2Z
#define CUFFT_PLAN CUFFT_Z2Z
#define CUFFT_PLAN_INIT CUFFT_D2Z
#endif

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

  double coretime;
  double ffttime;
  int iterate;
                                    // system specific 1d FFT info

#ifdef FFT_CUFFT
  //CUdeviceptr cudata;
  cufftData* cudata;
  cufftData* cudata2;
  unsigned int cudatasize;
  cufftHandle plan_fast;
  cufftHandle plan_mid;
  cufftHandle plan_slow;
  cufftHandle plan_3d;
  int nfast;
  int nmid;
  int nslow;
  int ihi_out,ilo_out,jhi_out,jlo_out,khi_out,klo_out;
  int me,nprocs;
#endif
  int init;
};

// function prototypes

void fft_3d_destroy_plan_cuda(struct fft_plan_3d *);
void factor_cuda(int, int *, int *);
void bifactor_cuda(int, int *, int *);
void fft_1d_only_cuda(FFT_DATA *, int, int, struct fft_plan_3d *);
void fft_3d_cudaA(FFT_DATA *, FFT_DATA *, int, struct fft_plan_3d *);
void fft_3d_cuda(FFT_DATA *, FFT_DATA *, int, struct fft_plan_3d *);
struct fft_plan_3d *fft_3d_create_plan_cuda(MPI_Comm, int, int, int,
  int, int, int, int, int, int, int, int, int, int, int, int,
  int, int, int *,bool init);
