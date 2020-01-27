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

/* ----------------------------------------------------------------------
   Contributing authors: Jim Shepherd (GA Tech) added SGI SCSL support
                         Axel Kohlmeyer (Temple U) added support for
                         FFTW3, KISS FFT, Dfti/MKL, and ACML.
                         Phil Blood (PSC) added single precision FFT.
                         Paul Coffman (IBM) added MPI collectives remap
------------------------------------------------------------------------- */

#include "fft3d.h"
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "remap.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef FFT_KISS
/* include kissfft implementation */
#include "kissfft.h"
#endif

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

/* ----------------------------------------------------------------------
   Data layout for 3d FFTs:

   data set of Nfast x Nmid x Nslow elements is owned by P procs
   on input, each proc owns a subsection of the elements
   on output, each proc will own a (possibly different) subsection
   my subsection must not overlap with any other proc's subsection,
     i.e. the union of all proc's input (or output) subsections must
     exactly tile the global Nfast x Nmid x Nslow data set
   when called from C, all subsection indices are
     C-style from 0 to N-1 where N = Nfast or Nmid or Nslow
   when called from F77, all subsection indices are
     F77-style from 1 to N where N = Nfast or Nmid or Nslow
   a proc can own 0 elements on input or output
     by specifying hi index < lo index
   on both input and output, data is stored contiguously on a processor
     with a fast-varying, mid-varying, and slow-varying index
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Perform 3d FFT

   Arguments:
   in           starting address of input data on this proc
   out          starting address of where output data for this proc
                  will be placed (can be same as in)
   flag         1 for forward FFT, -1 for inverse FFT
   plan         plan returned by previous call to fft_3d_create_plan
------------------------------------------------------------------------- */

void fft_3d(FFT_DATA *in, FFT_DATA *out, int flag, struct fft_plan_3d *plan)
{
  int i,total,length,offset,num;
  FFT_SCALAR norm;
#if defined(FFT_FFTW3)
  FFT_SCALAR *out_ptr;
#endif
  FFT_DATA *data,*copy;

  // system specific constants

#if defined(FFT_FFTW3)
  FFTW_API(plan) theplan;
#else
  // nothing to do for other FFTs
#endif

  // pre-remap to prepare for 1st FFTs if needed
  // copy = loc for remap result

  if (plan->pre_plan) {
    if (plan->pre_target == 0) copy = out;
    else copy = plan->copy;
    remap_3d((FFT_SCALAR *) in, (FFT_SCALAR *) copy,
             (FFT_SCALAR *) plan->scratch, plan->pre_plan);
    data = copy;
  }
  else
    data = in;

  // 1d FFTs along fast axis

  total = plan->total1;
  length = plan->length1;

#if defined(FFT_MKL)
  if (flag == -1)
    DftiComputeForward(plan->handle_fast,data);
  else
    DftiComputeBackward(plan->handle_fast,data);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    theplan=plan->plan_fast_forward;
  else
    theplan=plan->plan_fast_backward;
  FFTW_API(execute_dft)(theplan,data,data);
#else
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_fast_forward,&data[offset],&data[offset]);
  else
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_fast_backward,&data[offset],&data[offset]);
#endif

  // 1st mid-remap to prepare for 2nd FFTs
  // copy = loc for remap result

  if (plan->mid1_target == 0) copy = out;
  else copy = plan->copy;
  remap_3d((FFT_SCALAR *) data, (FFT_SCALAR *) copy,
           (FFT_SCALAR *) plan->scratch, plan->mid1_plan);
  data = copy;

  // 1d FFTs along mid axis

  total = plan->total2;
  length = plan->length2;

#if defined(FFT_MKL)
  if (flag == -1)
    DftiComputeForward(plan->handle_mid,data);
  else
    DftiComputeBackward(plan->handle_mid,data);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    theplan=plan->plan_mid_forward;
  else
    theplan=plan->plan_mid_backward;
  FFTW_API(execute_dft)(theplan,data,data);
#else
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_mid_forward,&data[offset],&data[offset]);
  else
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_mid_backward,&data[offset],&data[offset]);
#endif

  // 2nd mid-remap to prepare for 3rd FFTs
  // copy = loc for remap result

  if (plan->mid2_target == 0) copy = out;
  else copy = plan->copy;
  remap_3d((FFT_SCALAR *) data, (FFT_SCALAR *) copy,
           (FFT_SCALAR *) plan->scratch, plan->mid2_plan);
  data = copy;

  // 1d FFTs along slow axis

  total = plan->total3;
  length = plan->length3;

#if defined(FFT_MKL)
  if (flag == -1)
    DftiComputeForward(plan->handle_slow,data);
  else
    DftiComputeBackward(plan->handle_slow,data);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    theplan=plan->plan_slow_forward;
  else
    theplan=plan->plan_slow_backward;
  FFTW_API(execute_dft)(theplan,data,data);
#else
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_slow_forward,&data[offset],&data[offset]);
  else
    for (offset = 0; offset < total; offset += length)
      kiss_fft(plan->cfg_slow_backward,&data[offset],&data[offset]);
#endif

  // post-remap to put data in output format if needed
  // destination is always out

  if (plan->post_plan)
    remap_3d((FFT_SCALAR *) data, (FFT_SCALAR *) out,
             (FFT_SCALAR *) plan->scratch, plan->post_plan);

  // scaling if required
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->normnum;
#if defined(FFT_FFTW3)
    out_ptr = (FFT_SCALAR *)out;
#endif
    for (i = 0; i < num; i++) {
#if defined(FFT_FFTW3)
      *(out_ptr++) *= norm;
      *(out_ptr++) *= norm;
#elif defined(FFT_MKL)
      out[i] *= norm;
#else
      out[i].re *= norm;
      out[i].im *= norm;
#endif
    }
  }
}

/* ----------------------------------------------------------------------
   Create plan for performing a 3d FFT

   Arguments:
   comm                 MPI communicator for the P procs which own the data
   nfast,nmid,nslow     size of global 3d matrix
   in_ilo,in_ihi        input bounds of data I own in fast index
   in_jlo,in_jhi        input bounds of data I own in mid index
   in_klo,in_khi        input bounds of data I own in slow index
   out_ilo,out_ihi      output bounds of data I own in fast index
   out_jlo,out_jhi      output bounds of data I own in mid index
   out_klo,out_khi      output bounds of data I own in slow index
   scaled               0 = no scaling of result, 1 = scaling
   permute              permutation in storage order of indices on output
                          0 = no permutation
                          1 = permute once = mid->fast, slow->mid, fast->slow
                          2 = permute twice = slow->fast, fast->mid, mid->slow
   nbuf                 returns size of internal storage buffers used by FFT
   usecollective        use collective MPI operations for remapping data
------------------------------------------------------------------------- */

struct fft_plan_3d *fft_3d_create_plan(
       MPI_Comm comm, int nfast, int nmid, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int in_klo, int in_khi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int out_klo, int out_khi,
       int scaled, int permute, int *nbuf, int usecollective)
{
  struct fft_plan_3d *plan;
  int me,nprocs,nthreads;
  int flag,remapflag;
  int first_ilo,first_ihi,first_jlo,first_jhi,first_klo,first_khi;
  int second_ilo,second_ihi,second_jlo,second_jhi,second_klo,second_khi;
  int third_ilo,third_ihi,third_jlo,third_jhi,third_klo,third_khi;
  int out_size,first_size,second_size,third_size,copy_size,scratch_size;
  int np1,np2,ip1,ip2;

  // query MPI info

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

#if defined(_OPENMP)
  // query OpenMP info.
  // should have been initialized systemwide in Comm class constructor
  nthreads = omp_get_max_threads();
#else
  nthreads = 1;
#endif

  // compute division of procs in 2 dimensions not on-processor

  bifactor(nprocs,&np1,&np2);
  ip1 = me % np1;
  ip2 = me/np1;

  // allocate memory for plan data struct

  plan = (struct fft_plan_3d *) malloc(sizeof(struct fft_plan_3d));
  if (plan == NULL) return NULL;

  // remap from initial distribution to layout needed for 1st set of 1d FFTs
  // not needed if all procs own entire fast axis initially
  // first indices = distribution after 1st set of FFTs

  if (in_ilo == 0 && in_ihi == nfast-1)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    first_ilo = in_ilo;
    first_ihi = in_ihi;
    first_jlo = in_jlo;
    first_jhi = in_jhi;
    first_klo = in_klo;
    first_khi = in_khi;
    plan->pre_plan = NULL;
  } else {
    first_ilo = 0;
    first_ihi = nfast - 1;
    first_jlo = ip1*nmid/np1;
    first_jhi = (ip1+1)*nmid/np1 - 1;
    first_klo = ip2*nslow/np2;
    first_khi = (ip2+1)*nslow/np2 - 1;
    plan->pre_plan =
      remap_3d_create_plan(comm,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
                           first_ilo,first_ihi,first_jlo,first_jhi,
                           first_klo,first_khi,2,0,0,FFT_PRECISION,0);
    if (plan->pre_plan == NULL) return NULL;
  }

  // 1d FFTs along fast axis

  plan->length1 = nfast;
  plan->total1 = nfast * (first_jhi-first_jlo+1) * (first_khi-first_klo+1);

  // remap from 1st to 2nd FFT
  // choose which axis is split over np1 vs np2 to minimize communication
  // second indices = distribution after 2nd set of FFTs

  second_ilo = ip1*nfast/np1;
  second_ihi = (ip1+1)*nfast/np1 - 1;
  second_jlo = 0;
  second_jhi = nmid - 1;
  second_klo = ip2*nslow/np2;
  second_khi = (ip2+1)*nslow/np2 - 1;
  plan->mid1_plan =
      remap_3d_create_plan(comm,
                           first_ilo,first_ihi,first_jlo,first_jhi,
                           first_klo,first_khi,
                           second_ilo,second_ihi,second_jlo,second_jhi,
                           second_klo,second_khi,2,1,0,FFT_PRECISION,
                           usecollective);
  if (plan->mid1_plan == NULL) return NULL;

  // 1d FFTs along mid axis

  plan->length2 = nmid;
  plan->total2 = (second_ihi-second_ilo+1) * nmid * (second_khi-second_klo+1);

  // remap from 2nd to 3rd FFT
  // if final distribution is permute=2 with all procs owning entire slow axis
  //   then this remapping goes directly to final distribution
  //  third indices = distribution after 3rd set of FFTs

  if (permute == 2 && out_klo == 0 && out_khi == nslow-1)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0) {
    third_ilo = out_ilo;
    third_ihi = out_ihi;
    third_jlo = out_jlo;
    third_jhi = out_jhi;
    third_klo = out_klo;
    third_khi = out_khi;
  } else {
    third_ilo = ip1*nfast/np1;
    third_ihi = (ip1+1)*nfast/np1 - 1;
    third_jlo = ip2*nmid/np2;
    third_jhi = (ip2+1)*nmid/np2 - 1;
    third_klo = 0;
    third_khi = nslow - 1;
  }

  plan->mid2_plan =
    remap_3d_create_plan(comm,
                         second_jlo,second_jhi,second_klo,second_khi,
                         second_ilo,second_ihi,
                         third_jlo,third_jhi,third_klo,third_khi,
                         third_ilo,third_ihi,2,1,0,FFT_PRECISION,usecollective);
  if (plan->mid2_plan == NULL) return NULL;

  // 1d FFTs along slow axis

  plan->length3 = nslow;
  plan->total3 = (third_ihi-third_ilo+1) * (third_jhi-third_jlo+1) * nslow;

  // remap from 3rd FFT to final distribution
  //  not needed if permute = 2 and third indices = out indices on all procs

  if (permute == 2 &&
      out_ilo == third_ilo && out_ihi == third_ihi &&
      out_jlo == third_jlo && out_jhi == third_jhi &&
      out_klo == third_klo && out_khi == third_khi)
    flag = 0;
  else
    flag = 1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0)
    plan->post_plan = NULL;
  else {
    plan->post_plan =
      remap_3d_create_plan(comm,
                           third_klo,third_khi,third_ilo,third_ihi,
                           third_jlo,third_jhi,
                           out_klo,out_khi,out_ilo,out_ihi,
                           out_jlo,out_jhi,2,(permute+1)%3,0,FFT_PRECISION,0);
    if (plan->post_plan == NULL) return NULL;
  }

  // configure plan memory pointers and allocate work space
  // out_size = amount of memory given to FFT by user
  // first/second/third_size =
  //   amount of memory needed after pre,mid1,mid2 remaps
  // copy_size = amount needed internally for extra copy of data
  // scratch_size = amount needed internally for remap scratch space
  // for each remap:
  //   out space used for result if big enough, else require copy buffer
  //   accumulate largest required remap scratch space

  out_size = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) * (out_khi-out_klo+1);
  first_size = (first_ihi-first_ilo+1) * (first_jhi-first_jlo+1) *
    (first_khi-first_klo+1);
  second_size = (second_ihi-second_ilo+1) * (second_jhi-second_jlo+1) *
    (second_khi-second_klo+1);
  third_size = (third_ihi-third_ilo+1) * (third_jhi-third_jlo+1) *
    (third_khi-third_klo+1);

  copy_size = 0;
  scratch_size = 0;

  if (plan->pre_plan) {
    if (first_size <= out_size)
      plan->pre_target = 0;
    else {
      plan->pre_target = 1;
      copy_size = MAX(copy_size,first_size);
    }
    scratch_size = MAX(scratch_size,first_size);
  }

  if (plan->mid1_plan) {
    if (second_size <= out_size)
      plan->mid1_target = 0;
    else {
      plan->mid1_target = 1;
      copy_size = MAX(copy_size,second_size);
    }
    scratch_size = MAX(scratch_size,second_size);
  }

  if (plan->mid2_plan) {
    if (third_size <= out_size)
      plan->mid2_target = 0;
    else {
      plan->mid2_target = 1;
      copy_size = MAX(copy_size,third_size);
    }
    scratch_size = MAX(scratch_size,third_size);
  }

  if (plan->post_plan)
    scratch_size = MAX(scratch_size,out_size);

  *nbuf = copy_size + scratch_size;

  if (copy_size) {
    plan->copy = (FFT_DATA *) malloc(copy_size*sizeof(FFT_DATA));
    if (plan->copy == NULL) return NULL;
  }
  else plan->copy = NULL;

  if (scratch_size) {
    plan->scratch = (FFT_DATA *) malloc(scratch_size*sizeof(FFT_DATA));
    if (plan->scratch == NULL) return NULL;
  }
  else plan->scratch = NULL;

  // system specific pre-computation of 1d FFT coeffs
  // and scaling normalization

#if defined(FFT_MKL)
  DftiCreateDescriptor( &(plan->handle_fast), FFT_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total1/nfast);
  DftiSetValue(plan->handle_fast, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_fast, DFTI_INPUT_DISTANCE, (MKL_LONG)nfast);
  DftiSetValue(plan->handle_fast, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nfast);
#if defined(FFT_MKL_THREADS)
  DftiSetValue(plan->handle_fast, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_fast);

  DftiCreateDescriptor( &(plan->handle_mid), FFT_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nmid);
  DftiSetValue(plan->handle_mid, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total2/nmid);
  DftiSetValue(plan->handle_mid, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_mid, DFTI_INPUT_DISTANCE, (MKL_LONG)nmid);
  DftiSetValue(plan->handle_mid, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nmid);
#if defined(FFT_MKL_THREADS)
  DftiSetValue(plan->handle_mid, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_mid);

  DftiCreateDescriptor( &(plan->handle_slow), FFT_MKL_PREC, DFTI_COMPLEX, 1,
                        (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_NUMBER_OF_TRANSFORMS,
               (MKL_LONG)plan->total3/nslow);
  DftiSetValue(plan->handle_slow, DFTI_PLACEMENT,DFTI_INPLACE);
  DftiSetValue(plan->handle_slow, DFTI_INPUT_DISTANCE, (MKL_LONG)nslow);
  DftiSetValue(plan->handle_slow, DFTI_OUTPUT_DISTANCE, (MKL_LONG)nslow);
#if defined(FFT_MKL_THREADS)
  DftiSetValue(plan->handle_slow, DFTI_NUMBER_OF_USER_THREADS, nthreads);
#endif
  DftiCommitDescriptor(plan->handle_slow);

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

#elif defined(FFT_FFTW3)
#if defined(FFT_FFTW_THREADS)
  if (nthreads > 1) {
    FFTW_API(init_threads)();
    FFTW_API(plan_with_nthreads)(nthreads);
  }
#endif

  plan->plan_fast_forward =
    FFTW_API(plan_many_dft)(1, &nfast,plan->total1/plan->length1,
                            NULL,&nfast,1,plan->length1,
                            NULL,&nfast,1,plan->length1,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  plan->plan_fast_backward =
    FFTW_API(plan_many_dft)(1, &nfast,plan->total1/plan->length1,
                            NULL,&nfast,1,plan->length1,
                            NULL,&nfast,1,plan->length1,
                            FFTW_BACKWARD,FFTW_ESTIMATE);
  plan->plan_mid_forward =
    FFTW_API(plan_many_dft)(1, &nmid,plan->total2/plan->length2,
                            NULL,&nmid,1,plan->length2,
                            NULL,&nmid,1,plan->length2,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  plan->plan_mid_backward =
    FFTW_API(plan_many_dft)(1, &nmid,plan->total2/plan->length2,
                            NULL,&nmid,1,plan->length2,
                            NULL,&nmid,1,plan->length2,
                            FFTW_BACKWARD,FFTW_ESTIMATE);
  plan->plan_slow_forward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total3/plan->length3,
                            NULL,&nslow,1,plan->length3,
                            NULL,&nslow,1,plan->length3,
                            FFTW_FORWARD,FFTW_ESTIMATE);
  plan->plan_slow_backward =
    FFTW_API(plan_many_dft)(1, &nslow,plan->total3/plan->length3,
                            NULL,&nslow,1,plan->length3,
                            NULL,&nslow,1,plan->length3,
                            FFTW_BACKWARD,FFTW_ESTIMATE);

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }
#else
  plan->cfg_fast_forward = kiss_fft_alloc(nfast,0,NULL,NULL);
  plan->cfg_fast_backward = kiss_fft_alloc(nfast,1,NULL,NULL);

  if (nmid == nfast) {
    plan->cfg_mid_forward = plan->cfg_fast_forward;
    plan->cfg_mid_backward = plan->cfg_fast_backward;
  }
  else {
    plan->cfg_mid_forward = kiss_fft_alloc(nmid,0,NULL,NULL);
    plan->cfg_mid_backward = kiss_fft_alloc(nmid,1,NULL,NULL);
  }

  if (nslow == nfast) {
    plan->cfg_slow_forward = plan->cfg_fast_forward;
    plan->cfg_slow_backward = plan->cfg_fast_backward;
  }
  else if (nslow == nmid) {
    plan->cfg_slow_forward = plan->cfg_mid_forward;
    plan->cfg_slow_backward = plan->cfg_mid_backward;
  }
  else {
    plan->cfg_slow_forward = kiss_fft_alloc(nslow,0,NULL,NULL);
    plan->cfg_slow_backward = kiss_fft_alloc(nslow,1,NULL,NULL);
  }

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

#endif

  return plan;
}

/* ----------------------------------------------------------------------
   Destroy a 3d fft plan
------------------------------------------------------------------------- */

void fft_3d_destroy_plan(struct fft_plan_3d *plan)
{
  if (plan->pre_plan) remap_3d_destroy_plan(plan->pre_plan);
  if (plan->mid1_plan) remap_3d_destroy_plan(plan->mid1_plan);
  if (plan->mid2_plan) remap_3d_destroy_plan(plan->mid2_plan);
  if (plan->post_plan) remap_3d_destroy_plan(plan->post_plan);

  if (plan->copy) free(plan->copy);
  if (plan->scratch) free(plan->scratch);

#if defined(FFT_MKL)
  DftiFreeDescriptor(&(plan->handle_fast));
  DftiFreeDescriptor(&(plan->handle_mid));
  DftiFreeDescriptor(&(plan->handle_slow));
#elif defined(FFT_FFTW3)
  FFTW_API(destroy_plan)(plan->plan_slow_forward);
  FFTW_API(destroy_plan)(plan->plan_slow_backward);
  FFTW_API(destroy_plan)(plan->plan_mid_forward);
  FFTW_API(destroy_plan)(plan->plan_mid_backward);
  FFTW_API(destroy_plan)(plan->plan_fast_forward);
  FFTW_API(destroy_plan)(plan->plan_fast_backward);
#if defined(FFT_FFTW_THREADS)
  FFTW_API(cleanup_threads)();
#endif
#else
  if (plan->cfg_slow_forward != plan->cfg_fast_forward &&
      plan->cfg_slow_forward != plan->cfg_mid_forward) {
    free(plan->cfg_slow_forward);
    free(plan->cfg_slow_backward);
  }
  if (plan->cfg_mid_forward != plan->cfg_fast_forward) {
    free(plan->cfg_mid_forward);
    free(plan->cfg_mid_backward);
  }
  free(plan->cfg_fast_forward);
  free(plan->cfg_fast_backward);
#endif

  free(plan);
}

/* ----------------------------------------------------------------------
   recursively divide n into small factors, return them in list
------------------------------------------------------------------------- */

void factor(int n, int *num, int *list)
{
  if (n == 1) {
    return;
  } else if (n % 2 == 0) {
    *list = 2;
    (*num)++;
    factor(n/2,num,list+1);
  } else if (n % 3 == 0) {
    *list = 3;
    (*num)++;
    factor(n/3,num,list+1);
  } else if (n % 5 == 0) {
    *list = 5;
    (*num)++;
    factor(n/5,num,list+1);
  } else if (n % 7 == 0) {
    *list = 7;
    (*num)++;
    factor(n/7,num,list+1);
  } else if (n % 11 == 0) {
    *list = 11;
    (*num)++;
    factor(n/11,num,list+1);
  } else if (n % 13 == 0) {
    *list = 13;
    (*num)++;
    factor(n/13,num,list+1);
  } else {
    *list = n;
    (*num)++;
    return;
  }
}

/* ----------------------------------------------------------------------
   divide n into 2 factors of as equal size as possible
------------------------------------------------------------------------- */

void bifactor(int n, int *factor1, int *factor2)
{
  int n1,n2,facmax;

  facmax = static_cast<int> (sqrt((double) n));

  for (n1 = facmax; n1 > 0; n1--) {
    n2 = n/n1;
    if (n1*n2 == n) {
      *factor1 = n1;
      *factor2 = n2;
      return;
    }
  }
}

/* ----------------------------------------------------------------------
   perform just the 1d FFTs needed by a 3d FFT, no data movement
   used for timing purposes

   Arguments:
   in           starting address of input data on this proc, all set to 0.0
   nsize        size of in
   flag         1 for forward FFT, -1 for inverse FFT
   plan         plan returned by previous call to fft_3d_create_plan
------------------------------------------------------------------------- */

void fft_1d_only(FFT_DATA *data, int nsize, int flag, struct fft_plan_3d *plan)
{
  int i,num;
  FFT_SCALAR norm;
#if defined(FFT_FFTW3)
  FFT_SCALAR *data_ptr;
#endif

  // total = size of data needed in each dim
  // length = length of 1d FFT in each dim
  // total/length = # of 1d FFTs in each dim
  // if total > nsize, limit # of 1d FFTs to available size of data

  int total1 = plan->total1;
  int length1 = plan->length1;
  int total2 = plan->total2;
  int length2 = plan->length2;
  int total3 = plan->total3;
  int length3 = plan->length3;

// fftw3 and Dfti in MKL encode the number of transforms
// into the plan, so we cannot operate on a smaller data set.
#if defined(FFT_MKL) || defined(FFT_FFTW3)
  if ((total1 > nsize) || (total2 > nsize) || (total3 > nsize))
    return;
#endif
  if (total1 > nsize) total1 = (nsize/length1) * length1;
  if (total2 > nsize) total2 = (nsize/length2) * length2;
  if (total3 > nsize) total3 = (nsize/length3) * length3;

  // perform 1d FFTs in each of 3 dimensions
  // data is just an array of 0.0

#if defined(FFT_MKL)
  if (flag == -1) {
    DftiComputeForward(plan->handle_fast,data);
    DftiComputeForward(plan->handle_mid,data);
    DftiComputeForward(plan->handle_slow,data);
  } else {
    DftiComputeBackward(plan->handle_fast,data);
    DftiComputeBackward(plan->handle_mid,data);
    DftiComputeBackward(plan->handle_slow,data);
  }
#elif defined(FFT_FFTW3)
  FFTW_API(plan) theplan;
  if (flag == -1)
    theplan=plan->plan_fast_forward;
  else
    theplan=plan->plan_fast_backward;
  FFTW_API(execute_dft)(theplan,data,data);
  if (flag == -1)
    theplan=plan->plan_mid_forward;
  else
    theplan=plan->plan_mid_backward;
  FFTW_API(execute_dft)(theplan,data,data);
  if (flag == -1)
    theplan=plan->plan_slow_forward;
  else
    theplan=plan->plan_slow_backward;
  FFTW_API(execute_dft)(theplan,data,data);
#else
  if (flag == -1) {
    for (int offset = 0; offset < total1; offset += length1)
      kiss_fft(plan->cfg_fast_forward,&data[offset],&data[offset]);
    for (int offset = 0; offset < total2; offset += length2)
      kiss_fft(plan->cfg_mid_forward,&data[offset],&data[offset]);
    for (int offset = 0; offset < total3; offset += length3)
      kiss_fft(plan->cfg_slow_forward,&data[offset],&data[offset]);
  } else {
    for (int offset = 0; offset < total1; offset += length1)
      kiss_fft(plan->cfg_fast_backward,&data[offset],&data[offset]);
    for (int offset = 0; offset < total2; offset += length2)
      kiss_fft(plan->cfg_mid_backward,&data[offset],&data[offset]);
    for (int offset = 0; offset < total3; offset += length3)
      kiss_fft(plan->cfg_slow_backward,&data[offset],&data[offset]);
  }
#endif

  // scaling if required
  // limit num to size of data

  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = MIN(plan->normnum,nsize);
#if defined(FFT_FFTW3)
    data_ptr = (FFT_SCALAR *)data;
#endif
    for (i = 0; i < num; i++) {
#if defined(FFT_FFTW3)
      *(data_ptr++) *= norm;
      *(data_ptr++) *= norm;
#elif defined(FFT_MKL)
      data[i] *= norm;
#else
      data[i].re *= norm;
      data[i].im *= norm;
#endif
    }
  }
}
