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
                         Axel Kohlmeyer (Temple U) added FFTW3, KISSFFT support
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "fft3d.h"
#include "remap.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

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
  double norm;
  FFT_DATA *data,*copy;

  // system specific constants 

#if defined(FFT_SCSL)
  int isys = 0;
  FFT_PREC scalef = 1.0;
#elif defined(FFT_DEC)
  char c = 'C';
  char f = 'F';
  char b = 'B';
  int one = 1;
#elif defined(FFT_T3E)
  int isys = 0;
  double scalef = 1.0;
#else
  // nothing to do for other FFTs.
#endif

  // pre-remap to prepare for 1st FFTs if needed
  // copy = loc for remap result 

  if (plan->pre_plan) {
    if (plan->pre_target == 0) copy = out;
    else copy = plan->copy;
    remap_3d((double *) in, (double *) copy, (double *) plan->scratch,
	     plan->pre_plan);
    data = copy;
  }
  else
    data = in;

  // 1d FFTs along fast axis 

  total = plan->total1;
  length = plan->length1;

#if defined(FFT_SGI)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff1);
#elif defined(FFT_SCSL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,scalef,&data[offset],&data[offset],plan->coeff1,
	   plan->work1,&isys);
#elif defined(FFT_INTEL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff1);
#elif defined(FFT_DEC)
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#elif defined(FFT_T3E)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff1,
	   plan->work1,&isys);
#elif defined(FFT_FFTW2)
  if (flag == -1)
    fftw(plan->plan_fast_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_fast_backward,total/length,data,1,length,NULL,0,0);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    FFTW_API(execute_dft)(plan->plan_fast_forward,data,data);
  else
    FFTW_API(execute_dft)(plan->plan_fast_backward,data,data);
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
  remap_3d((double *) data, (double *) copy, (double *) plan->scratch,
	   plan->mid1_plan);
  data = copy;

  // 1d FFTs along mid axis 

  total = plan->total2;
  length = plan->length2;

#if defined(FFT_SGI)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff2);
#elif defined(FFT_SCSL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,scalef,&data[offset],&data[offset],plan->coeff2,
	   plan->work2,&isys);
#elif defined(FFT_INTEL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff2);
#elif defined(FFT_DEC)
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#elif defined(FFT_T3E)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff2,
	   plan->work2,&isys);
#elif defined(FFT_FFTW2)
  if (flag == -1)
    fftw(plan->plan_mid_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_mid_backward,total/length,data,1,length,NULL,0,0);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    FFTW_API(execute_dft)(plan->plan_mid_forward,data,data);
  else
    FFTW_API(execute_dft)(plan->plan_mid_backward,data,data);
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
  remap_3d((double *) data, (double *) copy, (double *) plan->scratch,
	   plan->mid2_plan);
  data = copy;

  // 1d FFTs along slow axis 

  total = plan->total3;
  length = plan->length3;

#if defined(FFT_SGI)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,&data[offset],1,plan->coeff3);
#elif defined(FFT_SCSL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(flag,length,scalef,&data[offset],&data[offset],plan->coeff3,
	   plan->work3,&isys);
#elif defined(FFT_INTEL)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&data[offset],&length,&flag,plan->coeff3);
#elif defined(FFT_DEC)
  if (flag == -1)
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length,&one);
  else
    for (offset = 0; offset < total; offset += length)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length,&one);
#elif defined(FFT_T3E)
  for (offset = 0; offset < total; offset += length)
    FFT_1D(&flag,&length,&scalef,&data[offset],&data[offset],plan->coeff3,
	   plan->work3,&isys);
#elif defined(FFT_FFTW2)
  if (flag == -1)
    fftw(plan->plan_slow_forward,total/length,data,1,length,NULL,0,0);
  else
    fftw(plan->plan_slow_backward,total/length,data,1,length,NULL,0,0);
#elif defined(FFT_FFTW3)
  if (flag == -1)
    FFTW_API(execute_dft)(plan->plan_slow_forward,data,data);
  else
    FFTW_API(execute_dft)(plan->plan_slow_backward,data,data);
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
    remap_3d((double *) data, (double *) out, (double *) plan->scratch,
	     plan->post_plan);

  // scaling if required 

#ifndef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->normnum;
    for (i = 0; i < num; i++) {
#ifdef FFT_FFTW3
      out[i][0] *= norm;
      out[i][1] *= norm;
#else
      out[i].re *= norm;
      out[i].im *= norm;
#endif
    }
  }
#endif

#ifdef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = plan->normnum;
    for (i = 0; i < num; i++) out[i] *= (norm,norm);
  }
#endif
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
------------------------------------------------------------------------- */

struct fft_plan_3d *fft_3d_create_plan(
       MPI_Comm comm, int nfast, int nmid, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int in_klo, int in_khi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int out_klo, int out_khi,
       int scaled, int permute, int *nbuf)
{
  struct fft_plan_3d *plan;
  int me,nprocs;
  int i,num,flag,remapflag,fftflag;
  int first_ilo,first_ihi,first_jlo,first_jhi,first_klo,first_khi;
  int second_ilo,second_ihi,second_jlo,second_jhi,second_klo,second_khi;
  int third_ilo,third_ihi,third_jlo,third_jhi,third_klo,third_khi;
  int out_size,first_size,second_size,third_size,copy_size,scratch_size;
  int np1,np2,ip1,ip2;
  int list[50];

  // system specific variables 

#ifdef FFT_SCSL
  FFT_DATA dummy_d[5];
  FFT_PREC dummy_p[5];
  int isign,isys;
  FFT_PREC scalef;
#endif
#ifdef FFT_INTEL
  FFT_DATA dummy;
#endif
#ifdef FFT_T3E
  FFT_DATA dummy[5];
  int isign,isys;
  double scalef;
#endif

  // query MPI info 

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

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
  }
  else {
    first_ilo = 0;
    first_ihi = nfast - 1;
    first_jlo = ip1*nmid/np1;
    first_jhi = (ip1+1)*nmid/np1 - 1;
    first_klo = ip2*nslow/np2;
    first_khi = (ip2+1)*nslow/np2 - 1;
    plan->pre_plan =
      remap_3d_create_plan(comm,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
			   first_ilo,first_ihi,first_jlo,first_jhi,
			   first_klo,first_khi,2,0,0,FFT_PRECISION);
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
			   second_klo,second_khi,2,1,0,FFT_PRECISION);
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
  }
  else {
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
			 third_ilo,third_ihi,2,1,0,FFT_PRECISION);
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
			   out_jlo,out_jhi,2,(permute+1)%3,0,FFT_PRECISION);
    if (plan->post_plan == NULL) return NULL;
  }

  // configure plan memory pointers and allocate work space
  // out_size = amount of memory given to FFT by user
  // first/second/third_size = amount of memory needed after pre,mid1,mid2 remaps
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

#if defined(FFT_SGI)

  plan->coeff1 = (FFT_DATA *) malloc((nfast+15)*sizeof(FFT_DATA));
  plan->coeff2 = (FFT_DATA *) malloc((nmid+15)*sizeof(FFT_DATA));
  plan->coeff3 = (FFT_DATA *) malloc((nslow+15)*sizeof(FFT_DATA));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL ||
      plan->coeff3 == NULL) return NULL;

  FFT_1D_INIT(nfast,plan->coeff1);
  FFT_1D_INIT(nmid,plan->coeff2);
  FFT_1D_INIT(nslow,plan->coeff3);

  if (scaled == 0) 
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

#elif defined(FFT_SCSL)

  plan->coeff1 = (FFT_PREC *) malloc((2*nfast+30)*sizeof(FFT_PREC));
  plan->coeff2 = (FFT_PREC *) malloc((2*nmid+30)*sizeof(FFT_PREC));
  plan->coeff3 = (FFT_PREC *) malloc((2*nslow+30)*sizeof(FFT_PREC));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL || 
      plan->coeff3 == NULL) return NULL;

  plan->work1 = (FFT_PREC *) malloc((2*nfast)*sizeof(FFT_PREC));
  plan->work2 = (FFT_PREC *) malloc((2*nmid)*sizeof(FFT_PREC));
  plan->work3 = (FFT_PREC *) malloc((2*nslow)*sizeof(FFT_PREC));

  if (plan->work1 == NULL || plan->work2 == NULL || 
      plan->work3 == NULL) return NULL;

  isign = 0;
  scalef = 1.0;
  isys = 0;

  FFT_1D_INIT(isign,nfast,scalef,dummy_d,dummy_d,plan->coeff1,dummy_p,&isys);
  FFT_1D_INIT(isign,nmid,scalef,dummy_d,dummy_d,plan->coeff2,dummy_p,&isys);
  FFT_1D_INIT(isign,nslow,scalef,dummy_d,dummy_d,plan->coeff3,dummy_p,&isys);

  if (scaled == 0) 
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

#elif defined(FFT_INTEL)

  flag = 0;

  num = 0;
  factor(nfast,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;
  num = 0;
  factor(nmid,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;
  num = 0;
  factor(nslow,&num,list);
  for (i = 0; i < num; i++)
    if (list[i] != 2 && list[i] != 3 && list[i] != 5) flag = 1;

  MPI_Allreduce(&flag,&fftflag,1,MPI_INT,MPI_MAX,comm);
  if (fftflag) {
    if (me == 0) printf("ERROR: FFTs are not power of 2,3,5\n");
    return NULL;
  }

  plan->coeff1 = (FFT_DATA *) malloc((3*nfast/2+1)*sizeof(FFT_DATA));
  plan->coeff2 = (FFT_DATA *) malloc((3*nmid/2+1)*sizeof(FFT_DATA));
  plan->coeff3 = (FFT_DATA *) malloc((3*nslow/2+1)*sizeof(FFT_DATA));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL || 
      plan->coeff3 == NULL) return NULL;

  flag = 0;
  FFT_1D_INIT(&dummy,&nfast,&flag,plan->coeff1);
  FFT_1D_INIT(&dummy,&nmid,&flag,plan->coeff2);
  FFT_1D_INIT(&dummy,&nslow,&flag,plan->coeff3);

  if (scaled == 0) {
    plan->scaled = 1;
    plan->norm = nfast*nmid*nslow;
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }
  else
    plan->scaled = 0;

#elif defined(FFT_DEC)

  if (scaled == 0) {
    plan->scaled = 1;
    plan->norm = nfast*nmid*nslow;
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }
  else
    plan->scaled = 0;

#elif defined(FFT_T3E)

  plan->coeff1 = (double *) malloc((12*nfast)*sizeof(double));
  plan->coeff2 = (double *) malloc((12*nmid)*sizeof(double));
  plan->coeff3 = (double *) malloc((12*nslow)*sizeof(double));

  if (plan->coeff1 == NULL || plan->coeff2 == NULL || 
      plan->coeff3 == NULL) return NULL;

  plan->work1 = (double *) malloc((8*nfast)*sizeof(double));
  plan->work2 = (double *) malloc((8*nmid)*sizeof(double));
  plan->work3 = (double *) malloc((8*nslow)*sizeof(double));

  if (plan->work1 == NULL || plan->work2 == NULL || 
      plan->work3 == NULL) return NULL;

  isign = 0;
  scalef = 1.0;
  isys = 0;

  FFT_1D_INIT(&isign,&nfast,&scalef,dummy,dummy,plan->coeff1,dummy,&isys);
  FFT_1D_INIT(&isign,&nmid,&scalef,dummy,dummy,plan->coeff2,dummy,&isys);
  FFT_1D_INIT(&isign,&nslow,&scalef,dummy,dummy,plan->coeff3,dummy,&isys);

  if (scaled == 0) 
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

#elif defined(FFT_FFTW2)

  plan->plan_fast_forward = 
    fftw_create_plan(nfast,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  plan->plan_fast_backward = 
    fftw_create_plan(nfast,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);

  if (nmid == nfast) {
    plan->plan_mid_forward = plan->plan_fast_forward;
    plan->plan_mid_backward = plan->plan_fast_backward;
  }
  else {
    plan->plan_mid_forward = 
      fftw_create_plan(nmid,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
    plan->plan_mid_backward = 
      fftw_create_plan(nmid,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  }

  if (nslow == nfast) {
    plan->plan_slow_forward = plan->plan_fast_forward;
    plan->plan_slow_backward = plan->plan_fast_backward;
  }
  else if (nslow == nmid) {
    plan->plan_slow_forward = plan->plan_mid_forward;
    plan->plan_slow_backward = plan->plan_mid_backward;
  }
  else {
    plan->plan_slow_forward = 
      fftw_create_plan(nslow,FFTW_FORWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
    plan->plan_slow_backward = 
      fftw_create_plan(nslow,FFTW_BACKWARD,FFTW_ESTIMATE | FFTW_IN_PLACE);
  }

  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

#elif defined(FFT_FFTW3)
  int maxtotal = MAX(MAX(plan->total1,plan->total2),plan->total3);
  plan->dummy_data = (FFT_DATA *) malloc(maxtotal*sizeof(FFT_DATA));
  
  plan->plan_fast_forward = 
    FFTW_API(plan_many_dft)(1,&nfast,plan->total1/nfast,
			    plan->dummy_data,NULL,1,nfast,
			    plan->dummy_data,NULL,1,nfast,
			    FFTW_FORWARD,FFTW_PATIENT);

  plan->plan_fast_backward = 
    FFTW_API(plan_many_dft)(1,&nfast,plan->total1/nfast,
			    plan->dummy_data,NULL,1,nfast,
			    plan->dummy_data,NULL,1,nfast,
			    FFTW_BACKWARD,FFTW_PATIENT);

  if (nmid == nfast) {
    plan->plan_mid_forward = plan->plan_fast_forward;
    plan->plan_mid_backward = plan->plan_fast_backward;
  }
  else {
    plan->plan_mid_forward = 
      FFTW_API(plan_many_dft)(1,&nmid,plan->total2/nmid,
			      plan->dummy_data,NULL,1,nmid,
			      plan->dummy_data,NULL,1,nmid,
			      FFTW_FORWARD,FFTW_PATIENT);

    plan->plan_mid_backward = 
      FFTW_API(plan_many_dft)(1,&nmid,plan->total2/nmid,
			      plan->dummy_data,NULL,1,nmid,
			      plan->dummy_data,NULL,1,nmid,
			      FFTW_BACKWARD,FFTW_PATIENT);
  }

  if (nslow == nfast) {
    plan->plan_slow_forward = plan->plan_fast_forward;
    plan->plan_slow_backward = plan->plan_fast_backward;
  }
  else if (nslow == nmid) {
    plan->plan_slow_forward = plan->plan_mid_forward;
    plan->plan_slow_backward = plan->plan_mid_backward;
  }
  else {
    plan->plan_slow_forward = 
      FFTW_API(plan_many_dft)(1,&nslow,plan->total3/nslow,
			      plan->dummy_data,NULL,1,nslow,
			      plan->dummy_data,NULL,1,nslow,
			      FFTW_FORWARD,FFTW_PATIENT);

    plan->plan_slow_backward = 
      FFTW_API(plan_many_dft)(1,&nslow,plan->total3/nslow,
			      plan->dummy_data,NULL,1,nslow,
			      plan->dummy_data,NULL,1,nslow,
			      FFTW_BACKWARD,FFTW_PATIENT);
  }

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

#if defined(FFT_SGI)
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->coeff3);
#elif defined(FFT_SCSL)
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->coeff3);
  free(plan->work1);
  free(plan->work2);
  free(plan->work3);
#elif defined(FFT_INTEL)
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->coeff3);
#elif defined(FFT_T3E)
  free(plan->coeff1);
  free(plan->coeff2);
  free(plan->coeff3);
  free(plan->work1);
  free(plan->work2);
  free(plan->work3);
#elif defined(FFT_FFTW2)
  if (plan->plan_slow_forward != plan->plan_fast_forward &&
      plan->plan_slow_forward != plan->plan_mid_forward) {
    fftw_destroy_plan(plan->plan_slow_forward);
    fftw_destroy_plan(plan->plan_slow_backward);
  }
  if (plan->plan_mid_forward != plan->plan_fast_forward) {
    fftw_destroy_plan(plan->plan_mid_forward);
    fftw_destroy_plan(plan->plan_mid_backward);
  }
  fftw_destroy_plan(plan->plan_fast_forward);
  fftw_destroy_plan(plan->plan_fast_backward);
#elif defined(FFT_FFTW3)
  if (plan->plan_slow_forward != plan->plan_fast_forward &&
      plan->plan_slow_forward != plan->plan_mid_forward) {
    FFTW_API(destroy_plan)(plan->plan_slow_forward);
    FFTW_API(destroy_plan)(plan->plan_slow_backward);
  }
  if (plan->plan_mid_forward != plan->plan_fast_forward) {
    FFTW_API(destroy_plan)(plan->plan_mid_forward);
    FFTW_API(destroy_plan)(plan->plan_mid_backward);
  }
  FFTW_API(destroy_plan)(plan->plan_fast_forward);
  FFTW_API(destroy_plan)(plan->plan_fast_backward);
  free(plan->dummy_data);
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
  }
  else if (n % 2 == 0) {
    *list = 2;
    (*num)++;
    factor(n/2,num,list+1);
  }
  else if (n % 3 == 0) {
    *list = 3;
    (*num)++;
    factor(n/3,num,list+1);
  }
  else if (n % 5 == 0) {
    *list = 5;
    (*num)++;
    factor(n/5,num,list+1);
  }
  else if (n % 7 == 0) {
    *list = 7;
    (*num)++;
    factor(n/7,num,list+1);
  }
  else if (n % 11 == 0) {
    *list = 11;
    (*num)++;
    factor(n/11,num,list+1);
  }
  else if (n % 13 == 0) {
    *list = 13;
    (*num)++;
    factor(n/13,num,list+1);
  }
  else {
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
  int i,total,length,offset,num;
  double norm;

  // system specific constants 

#ifdef FFT_SCSL
  int isys = 0;
  FFT_PREC scalef = 1.0;
#endif
#ifdef FFT_DEC
  char c = 'C';
  char f = 'F';
  char b = 'B';
  int one = 1;
#endif
#ifdef FFT_T3E
  int isys = 0;
  double scalef = 1.0;
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

  if (total1 > nsize) total1 = (nsize/length1) * length1;
  if (total2 > nsize) total2 = (nsize/length2) * length2;
  if (total3 > nsize) total3 = (nsize/length3) * length3;

  // perform 1d FFTs in each of 3 dimensions
  // data is just an array of 0.0

#ifdef FFT_SGI
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(flag,length1,&data[offset],1,plan->coeff1);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(flag,length2,&data[offset],1,plan->coeff2);
  for (offset = 0; offset < total3; offset += length3)
    FFT_1D(flag,length3,&data[offset],1,plan->coeff3);
#elif defined(FFT_SCSL)
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(flag,length1,scalef,&data[offset],&data[offset],plan->coeff1,
	   plan->work1,&isys);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(flag,length2,scalef,&data[offset],&data[offset],plan->coeff2,
	   plan->work2,&isys);
  for (offset = 0; offset < total3; offset += length3)
    FFT_1D(flag,length3,scalef,&data[offset],&data[offset],plan->coeff3,
	   plan->work3,&isys);
#elif defined(FFT_INTEL)
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(&data[offset],&length1,&flag,plan->coeff1);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(&data[offset],&length2,&flag,plan->coeff2);
  for (offset = 0; offset < total3; offset += length3)
    FFT_1D(&data[offset],&length3,&flag,plan->coeff3);
#elif defined(FFT_DEC)
  if (flag == -1) {
    for (offset = 0; offset < total1; offset += length1)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length1,&one);
    for (offset = 0; offset < total2; offset += length2)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length2,&one);
    for (offset = 0; offset < total3; offset += length3)
      FFT_1D(&c,&c,&f,&data[offset],&data[offset],&length3,&one);
  } else {
    for (offset = 0; offset < total1; offset += length1)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length1,&one);
    for (offset = 0; offset < total2; offset += length2)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length2,&one);
    for (offset = 0; offset < total3; offset += length3)
      FFT_1D(&c,&c,&b,&data[offset],&data[offset],&length3,&one);
  }
#elif defined(FFT_T3E)
  for (offset = 0; offset < total1; offset += length1)
    FFT_1D(&flag,&length1,&scalef,&data[offset],&data[offset],plan->coeff1,
	   plan->work1,&isys);
  for (offset = 0; offset < total2; offset += length2)
    FFT_1D(&flag,&length2,&scalef,&data[offset],&data[offset],plan->coeff2,
	   plan->work2,&isys);
  for (offset = 0; offset < total3; offset += length3)
    FFT_1D(&flag,&length3,&scalef,&data[offset],&data[offset],plan->coeff3,
	   plan->work3,&isys);
#elif defined(FFT_FFTW2)
  if (flag == -1) {
    fftw(plan->plan_fast_forward,total1/length1,data,1,0,NULL,0,0);
    fftw(plan->plan_mid_forward,total2/length2,data,1,0,NULL,0,0);
    fftw(plan->plan_slow_forward,total3/length3,data,1,0,NULL,0,0);
  } else {
    fftw(plan->plan_fast_backward,total1/length1,data,1,0,NULL,0,0);
    fftw(plan->plan_mid_backward,total2/length2,data,1,0,NULL,0,0);
    fftw(plan->plan_slow_backward,total3/length3,data,1,0,NULL,0,0);
  }
#elif defined(FFT_FFTW3)
  if (flag == -1) {
    FFTW_API(execute_dft)(plan->plan_fast_forward,data,data);
    FFTW_API(execute_dft)(plan->plan_mid_forward,data,data);
    FFTW_API(execute_dft)(plan->plan_slow_forward,data,data);
  } else {
    FFTW_API(execute_dft)(plan->plan_fast_backward,data,data);
    FFTW_API(execute_dft)(plan->plan_mid_backward,data,data);
    FFTW_API(execute_dft)(plan->plan_slow_backward,data,data);
  }
#else
  if (flag == -1) {
    for (offset = 0; offset < total1; offset += length1)
      kiss_fft(plan->cfg_fast_forward,&data[offset],&data[offset]);
    for (offset = 0; offset < total2; offset += length2)
      kiss_fft(plan->cfg_mid_forward,&data[offset],&data[offset]);
    for (offset = 0; offset < total3; offset += length3)
      kiss_fft(plan->cfg_slow_forward,&data[offset],&data[offset]);
  } else {
    for (offset = 0; offset < total1; offset += length1)
      kiss_fft(plan->cfg_fast_backward,&data[offset],&data[offset]);
    for (offset = 0; offset < total2; offset += length2)
      kiss_fft(plan->cfg_mid_backward,&data[offset],&data[offset]);
    for (offset = 0; offset < total3; offset += length3)
      kiss_fft(plan->cfg_slow_backward,&data[offset],&data[offset]);
  }
#endif

  // scaling if required 
  // limit num to size of data

#ifndef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = MIN(plan->normnum,nsize);
    for (i = 0; i < num; i++) {
#ifdef FFT_FFTW3
      data[i][0] *= norm;
      data[i][1] *= norm;
#else
      data[i].re *= norm;
      data[i].im *= norm;
#endif
    }
  }
#endif

#ifdef FFT_T3E
  if (flag == 1 && plan->scaled) {
    norm = plan->norm;
    num = MIN(plan->normnum,nsize);
    for (i = 0; i < num; i++) data[i] *= (norm,norm);
  }
#endif
}

#ifdef FFT_KISSFFT
/*
  Copyright (c) 2003-2008, Mark Borgerding

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
   this code is adapted from kiss_fft_v1_2_8 
   homepage: http://kissfft.sf.net/ 
   
   changes 2008-2010 by Axel Kohlmeyer <akohlmey@gmail.com>
*/
#if defined(_OPENMP)
#include <omp.h>
#endif
 
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944
#endif

/*
  Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
 * */

#define S_MUL(a,b) ( (a)*(b) )
#define C_MUL(m,a,b) \
    do{ (m).re = (a).re*(b).re - (a).im*(b).im;\
        (m).im = (a).re*(b).im + (a).im*(b).re; }while(0)
#define C_FIXDIV(c,div) /* NOOP */
#define C_MULBYSCALAR( c, s ) \
    do{ (c).re *= (s);\
        (c).im *= (s); }while(0)

#ifndef CHECK_OVERFLOW_OP
#  define CHECK_OVERFLOW_OP(a,op,b) /* noop */
#endif

#define  C_ADD( res, a,b)\
    do { \
	    CHECK_OVERFLOW_OP((a).re,+,(b).re)\
	    CHECK_OVERFLOW_OP((a).im,+,(b).im)\
	    (res).re=(a).re+(b).re;  (res).im=(a).im+(b).im; \
    }while(0)
#define  C_SUB( res, a,b)\
    do { \
	    CHECK_OVERFLOW_OP((a).re,-,(b).re)\
	    CHECK_OVERFLOW_OP((a).im,-,(b).im)\
	    (res).re=(a).re-(b).re;  (res).im=(a).im-(b).im; \
    }while(0)
#define C_ADDTO( res , a)\
    do { \
	    CHECK_OVERFLOW_OP((res).re,+,(a).re)\
	    CHECK_OVERFLOW_OP((res).im,+,(a).im)\
	    (res).re += (a).re;  (res).im += (a).im;\
    }while(0)

#define C_SUBFROM( res , a)\
    do {\
	    CHECK_OVERFLOW_OP((res).re,-,(a).re)\
	    CHECK_OVERFLOW_OP((res).im,-,(a).im)\
	    (res).re -= (a).re;  (res).im -= (a).im; \
    }while(0)


#define KISS_FFT_COS(phase) (kiss_fft_scalar) cos(phase)
#define KISS_FFT_SIN(phase) (kiss_fft_scalar) sin(phase)
#define HALF_OF(x) ((x)*.5)

#define  kf_cexp(x,phase) \
	do{ \
		(x)->re = KISS_FFT_COS(phase);\
		(x)->im = KISS_FFT_SIN(phase);\
	}while(0)


/* a debugging function */
#define pcpx(c)\
    fprintf(stderr,"%g + %gi\n",(double)((c)->re),(double)((c)->im) )

/* The guts header contains all the multiplication and addition macros that are defined for
 fixed or floating point complex numbers.  It also delares the kf_ internal functions.
 */

static FFT_DATA *scratchbuf=NULL;
static size_t nscratchbuf=0;
static FFT_DATA *tmpbuf=NULL;
static size_t ntmpbuf=0;

#define CHECKBUF(buf,nbuf,n) \
    do { \
        if ( nbuf < (size_t)(n) ) {\
            free(buf); \
            buf = (FFT_DATA*)KISS_FFT_MALLOC(sizeof(FFT_DATA)*(n)); \
            nbuf = (size_t)(n); \
        } \
   }while(0)


static void kf_bfly2(FFT_DATA *Fout, const size_t fstride,
                     const kiss_fft_cfg st, int m)
{
    FFT_DATA *Fout2;
    FFT_DATA *tw1 = st->twiddles;
    FFT_DATA t;

    Fout2 = Fout + m;
    do {
        C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
        ++Fout2;
        ++Fout;
    } while(--m);
}

static void kf_bfly4(FFT_DATA * Fout, const size_t fstride,
                     const kiss_fft_cfg st, const size_t m)
{
    FFT_DATA *tw1, *tw2, *tw3;
    FFT_DATA scratch[6];
    size_t k=m;
    const size_t m2=2*m;
    const size_t m3=3*m;

    tw3 = tw2 = tw1 = st->twiddles;

    do {
        C_FIXDIV(*Fout,4); C_FIXDIV(Fout[m],4); C_FIXDIV(Fout[m2],4); C_FIXDIV(Fout[m3],4);

        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );

        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if (st->inverse) {
            Fout[m].re = scratch[5].re - scratch[4].im;
            Fout[m].im = scratch[5].im + scratch[4].re;
            Fout[m3].re = scratch[5].re + scratch[4].im;
            Fout[m3].im = scratch[5].im - scratch[4].re;
        } else{
            Fout[m].re = scratch[5].re + scratch[4].im;
            Fout[m].im = scratch[5].im - scratch[4].re;
            Fout[m3].re = scratch[5].re - scratch[4].im;
            Fout[m3].im = scratch[5].im + scratch[4].re;
        }
        ++Fout;
    } while(--k);
}

static void kf_bfly3(FFT_DATA * Fout, const size_t fstride,
                     const kiss_fft_cfg st, size_t m)
{
    size_t k=m;
    const size_t m2 = 2*m;
    FFT_DATA *tw1, *tw2;
    FFT_DATA scratch[5];
    FFT_DATA epi3;
    epi3 = st->twiddles[fstride*m];

    tw1=tw2=st->twiddles;

    do {
        C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

        C_MUL(scratch[1],Fout[m] , *tw1);
        C_MUL(scratch[2],Fout[m2] , *tw2);

        C_ADD(scratch[3],scratch[1],scratch[2]);
        C_SUB(scratch[0],scratch[1],scratch[2]);
        tw1 += fstride;
        tw2 += fstride*2;

        Fout[m].re = Fout->re - HALF_OF(scratch[3].re);
        Fout[m].im = Fout->im - HALF_OF(scratch[3].im);

        C_MULBYSCALAR( scratch[0] , epi3.im );

        C_ADDTO(*Fout,scratch[3]);

        Fout[m2].re = Fout[m].re + scratch[0].im;
        Fout[m2].im = Fout[m].im - scratch[0].re;

        Fout[m].re -= scratch[0].im;
        Fout[m].im += scratch[0].re;

        ++Fout;
    } while(--k);
}

static void kf_bfly5(FFT_DATA * Fout, const size_t fstride,
                     const kiss_fft_cfg st, int m)
{
    FFT_DATA *Fout0, *Fout1, *Fout2, *Fout3, *Fout4;
    int u;
    FFT_DATA scratch[13];
    FFT_DATA * twiddles = st->twiddles;
    FFT_DATA *tw;
    FFT_DATA ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

        Fout0->re += scratch[7].re + scratch[8].re;
        Fout0->im += scratch[7].im + scratch[8].im;

        scratch[5].re = scratch[0].re + S_MUL(scratch[7].re,ya.re) + S_MUL(scratch[8].re,yb.re);
        scratch[5].im = scratch[0].im + S_MUL(scratch[7].im,ya.re) + S_MUL(scratch[8].im,yb.re);

        scratch[6].re =  S_MUL(scratch[10].im,ya.im) + S_MUL(scratch[9].im,yb.im);
        scratch[6].im = -S_MUL(scratch[10].re,ya.im) - S_MUL(scratch[9].re,yb.im);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);

        scratch[11].re = scratch[0].re + S_MUL(scratch[7].re,yb.re) + S_MUL(scratch[8].re,ya.re);
        scratch[11].im = scratch[0].im + S_MUL(scratch[7].im,yb.re) + S_MUL(scratch[8].im,ya.re);
        scratch[12].re = - S_MUL(scratch[10].im,yb.im) + S_MUL(scratch[9].im,ya.im);
        scratch[12].im = S_MUL(scratch[10].re,yb.im) - S_MUL(scratch[9].re,ya.im);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(FFT_DATA * Fout, const size_t fstride,
                            const kiss_fft_cfg st, int m, int p)
{
    int u,k,q1,q;
    FFT_DATA * twiddles = st->twiddles;
    FFT_DATA t;
    int Norig = st->nfft;

    CHECKBUF(scratchbuf,nscratchbuf,p);

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratchbuf[q1] = Fout[ k  ];
            C_FIXDIV(scratchbuf[q1],p);
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratchbuf[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratchbuf[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
}

static void kf_work(FFT_DATA * Fout, const FFT_DATA *f,
                    const size_t fstride, int in_stride,
                    int * factors, const kiss_fft_cfg st)
{
    FFT_DATA * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const FFT_DATA * Fout_end = Fout + p*m;

/* add OpenMP support here */
#if defined(_OPENMP) && 0
    /* use openmp extensions at the top-level (not recursive).
     * NOTE: we can only use this strategy if the size of the 
     * data set cannot be factorized. <g>.
     */
    if (fstride==1 && m != 1) {
        int k;

/* execute the p different work units in different threads */
#pragma omp parallel for private(k) schedule(static)
        for (k=0;k<p;++k) 
            kf_work( Fout +k*m, f+ fstride*in_stride*k,fstride*p,in_stride,factors,st);
        /* all threads have joined by this point */

        switch (p) {
            case 2: kf_bfly2(Fout,fstride,st,m); break;
            case 3: kf_bfly3(Fout,fstride,st,m); break; 
            case 4: kf_bfly4(Fout,fstride,st,m); break;
            case 5: kf_bfly5(Fout,fstride,st,m); break; 
            default: kf_bfly_generic(Fout,fstride,st,m,p); break;
        }
        return;
    }
#endif

    if (m==1) {
        do {
            *Fout = *f;
            f += fstride*in_stride;
        } while (++Fout != Fout_end);
    } else {
        do {
            /* recursive call:
               DFT of size m*p performed by doing
               p instances of smaller DFTs of size m, 
               each one takes a decimated version of the input */
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        } while( (Fout += m) != Fout_end);
    }

    Fout=Fout_beg;

    /* recombine the p smaller DFTs */
    switch (p) {
      case 2: kf_bfly2(Fout,fstride,st,m); break;
      case 3: kf_bfly3(Fout,fstride,st,m); break; 
      case 4: kf_bfly4(Fout,fstride,st,m); break;
      case 5: kf_bfly5(Fout,fstride,st,m); break; 
      default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where 
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static void kf_factor(int n, int *facbuf)
{
    int p=4, nf=0;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

#if defined(_OPENMP) && 0
/* to maximally utilize the available cpus we try to match the
   first radix to the number of threads as closely as possible. */
    p=omp_get_num_threads();
    while (n % p) --p;
    if (p > 1) {
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
        ++nf;
    }
    p=4; /* back to normal */
#endif

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
kiss_fft_cfg kiss_fft_alloc(int nfft, int inverse_fft, void *mem, size_t *lenmem)
{
    kiss_fft_cfg st=NULL;
    size_t memneeded = sizeof(struct kiss_fft_state)
        + sizeof(FFT_DATA)*(nfft-1); /* twiddle factors */

    if (lenmem==NULL) {
        st=(kiss_fft_cfg)KISS_FFT_MALLOC( memneeded );
    } else {
        if (mem != NULL && *lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }

    if (st) {
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;

        for (i=0;i<nfft;++i) {
            const double phase = (st->inverse ? 2.0*M_PI:-2.0*M_PI)*i / nfft;
            kf_cexp(st->twiddles+i, phase );
        }

        kf_factor(nfft,st->factors);
    }
    return st;
}
   
void kiss_fft_stride(kiss_fft_cfg st, const FFT_DATA *fin, FFT_DATA *fout, int in_stride)
{
    if (fin == fout) {
        CHECKBUF(tmpbuf,ntmpbuf,st->nfft);
        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(FFT_DATA)*st->nfft);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

void kiss_fft(kiss_fft_cfg cfg, const FFT_DATA *fin, FFT_DATA *fout)
{
    kiss_fft_stride(cfg,fin,fout,1);
}


/* not really necessary to call, but if someone is doing in-place ffts,
   they may want to free the buffers from CHECKBUF
 */ 
void kiss_fft_cleanup(void)
{
    free(scratchbuf);
    scratchbuf = NULL;
    nscratchbuf=0;
    free(tmpbuf);
    tmpbuf=NULL;
    ntmpbuf=0;
}

int kiss_fft_next_fast_size(int n)
{
    while(1) {
        int m=n;
        while ( (m%2) == 0 ) m/=2;
        while ( (m%3) == 0 ) m/=3;
        while ( (m%5) == 0 ) m/=5;
        if (m<=1)
            break; /* n is completely factorable by twos, threes, and fives */
        n++;
    }
    return n;
}

#endif
