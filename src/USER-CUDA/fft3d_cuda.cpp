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

/* ----------------------------------------------------------------------
   Contributing authors: Jim Shepherd (GA Tech) added SGI SCSL support
------------------------------------------------------------------------- */

#include "mpi.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "fft3d_cuda.h"
#include "fft3d_cuda_cu.h"
#include "remap.h"
#include <ctime>
#include "cuda_wrapper_cu.h"

#ifdef FFT_CUFFT
#endif

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

void fft_3d_cuda(FFT_DATA *in, FFT_DATA *out, int flag, struct fft_plan_3d *plan)
{
#ifdef FFT_CUFFT
  plan->iterate++;
  timespec starttime,starttime2;
  timespec endtime,endtime2;
	
  int i,total,length,offset,num;
  double norm;
  FFT_DATA *data,*copy;
  // system specific constants 


  // pre-remap to prepare for 1st FFTs if needed
  // copy = loc for remap result 
  int nprocs=plan->nprocs;
if(nprocs>1)
{
  if(plan->init)
  clock_gettime(CLOCK_REALTIME,&starttime);
  if (plan->pre_plan) {
    if (plan->pre_target == 0) copy = out;
    else copy = plan->copy;
    if(plan->init) remap_3d((double *) in, (double *) out, (double *) plan->scratch,plan->pre_plan);
    data = out;
  }
  else
    data = in;
}
  cufftResult retvalc;
  if(plan->init)
  {
	if(nprocs>1)
	{
      if(sizeof(FFT_FLOAT)==sizeof(double))cudaMemcpy((void*) (plan->cudata2), (void*) data, plan->cudatasize/2,cudaMemcpyHostToDevice);
      if(sizeof(FFT_FLOAT)==sizeof(float)) cudaMemcpy((void*) (plan->cudata2), (void*) data, plan->cudatasize,cudaMemcpyHostToDevice);
      initfftdata((double*)plan->cudata2,(FFT_FLOAT*)plan->cudata,plan->nfast,plan->nmid,plan->nslow);
    }
  }
    if (flag == -1)
    {
      retvalc=cufft(plan->plan_3d, plan->cudata, plan->cudata2,CUFFT_FORWARD);
    }
    else
    {
      retvalc=cufft(plan->plan_3d, plan->cudata, plan->cudata2,CUFFT_INVERSE);
    }
    if(retvalc!=CUFFT_SUCCESS) {printf("ErrorCUFFT: %i\n",retvalc);exit(EXIT_FAILURE);}

    FFTsyncthreads();
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

struct fft_plan_3d *fft_3d_create_plan_cuda(
       MPI_Comm comm, int nfast, int nmid, int nslow,
       int in_ilo, int in_ihi, int in_jlo, int in_jhi,
       int in_klo, int in_khi,
       int out_ilo, int out_ihi, int out_jlo, int out_jhi,
       int out_klo, int out_khi,
       int scaled, int permute, int *nbuf,bool ainit)
{
#ifdef FFT_CUFFT
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

  // query MPI info 

  MPI_Comm_rank(comm,&me);
  MPI_Comm_size(comm,&nprocs);

#ifndef FFT_CUFFT
    error->all(FLERR,"ERROR: Trying to use cuda fft without FFT_CUFFT set. Recompile with make option 'cufft=1'.");
#endif
  // compute division of procs in 2 dimensions not on-processor 
  bifactor_cuda(nprocs,&np1,&np2);
  ip1 = me % np1;
  ip2 = me/np1;

  // in case of CUDA FFT every proc does the full FFT in order to avoid data transfers (the problem is other wise heavily bandwidth limited)

  int ip1out = ip1;
  int ip2out = ip2;
  int np1out = np1;
  int np2out = np2;
  
  ip1 = 0;
  ip2 = 0;
  np1 = 1;
  np2 = 1;

  // allocate memory for plan data struct 

  plan = (struct fft_plan_3d *) malloc(sizeof(struct fft_plan_3d));
  if (plan == NULL) return NULL;
  plan->init=ainit;

  // remap from initial distribution to layout needed for 1st set of 1d FFTs
  // not needed if all procs own entire fast axis initially
  // first indices = distribution after 1st set of FFTs 

  if (in_ilo == 0 && in_ihi == nfast-1)
    flag = 0;
  else
    flag = 1;

  if(nprocs>1)flag=1;

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
    int members=2;
    if(plan->init) members=1;
    plan->pre_plan =
      remap_3d_create_plan(comm,in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
			   first_ilo,first_ihi,first_jlo,first_jhi,
			   first_klo,first_khi,
			   members,0,0,2);
    if (plan->pre_plan == NULL) return NULL;
  }

  // 1d FFTs along fast axis 

  plan->length1 = nfast;
  plan->total1 = nfast * nmid * nslow;

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
			   second_klo,second_khi,
			   2,1,0,2);
  if (plan->mid1_plan == NULL) return NULL;

  // 1d FFTs along mid axis 

  plan->length2 = nmid;
  plan->total2 = nfast * nmid * nslow;

  // remap from 2nd to 3rd FFT
  // if final distribution is permute=2 with all procs owning entire slow axis
  //   then this remapping goes directly to final distribution
  //  third indices = distribution after 3rd set of FFTs 

  flag=1;

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
			 third_ilo,third_ihi,
			 2,1,0,2);
  if (plan->mid2_plan == NULL) return NULL;

  // 1d FFTs along slow axis 

  plan->length3 = nslow;
  plan->total3 = nfast * nmid * nslow;

  // remap from 3rd FFT to final distribution
  //  not needed if permute = 2 and third indices = out indices on all procs 

  flag=1;

  MPI_Allreduce(&flag,&remapflag,1,MPI_INT,MPI_MAX,comm);

  if (remapflag == 0)
    plan->post_plan = NULL;
  else {
    plan->post_plan =
      remap_3d_create_plan(comm,
			   third_klo,third_khi,third_ilo,third_ihi,
			   third_jlo,third_jhi,
			   out_klo,out_khi,out_ilo,out_ihi,
			   out_jlo,out_jhi,
			   2,(permute+1)%3,0,2);
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

  plan->ihi_out=out_ihi;
  plan->ilo_out=out_ilo;
  plan->jhi_out=out_jhi;
  plan->jlo_out=out_jlo;
  plan->khi_out=out_khi;
  plan->klo_out=out_klo;

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

  cufftResult retvalc;
  int nfft = (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) *
    (in_khi-in_klo+1);
  int nfft_brick = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
    (out_khi-out_klo+1);
    
  int nfft_both = MAX(nfft,nfft_brick);
  nfft_both=nfast*nmid*nslow;

  plan->cudatasize=nfft_both*sizeof(FFT_DATA);

  //retvalc=cufftPlan1d(&(plan->plan_fast), nfast, CUFFT_PLAN,plan->total1/nfast);
  //if(retvalc!=CUFFT_SUCCESS) printf("ErrorCUFFT1: %i\n",retvalc);
  plan->nfast=nfast;

  //retvalc=cufftPlan1d(&(plan->plan_mid), nmid, CUFFT_PLAN,plan->total2/nmid);
  //if(retvalc!=CUFFT_SUCCESS) printf("ErrorCUFFT2: %i\n",retvalc);
  plan->nmid=nmid;

  //retvalc=cufftPlan1d(&(plan->plan_slow), nslow, CUFFT_PLAN,plan->total3/nslow);
  //if(retvalc!=CUFFT_SUCCESS) printf("ErrorCUFFT3: %i\n",retvalc);
  plan->nslow=nslow;

  retvalc=cufftPlan3d(&(plan->plan_3d), nslow,nmid,nfast, CUFFT_PLAN);
  if(retvalc!=CUFFT_SUCCESS) printf("ErrorCUFFT3: %i\n",retvalc);

  plan->nprocs=nprocs;
  plan->me=me;
  if (scaled == 0)
    plan->scaled = 0;
  else {
    plan->scaled = 1;
    plan->norm = 1.0/(nfast*nmid*nslow);
    plan->normnum = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
      (out_khi-out_klo+1);
  }

  plan->coretime=0;
  plan->iterate=0;
  plan->ffttime=0;
  return plan;
  #endif
}

/* ----------------------------------------------------------------------
   Destroy a 3d fft plan 
------------------------------------------------------------------------- */

void fft_3d_destroy_plan_cuda(struct fft_plan_3d *plan)
{
#ifdef FFT_CUFFT
  if (plan->pre_plan) remap_3d_destroy_plan(plan->pre_plan);
  if (plan->mid1_plan) remap_3d_destroy_plan(plan->mid1_plan);
  if (plan->mid2_plan) remap_3d_destroy_plan(plan->mid2_plan);
  if (plan->post_plan) remap_3d_destroy_plan(plan->post_plan);

  if (plan->copy) free(plan->copy);
  if (plan->scratch) free(plan->scratch);


  //cufftDestroy(plan->plan_fast);
  //cufftDestroy(plan->plan_mid);
  //cufftDestroy(plan->plan_slow);
  cufftDestroy(plan->plan_3d);
  free(plan);
#endif
}

/* ----------------------------------------------------------------------
   recursively divide n into small factors, return them in list
------------------------------------------------------------------------- */

void factor_cuda(int n, int *num, int *list)
{
  if (n == 1) {
    return;
  }
  else if (n % 2 == 0) {
    *list = 2;
    (*num)++;
    factor_cuda(n/2,num,list+1);
  }
  else if (n % 3 == 0) {
    *list = 3;
    (*num)++;
    factor_cuda(n/3,num,list+1);
  }
  else if (n % 5 == 0) {
    *list = 5;
    (*num)++;
    factor_cuda(n/5,num,list+1);
  }
  else if (n % 7 == 0) {
    *list = 7;
    (*num)++;
    factor_cuda(n/7,num,list+1);
  }
  else if (n % 11 == 0) {
    *list = 11;
    (*num)++;
    factor_cuda(n/11,num,list+1);
  }
  else if (n % 13 == 0) {
    *list = 13;
    (*num)++;
    factor_cuda(n/13,num,list+1);
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

void bifactor_cuda(int n, int *factor1, int *factor2)
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

void fft_1d_only_cuda(FFT_DATA *data, int nsize, int flag, struct fft_plan_3d *plan)
{
#ifdef FFT_CUFFT
  int i,total,length,offset,num;
  double norm;

  // system specific constants 



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


  cudaMemcpy((void**) &(plan->cudata), (void*) data, plan->cudatasize,cudaMemcpyHostToDevice);
  if (flag == -1) {
    cufft(plan->plan_3d, plan->cudata, plan->cudata,CUFFT_FORWARD);
    /*cufft(plan->plan_fast, plan->cudata, plan->cudata,CUFFT_FORWARD);
    cufft(plan->plan_mid, plan->cudata, plan->cudata,CUFFT_FORWARD);
    cufft(plan->plan_slow, plan->cudata, plan->cudata,CUFFT_FORWARD);*/
  } else {
    cufft(plan->plan_3d, plan->cudata, plan->cudata,CUFFT_FORWARD);
    /*cufft(plan->plan_fast, plan->cudata, plan->cudata,CUFFT_INVERSE);
    cufft(plan->plan_mid,plan->cudata, plan->cudata,CUFFT_INVERSE);
    cufft(plan->plan_slow, plan->cudata, plan->cudata,CUFFT_INVERSE);*/
  }
  cudaMemcpy((void*) data, (void**) &(plan->cudata), plan->cudatasize,cudaMemcpyDeviceToHost);

  // scaling if required 
  // limit num to size of data

#endif
}
