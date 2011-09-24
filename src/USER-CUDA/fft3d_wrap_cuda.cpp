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

#include "mpi.h"
#include "fft3d_wrap_cuda.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FFT3dCuda::FFT3dCuda(LAMMPS *lmp, MPI_Comm comm, int nfast, int nmid, int nslow,
	     int in_ilo, int in_ihi, int in_jlo, int in_jhi,
	     int in_klo, int in_khi,
	     int out_ilo, int out_ihi, int out_jlo, int out_jhi,
	     int out_klo, int out_khi,
	     int scaled, int permute, int *nbuf,bool init) : Pointers(lmp)
{
#ifdef FFT_CUFFT
  plan = fft_3d_create_plan_cuda(comm,nfast,nmid,nslow,
			    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
			    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
			    scaled,permute,nbuf,init);
#endif
#ifndef FFT_CUFFT
  plan = fft_3d_create_plan(comm,nfast,nmid,nslow,
			    in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
			    out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
			    scaled,permute,nbuf);
#endif
  if (plan == NULL) error->one(FLERR,"Could not create 3d FFT plan");
}

/* ---------------------------------------------------------------------- */

FFT3dCuda::~FFT3dCuda()
{
#ifdef FFT_CUFFT
  fft_3d_destroy_plan_cuda(plan);
#endif
#ifndef FFT_CUFFT
   fft_3d_destroy_plan(plan);
#endif
}

/* ---------------------------------------------------------------------- */

void FFT3dCuda::compute(double *in, double *out, int flag)
{
#ifdef FFT_CUFFT
  fft_3d_cuda((FFT_DATA *) in,(FFT_DATA *) out,flag,plan);
#endif
#ifndef FFT_CUFFT
  fft_3d((FFT_DATA *) in,(FFT_DATA *) out,flag,plan);
#endif
}

/* ---------------------------------------------------------------------- */

void FFT3dCuda::timing1d(double *in, int nsize, int flag)
{
#ifdef FFT_CUFFT
  fft_1d_only_cuda((FFT_DATA *) in,nsize,flag,plan);
#endif
#ifndef FFT_CUFFT
  fft_1d_only((FFT_DATA *) in,nsize,flag,plan);
#endif
}

#ifdef FFT_CUFFT
void FFT3dCuda::set_cudata(void* cudata,void* cudata2)
{ 
  
  plan->cudata=(cufftData*) cudata;
  plan->cudata2=(cufftData*) cudata2;
  
}
#endif
