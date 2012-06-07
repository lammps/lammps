/* -*- c++ -*- ----------------------------------------------------------
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

#ifndef FFT3D_WRAP_CUDA_H_
#define FFT3D_WRAP_CUDA_H_

#include "pointers.h"

#ifdef FFT_CUFFT
  #include "fft3d_cuda.h"
#endif
#ifndef FFT_CUFFT
  #include "fft3d.h"
#endif

namespace LAMMPS_NS {

class FFT3dCuda : protected Pointers {
 public:
  FFT3dCuda(class LAMMPS *, MPI_Comm,int,int,int,int,int,int,int,int,int,
        int,int,int,int,int,int,int,int,int *,bool);
  ~FFT3dCuda();
  void compute(double *, double *, int);
  void timing1d(double *, int, int);

#ifdef FFT_CUFFT
  void set_cudata(void* cudata,void* cudata2);
#endif
 private:
  struct fft_plan_3d *plan;
};

}

#endif /*FFT3D_WRAP_CUDA_H_*/
