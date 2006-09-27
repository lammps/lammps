/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FFT3D_WRAP_H
#define FFT3D_WRAP_H

#include "lammps.h"
#include "fft3d.h"

class FFT3d : public LAMMPS {
 public:
  FFT3d(MPI_Comm,int,int,int,int,int,int,int,int,int,
	int,int,int,int,int,int,int,int,int *);
  ~FFT3d();
  void compute(double *, double *, int);
  void timing1d(double *, int, int);

 private:
  struct fft_plan_3d *plan;
};

#endif
