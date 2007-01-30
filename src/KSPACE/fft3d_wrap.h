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

#ifndef FFT3D_WRAP_H
#define FFT3D_WRAP_H

#include "pointers.h"
#include "fft3d.h"

namespace LAMMPS_NS {

class FFT3d : protected Pointers {
 public:
  FFT3d(class LAMMPS *, MPI_Comm,int,int,int,int,int,int,int,int,int,
	int,int,int,int,int,int,int,int,int *);
  ~FFT3d();
  void compute(double *, double *, int);
  void timing1d(double *, int, int);

 private:
  struct fft_plan_3d *plan;
};

}

#endif
