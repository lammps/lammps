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

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/gpu,PPPMGPU)

#else

#ifndef LMP_PPPM_GPU_H
#define LMP_PPPM_GPU_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMGPU : public PPPM {
 public:
  PPPMGPU(class LAMMPS *, int, char **);
  virtual ~PPPMGPU();
  virtual void init();
  virtual void compute(int, int);
  virtual double memory_usage();

 protected:

  FFT_SCALAR ***density_brick_gpu, ***vd_brick;

  virtual void allocate();
  virtual void deallocate();
  virtual void brick2fft();
  virtual void fillbrick();
  virtual void poisson(int, int);

  double poisson_time;  

  FFT_SCALAR ***create_3d_offset(int, int, int, int, int, int, const char *,
			     FFT_SCALAR *, int);
  void destroy_3d_offset(FFT_SCALAR ***, int, int);
};

}

#endif
#endif
