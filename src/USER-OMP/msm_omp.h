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

#ifdef KSPACE_CLASS

KSpaceStyle(msm/omp,MSMOMP)

#else

#ifndef LMP_MSM_OMP_H
#define LMP_MSM_OMP_H

#include "msm.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

  class MSMOMP : public MSM, public ThrOMP {
 public:
  MSMOMP(class LAMMPS *);
  virtual ~MSMOMP () {};

 protected:
  virtual void direct(int);
  virtual void compute(int,int);

 private:
  template <int, int, int> void direct_eval(int);
  template <int> void direct_peratom(int);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Must use 'kspace_modify pressure/scalar no' with kspace_style msm/omp

The kspace scalar pressure option is not compatible with kspace_style msm/omp.

*/
