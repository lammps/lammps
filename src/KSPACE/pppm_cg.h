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

KSpaceStyle(pppm/cg,PPPMCG)

#else

#ifndef LMP_PPPM_CG_H
#define LMP_PPPM_CG_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMCG : public PPPM {
 public:
  PPPMCG(class LAMMPS *, int, char **);
  virtual ~PPPMCG();
  virtual void compute(int eflag, int vflag);
  virtual double memory_usage();

 protected:
  int num_charged;
  int *is_charged;
  double smallq;

  virtual void particle_map();
  virtual void make_rho();
  virtual void fieldforce();
  virtual void slabcorr(int);
};

}

#endif
#endif
