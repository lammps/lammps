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

#ifdef ANGLE_CLASS

AngleStyle(cg/cmm/old,AngleCGCMM)

#else

#ifndef LMP_ANGLE_CG_CMM_H
#define LMP_ANGLE_CG_CMM_H

#include "stdio.h"
#include "angle.h"
#include "cg_cmm_parms.h"

namespace LAMMPS_NS {

class AngleCGCMM : public Angle, public CGCMMParms {
 public:
  AngleCGCMM(class LAMMPS *);
  ~AngleCGCMM();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 protected:
  void ev_tally_lj13(int, int, int, int, double, double,
                     double, double, double);

 private:
  double *k,*theta0;
  int *cg_type;
  double *epsilon, *sigma, *rcut;

  void allocate();
};

}

#endif
#endif
