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

#ifndef FIX_NPT_ASPHERE_H
#define FIX_NPT_ASPHERE_H

#include "fix_npt.h"

namespace LAMMPS_NS {

class FixNPTAsphere : public FixNPT {
 public:
  FixNPTAsphere(class LAMMPS *, int, char **);
  ~FixNPTAsphere();
  void init();
  void initial_integrate(int);
  void final_integrate();

 private:
  double dtq;
  double factor_rotate;
  double **inertia;

  void richardson(double *, double *, double *);
  void omega_from_mq(double *, double *, double *, double *);
  void calculate_inertia();
};

}

#endif
