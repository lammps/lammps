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

#ifdef FIX_CLASS

FixStyle(restrain,FixRestrain)

#else

#ifndef LMP_FIX_RESTRAIN_H
#define LMP_FIX_RESTRAIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRestrain : public Fix {
 public:
  FixRestrain(class LAMMPS *, int, char **);
  ~FixRestrain();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();

 private:
  int nlevels_respa;
  int n_bonds, rstyle;
  double k_start, k_stop, energy, energy_all;
  int **atom_id;
  double *target, *cos_shift, *sin_shift;

  void restrain_dihedral();
};

}

#endif
#endif
