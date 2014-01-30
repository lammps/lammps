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

#ifdef PAIR_CLASS

PairStyle(lj/cut/dipole/long,PairLJCutDipoleLong)

#else

#ifndef LMP_PAIR_LJ_CUT_DIPOLE_LONG_H
#define LMP_PAIR_LJ_CUT_DIPOLE_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutDipoleLong : public Pair {
 public:
  double cut_coul;
  double **sigma;

  PairLJCutDipoleLong(class LAMMPS *);
  ~PairLJCutDipoleLong();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 private:
  double cut_lj_global;
  double **cut_lj,**cut_ljsq;
  double cut_coulsq;
  double **epsilon;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double g_ewald;
  int ewald_order;
  virtual void *extract(const char *, int &);

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_style command

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair dipole/long requires atom attributes q, mu, torque

The atom style defined does not have these attributes.

E: Cannot (yet) use 'electron' units with dipoles

This feature is not yet supported.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
