/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(spin/dipole/long,PairSpinDipoleLong)

#else

#ifndef LMP_PAIR_SPIN_DIPOLE_LONG_H
#define LMP_PAIR_SPIN_DIPOLE_LONG_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinDipoleLong : public PairSpin {
 public:
  double cut_coul;
  double **sigma;

  PairSpinDipoleLong(LAMMPS *);
  virtual ~PairSpinDipoleLong();
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void *extract(const char *, int &);

  void compute(int, int);
  void compute_single_pair(int, double *);

  void compute_long(int, int, double *, double *, double *,
      double *, double *);
  void compute_long_mech(int, int, double *, double *, double *,
      double *, double *);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  double cut_spin_long_global;  // global long cutoff distance

 protected:
  double hbar;                  // reduced Planck's constant
  double mub;                   // Bohr's magneton
  double mu_0;                  // vacuum permeability
  double mub2mu0;               // prefactor for mech force
  double mub2mu0hbinv;          // prefactor for mag force

  double **cut_spin_long;       // cutoff distance long

  double g_ewald;
  int ewald_order;

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

E: Can only use 'metal' units with spins

This feature is not yet supported.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
