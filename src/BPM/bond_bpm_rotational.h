/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(bpm/rotational,BondBPMRotational)
// clang-format on
#else

#ifndef LMP_BOND_BPM_ROTATIONAL_H
#define LMP_BOND_BPM_ROTATIONAL_H

#include "bond_bpm.h"

namespace LAMMPS_NS {

class BondBPMRotational : public BondBPM {
 public:
  BondBPMRotational(class LAMMPS *);
  virtual ~BondBPMRotational();
  virtual void compute(int, int);
  void coeff(int, char **);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);

 protected:
  double *Kr, *Ks, *Kt, *Kb, *gnorm, *gslide, *groll, *gtwist;
  double *Fcr, *Fcs, *Tct, *Tcb;
  double acos_limit(double);

  void allocate();
  void store_data();
  double store_bond(int, int, int);
  class FixBondHistory *fix_bond_history;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Atom missing in BPM bond

Bonded atom cannot be found

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

E: Bond bpm/rotational requires atom style sphere/bpm

Self-explanatory.

E: Bond style bpm requires 1-3 and 1-4 special weights of 1.0

Self-explanatory.

W: Bond style bpm/rotational not intended for 2d use, may be inefficient

This bond style will perform a lot of unnecessary calculations in 2d

*/
