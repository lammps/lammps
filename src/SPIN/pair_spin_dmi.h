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

#ifdef PAIR_CLASS

PairStyle(spin/dmi,PairSpinDmi)

#else

#ifndef LMP_PAIR_SPIN_DMI_H
#define LMP_PAIR_SPIN_DMI_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinDmi : public PairSpin {
 public:
  PairSpinDmi(LAMMPS *lmp) : PairSpin(lmp) {}
  virtual ~PairSpinDmi();
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void *extract(const char *, int &);

  void compute(int, int);
  void compute_single_pair(int, double *);

  void compute_dmi(int, int, double *, double *, double *);
  void compute_dmi_mech(int, int, double, double *, double *, double *, double *);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  double cut_spin_dmi_global;                   // short range pair cutoff

 protected:
  double **DM;                                  // dmi coeff in eV
  double **v_dmx, **v_dmy, **v_dmz;             // dmi direction
  double **vmech_dmx, **vmech_dmy, **vmech_dmz; // dmi mech direction
  double **cut_spin_dmi;                        // cutoff distance dmi

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args in pair_spin command

Self-explanatory.

E: Spin simulations require metal unit style

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair spin requires atom attribute spin

The atom style defined does not have these attributes.

*/
