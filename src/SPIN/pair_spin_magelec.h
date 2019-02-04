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

PairStyle(spin/magelec,PairSpinMagelec)

#else

#ifndef LMP_PAIR_SPIN_MAGELEC_H
#define LMP_PAIR_SPIN_MAGELEC_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinMagelec : public PairSpin {
 public:
  PairSpinMagelec(class LAMMPS *);
  virtual ~PairSpinMagelec();
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void *extract(const char *, int &);

  void compute(int, int);
  void compute_single_pair(int, double *);

  void compute_magelec(int, int, double *, double *, double *);
  void compute_magelec_mech(int, int, double *, double *, double *);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  double cut_spin_magelec_global;	// global me cutoff

 protected:
  double **ME, **ME_mech;		// magelec coeff in eV
  double **v_mex, **v_mey, **v_mez;	// magelec direction
  double **cut_spin_magelec;		// magelec cutoff distance

  int lattice_flag;                     // flag for mech force computation
  class FixNVESpin *lockfixnvespin;     // ptr to FixNVESpin for setups

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
