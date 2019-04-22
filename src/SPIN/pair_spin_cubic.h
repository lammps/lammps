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

PairStyle(spin/cubic,PairSpinCubic)

#else

#ifndef LMP_PAIR_SPIN_CUBIC_H
#define LMP_PAIR_SPIN_CUBIC_H

#include "pair_spin.h"

namespace LAMMPS_NS {

class PairSpinCubic : public PairSpin {
 public:
  PairSpinCubic(class LAMMPS *);
  virtual ~PairSpinCubic();
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void *extract(const char *, int &);

  void compute(int, int);
  void compute_single_pair(int, double *);

  void compute_cubic(int, double *, double *, double *, double *, double *);
  void compute_cubic_mech(int, double *, double *, double *, double *, double *, double *);
  double compute_cubic_energy(int, double *, double *, double *, double *);

  void set_axis(int, double *, double *, double *, double *);

  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

  double cut_spin_cubic_global;		// global cutoff distance

 protected:
  double **K1_mag;			// cubic aniso coeffs in eV
  double **K2_mag;			
  double **K1_mech;			// mech coeffs coeffs in
  double **K2_mech;			
  double **cut_spin_cubic;		// cutoff distance

  //double *ea1, *ea2, *ea3;

  int lattice_flag; 			// flag for mech force computation
  class FixNVESpin *lockfixnvespin;	// ptr to FixNVESpin for setups

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
