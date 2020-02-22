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

#ifndef LMP_PAIR_SPIN_H
#define LMP_PAIR_SPIN_H

#include "pair.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class PairSpin : public Pair {
friend class FixNVESpin;
 public:
  PairSpin(class LAMMPS *);
  virtual ~PairSpin();
  virtual void settings(int, char **);
  virtual void coeff(int, char **) {}
  virtual void init_style();
  virtual double init_one(int, int) {return 0.0;}
  virtual void *extract(const char *, int &) {return NULL;}

  virtual void compute(int, int) {}
  virtual void compute_single_pair(int, double *) {}
  
  // test emag list storing mag energies
  int nlocal_max;                       // max value of nlocal (for size of lists)
  double *emag;                         // energy list

 protected:
  double hbar;                          // Planck constant (eV.ps.rad-1)
  int lattice_flag;                     // flag for mech force computation

  virtual void allocate() {}
};

}

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
