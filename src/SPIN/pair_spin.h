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

#ifndef LMP_PAIR_SPIN_H
#define LMP_PAIR_SPIN_H

#include "pair.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class PairSpin : public Pair {
  friend class FixNVESpin;

 public:
  PairSpin(class LAMMPS *);

  void settings(int, char **) override;
  void coeff(int, char **) override {}
  void init_style() override;
  double init_one(int, int) override { return 0.0; }
  void *extract(const char *, int &) override { return nullptr; }

  void compute(int, int) override {}
  virtual void compute_single_pair(int, double *) {}

  // storing magnetic energies

  int nlocal_max;    // max nlocal (for list size)
  double *emag;      // energy list

 protected:
  double hbar;         // Planck constant (eV.ps.rad-1)
  int lattice_flag;    // flag for mech force computation

  virtual void allocate() {}
};

}    // namespace LAMMPS_NS

#endif
