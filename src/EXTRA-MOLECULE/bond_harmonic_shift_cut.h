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
BondStyle(harmonic/shift/cut,BondHarmonicShiftCut);
// clang-format on
#else

#ifndef LMP_BOND_HARMONIC_SHIFT_CUT_H
#define LMP_BOND_HARMONIC_SHIFT_CUT_H

#include "bond.h"

namespace LAMMPS_NS {

class BondHarmonicShiftCut : public Bond {
 public:
  BondHarmonicShiftCut(class LAMMPS *);
  virtual ~BondHarmonicShiftCut();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);

 protected:
  double *k, *r0, *r1;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
