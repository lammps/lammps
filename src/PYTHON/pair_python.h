/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Pair zero is a dummy pair interaction useful for requiring a
   force cutoff distance in the absence of pair-interactions or
   with hybrid/overlay if a larger force cutoff distance is required.

   This can be used in conjunction with bond/create to create bonds
   that are longer than the cutoff of a given force field, or to
   calculate radial distribution functions for models without
   pair interactions.

------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(python,PairPython);
// clang-format on
#else

#ifndef LMP_PAIR_PYTHON_H
#define LMP_PAIR_PYTHON_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPython : public Pair {
 public:
  PairPython(class LAMMPS *);
  ~PairPython() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double cut_global;
  void *py_potential;
  int *skip_types;

  virtual void allocate();
  void *get_member_function(const char *);
};

}    // namespace LAMMPS_NS

#endif
#endif
