/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mliap,PairMLIAP);
// clang-format on
#else

#ifndef LMP_PAIR_MLIAP_H
#define LMP_PAIR_MLIAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMLIAP : public Pair {
 public:
  PairMLIAP(class LAMMPS *);
  ~PairMLIAP() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void e_tally(class MLIAPData *);
  void v_tally(int, int, double *, double *);
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;
  int *map;    // mapping from atom types to elements

 protected:
  virtual void allocate();

  class MLIAPModel *model;
  class MLIAPDescriptor *descriptor;
  class MLIAPData *data;
  bool is_child;
};

}    // namespace LAMMPS_NS

#endif
#endif
