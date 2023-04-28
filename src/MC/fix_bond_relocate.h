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

#ifdef FIX_CLASS
// clang-format off
FixStyle(bond/relocate,FixBondRelocate);
// clang-format on
#else

#ifndef LMP_FIX_BOND_RELOCATE_H
#define LMP_FIX_BOND_RELOCATE_H

#include "fix.h"
#include "fix_bond_create.h"

namespace LAMMPS_NS {

class FixBondRelocate : public Fix {
 public:
  FixBondRelocate(class LAMMPS *, int, char **);
  ~FixBondRelocate() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void post_integrate() override;
  int modify_param(int, char **) override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  double fraction, locutsq, hicutsq;
  tagint relocate_btype;
  tagint decide_btype;
  int decide_functionality;
  int nmax, tflag;
  int *randomized_eligible_sourceatoms;
  int *randomized_eligible_targetatoms;
  int naccept, foursome;
  int angleflag;
  char *id_temp;
  int *type;
  double **x;
  //   FixBondCreate *create_fix;

  int maxpermute;
  int *permute;

  class NeighList *list;
  class Compute *temperature;
  class RanMars *random;

  double dist_rsq(int, int);
  double pair_eng(int, int);
  double bond_eng(int, int, int);
  double angle_eng(int, int, int, int);

  void neighbor_permutation(int);

  void break_bond(int, int);
  void create_bond(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
