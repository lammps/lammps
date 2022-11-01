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
FixStyle(bond/swap,FixBondSwap);
// clang-format on
#else

#ifndef LMP_FIX_BONDSWAP_H
#define LMP_FIX_BONDSWAP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondSwap : public Fix {
 public:
  FixBondSwap(class LAMMPS *, int, char **);
  ~FixBondSwap() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void post_integrate() override;
  int modify_param(int, char **) override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  double fraction, cutsq;
  int nmax, tflag;
  int *alist;
  int naccept, foursome;
  int angleflag;
  char *id_temp;
  int *type;
  double **x;

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
};

}    // namespace LAMMPS_NS

#endif
#endif
