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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(group/group,ComputeGroupGroup);
// clang-format on
#else

#ifndef LMP_COMPUTE_GROUP_GROUP_H
#define LMP_COMPUTE_GROUP_GROUP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGroupGroup : public Compute {
 public:
  ComputeGroupGroup(class LAMMPS *, int, char **);
  ~ComputeGroupGroup() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  double compute_scalar() override;
  void compute_vector() override;

 private:
  char *group2;
  int jgroup, jgroupbit, othergroupbit;
  double **cutsq;
  double e_self, e_correction;
  int pairflag, kspaceflag, boundaryflag, molflag;
  class Pair *pair;
  class NeighList *list;
  class KSpace *kspace;

  void pair_contribution();
  void kspace_contribution();
  void kspace_correction();
};

}    // namespace LAMMPS_NS

#endif
#endif
