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
  ~ComputeGroupGroup();
  void init();
  void init_list(int, class NeighList *);
  double compute_scalar();
  void compute_vector();

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute group/group group ID does not exist

Self-explanatory.

E: Compute group/group molecule requires molecule IDs

UNDOCUMENTED

E: No pair style defined for compute group/group

Cannot calculate group interactions without a pair style defined.

E: Pair style does not support compute group/group

The pair_style does not have a single() function, so it cannot be
invoked by the compute group/group command.

E: No Kspace style defined for compute group/group

Self-explanatory.

E: Kspace style does not support compute group/group

Self-explanatory.

W: Both groups in compute group/group have a net charge; the Kspace boundary correction to energy will be non-zero

Self-explanatory.

*/
