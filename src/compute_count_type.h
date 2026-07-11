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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(count/type,ComputeCountType);
// clang-format on
#else

#ifndef LMP_COMPUTE_COMPUTE_TYPE_H
#define LMP_COMPUTE_COMPUTE_TYPE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCountType : public Compute {
 public:
  ComputeCountType(class LAMMPS *, int, char **);
  ~ComputeCountType() override;
  void init() override {}
  double compute_scalar() override;
  void compute_vector() override;

 protected:
  int mode;

  int *count;
  bigint *bcount_me;
  bigint *bcount;

  int count_atoms();
  int count_bonds();
  int count_angles();
  int count_dihedrals();
  int count_impropers();
};

}    // namespace LAMMPS_NS

#endif
#endif
