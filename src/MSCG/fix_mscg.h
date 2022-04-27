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

#ifdef FIX_CLASS
// clang-format off
FixStyle(mscg,FixMSCG);
// clang-format on
#else

#ifndef LMP_FIX_MSCG_H
#define LMP_FIX_MSCG_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMSCG : public Fix {
 public:
  FixMSCG(class LAMMPS *, int, char **);
  ~FixMSCG() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void end_of_step() override;
  void post_run() override;

 private:
  int range_flag, name_flag, me, nprocs;
  int nframes, n_frames, block_size, n_cg_sites, n_cg_types, *cg_site_types;
  int max_partners_bond, max_partners_angle, max_partners_dihedral;
  unsigned *n_partners_bond, *n_partners_angle, *n_partners_dihedral;
  unsigned **partners_bond, **partners_angle, **partners_dihedral;
  double *x1d, *f1d, **f;
  double box_half_lengths[3];
  char **type_names;
  void *mscg_struct;
};

}    // namespace LAMMPS_NS

#endif
#endif
