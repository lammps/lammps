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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(eam/fs/opt,PairEAMFSOpt);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_FS_OPT_H
#define LMP_PAIR_EAM_FS_OPT_H

#include "pair_eam_fs.h"
#include "pair_eam_opt.h"

namespace LAMMPS_NS {

class PairEAMFSOpt : public PairEAMFS, public PairEAMOpt {
 public:
  PairEAMFSOpt(class LAMMPS *);
  virtual ~PairEAMFSOpt() {}
};

}    // namespace LAMMPS_NS

#endif
#endif
