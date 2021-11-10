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
PairStyle(eam/fs/omp,PairEAMFSOMP);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_FS_OMP_H
#define LMP_PAIR_EAM_FS_OMP_H

#include "pair_eam_omp.h"

namespace LAMMPS_NS {

// need virtual public b/c of how eam/fs/opt inherits from it

class PairEAMFSOMP : virtual public PairEAMOMP {
 public:
  PairEAMFSOMP(class LAMMPS *);
  virtual ~PairEAMFSOMP() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}    // namespace LAMMPS_NS

#endif
#endif
