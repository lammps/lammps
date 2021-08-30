// clang-format off
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
PairStyle(eam/alloy/intel,PairEAMAlloyIntel);
// clang-format on
#else

#ifndef LMP_PAIR_EAM_ALLOY_INTEL_H
#define LMP_PAIR_EAM_ALLOY_INTEL_H

#include "pair_eam_intel.h"

namespace LAMMPS_NS {

// need virtual public b/c of how eam/alloy/opt inherits from it

class PairEAMAlloyIntel : virtual public PairEAMIntel {
 public:
  PairEAMAlloyIntel(class LAMMPS *);
  virtual ~PairEAMAlloyIntel() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}    // namespace LAMMPS_NS

#endif
#endif
