/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(eam/fs/intel,PairEAMFSIntel)

#else

#ifndef LMP_PAIR_EAM_FS_INTEL_H
#define LMP_PAIR_EAM_FS_INTEL_H

#include "pair_eam_intel.h"

namespace LAMMPS_NS {

// need virtual public b/c of how eam/fs/opt inherits from it

class PairEAMFSIntel : virtual public PairEAMIntel {
 public:
  PairEAMFSIntel(class LAMMPS *);
  virtual ~PairEAMFSIntel() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}

#endif
#endif
