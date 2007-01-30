/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_EAM_ALLOY_H
#define PAIR_EAM_ALLOY_H

#include "pair_eam.h"

namespace LAMMPS_NS {

// use virtual public since this class is parent in multiple inheritance

class PairEAMAlloy : virtual public PairEAM {
 public:
  PairEAMAlloy(class LAMMPS *);
  virtual ~PairEAMAlloy() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}

#endif
