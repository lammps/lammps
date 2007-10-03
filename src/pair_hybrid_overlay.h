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

#ifndef PAIR_HYBRID_OVERLAY_H
#define PAIR_HYBRID_OVERLAY_H

#include "pair_hybrid.h"

namespace LAMMPS_NS {

class PairHybridOverlay : public PairHybrid {
 public:
  PairHybridOverlay(class LAMMPS *);
  void coeff(int, char **);

 private:
  void modify_requests();
};

}

#endif
