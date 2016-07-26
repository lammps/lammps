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

#ifdef NEIGH_PAIR_CLASS

NeighPairStyle(NEIGH_PAIR_FULL_NSQ_GHOST,NeighPairFullNsqGhost)

#else

#ifndef LMP_NEIGH_PAIR_FULL_NSQ_GHOST_H
#define LMP_NEIGH_PAIR_FULL_NSQ_GHOST_H

#include "neigh_pair.h"

namespace LAMMPS_NS {

class NeighPairFullNsqGhost : public NeighPair {
 public:
  NeighPairFullNsqGhost(class LAMMPS *);
  ~NeighPairFullNsqGhost() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
