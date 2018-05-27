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

#ifdef NPAIR_CLASS

NPairStyle(half/multi/newton/tri,
           NPairHalfMultiNewtonTri,
           NP_HALF | NP_MULTI | NP_NEWTON | NP_TRI)

#else

#ifndef LMP_NPAIR_HALF_MULTI_NEWTON_TRI_H
#define LMP_NPAIR_HALF_MULTI_NEWTON_TRI_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairHalfMultiNewtonTri : public NPair {
 public:
  NPairHalfMultiNewtonTri(class LAMMPS *);
  ~NPairHalfMultiNewtonTri() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Neighbor list overflow, boost neigh_modify one

UNDOCUMENTED

*/
