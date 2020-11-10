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

NPairStyle(half/bytype/newton/tri,
           NPairHalfBytypeNewtonTri,
           NP_HALF | NP_BYTYPE | NP_NEWTON | NP_TRI)

#else

#ifndef LMP_NPAIR_HALF_BYTYPE_NEWTON_TRI_H
#define LMP_NPAIR_HALF_BYTYPE_NEWTON_TRI_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairHalfBytypeNewtonTri : public NPair {
 public:
  NPairHalfBytypeNewtonTri(class LAMMPS *);
  ~NPairHalfBytypeNewtonTri() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
