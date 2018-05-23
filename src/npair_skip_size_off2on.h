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

NPairStyle(skip/size/off2on,
           NPairSkipSizeOff2on,
           NP_SKIP | NP_SIZE | NP_OFF2ON | NP_HALF |
           NP_NSQ | NP_BIN | NP_MULTI |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI)

#else

#ifndef LMP_NPAIR_SKIP_SIZE_OFF2ON_H
#define LMP_NPAIR_SKIP_SIZE_OFF2ON_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairSkipSizeOff2on : public NPair {
 public:
  NPairSkipSizeOff2on(class LAMMPS *);
  ~NPairSkipSizeOff2on() {}
  void build(class NeighList *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Neighbor list overflow, boost neigh_modify one

UNDOCUMENTED

*/
