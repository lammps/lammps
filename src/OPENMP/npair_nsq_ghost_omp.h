/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NPAIR_CLASS
// clang-format off
typedef NPairNsqGhostOmp<0> NPairFullNsqGhostOmp;
NPairStyle(full/nsq/ghost/omp,
           NPairFullNsqGhostOmp,
           NP_FULL | NP_NSQ | NP_NEWTON | NP_NEWTOFF | NP_GHOST | NP_OMP | NP_ORTHO | NP_TRI);

typedef NPairNsqGhostOmp<1> NPairHalfNsqNewtoffGhostOmp;
NPairStyle(half/nsq/newtoff/ghost/omp,
           NPairHalfNsqNewtoffGhostOmp,
           NP_HALF | NP_NSQ | NP_NEWTOFF | NP_GHOST | NP_OMP | NP_ORTHO | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_NSQ_GHOST_OMP_H
#define LMP_NPAIR_NSQ_GHOST_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF>
class NPairNsqGhostOmp : public NPair {
 public:
  NPairNsqGhostOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
