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

#ifdef NPAIR_CLASS
// clang-format off
typedef NPairBinAtomonly<0, 1, 0> NPairFullBinAtomonly;
NPairStyle(full/bin/atomonly,
           NPairFullBinAtomonly,
           NP_FULL | NP_BIN | NP_ATOMONLY |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairBinAtomonly<1, 1, 0> NPairHalfBinAtomonlyNewton;
NPairStyle(half/bin/atomonly/newton,
           NPairHalfBinAtomonlyNewton,
           NP_HALF | NP_BIN | NP_ATOMONLY | NP_NEWTON | NP_ORTHO);
// clang-format on
#else

#ifndef LMP_NPAIR_BIN_ATOMONLY_H
#define LMP_NPAIR_BIN_ATOMONLY_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI>
class NPairBinAtomonly : public NPair {
 public:
  NPairBinAtomonly(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Neighbor list overflow, boost neigh_modify one

UNDOCUMENTED

*/
