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
typedef NPairMultiOld<0, 1, 0, 0> NPairFullMultiOld;
NPairStyle(full/multi/old,
           NPairFullMultiOld,
           NP_FULL | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOld<1, 0, 0, 0> NPairHalfMultiOldNewtoff;
NPairStyle(half/multi/old/newtoff,
           NPairHalfMultiOldNewtoff,
           NP_HALF | NP_MULTI_OLD | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOld<1, 1, 0, 0> NPairHalfMultiOldNewton;
NPairStyle(half/multi/old/newton,
           NPairHalfMultiOldNewton,
           NP_HALF | NP_MULTI_OLD | NP_NEWTON | NP_ORTHO);

typedef NPairMultiOld<1, 1, 1, 0> NPairHalfMultiOldNewtonTri;
NPairStyle(half/multi/old/newton/tri,
           NPairHalfMultiOldNewtonTri,
           NP_HALF | NP_MULTI_OLD | NP_NEWTON | NP_TRI);

typedef NPairMultiOld<0, 1, 0, 1> NPairFullSizeMultiOld;
NPairStyle(full/size/multi/old,
           NPairFullSizeMultiOld,
           NP_FULL | NP_SIZE | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOld<1, 0, 0, 1> NPairHalfSizeMultiOldNewtoff;
NPairStyle(half/size/multi/old/newtoff,
           NPairHalfSizeMultiOldNewtoff,
           NP_HALF | NP_SIZE | NP_MULTI_OLD | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOld<1, 1, 0, 1> NPairHalfSizeMultiOldNewton;
NPairStyle(half/size/multi/old/newton,
           NPairHalfSizeMultiOldNewton,
           NP_HALF | NP_SIZE | NP_MULTI_OLD | NP_NEWTON | NP_ORTHO);

typedef NPairMultiOld<1, 1, 1, 1> NPairHalfSizeMultiOldNewtonTri;
NPairStyle(half/size/multi/old/newton/tri,
           NPairHalfSizeMultiOldNewtonTri,
           NP_HALF | NP_SIZE | NP_MULTI_OLD | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_MULTI_OLD_H
#define LMP_NPAIR_MULTI_OLD_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE>
class NPairMultiOld : public NPair {
 public:
  NPairMultiOld(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
