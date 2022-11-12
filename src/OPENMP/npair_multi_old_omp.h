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
typedef NPairMultiOldOmp<0, 1, 0, 0> NPairFullMultiOldOmp;
NPairStyle(full/multi/old/omp,
           NPairFullMultiOldOmp,
           NP_FULL | NP_MULTI_OLD | NP_OMP |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOldOmp<1, 0, 0, 0> NPairHalfMultiOldNewtoffOmp;
NPairStyle(half/multi/old/newtoff/omp,
           NPairHalfMultiOldNewtoffOmp,
           NP_HALF | NP_MULTI_OLD | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOldOmp<1, 1, 0, 0> NPairHalfMultiOldNewtonOmp;
NPairStyle(half/multi/old/newton/omp,
           NPairHalfMultiOldNewtonOmp,
           NP_HALF | NP_MULTI_OLD | NP_OMP | NP_NEWTON | NP_ORTHO);

typedef NPairMultiOldOmp<1, 1, 1, 0> NPairHalfMultiOldNewtonTriOmp;
NPairStyle(half/multi/old/newton/tri/omp,
           NPairHalfMultiOldNewtonTriOmp,
           NP_HALF | NP_MULTI_OLD | NP_OMP | NP_NEWTON | NP_TRI);

typedef NPairMultiOldOmp<0, 1, 0, 1> NPairFullSizeMultiOldOmp;
NPairStyle(full/size/multi/old/omp,
           NPairFullSizeMultiOldOmp,
           NP_FULL | NP_SIZE | NP_MULTI_OLD | NP_OMP |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOldOmp<1, 0, 0, 1> NPairHalfSizeMultiOldNewtoffOmp;
NPairStyle(half/size/multi/old/newtoff/omp,
           NPairHalfSizeMultiOldNewtoffOmp,
           NP_HALF | NP_SIZE | NP_MULTI_OLD | NP_OMP | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairMultiOldOmp<1, 1, 0, 1> NPairHalfSizeMultiOldNewtonOmp;
NPairStyle(half/size/multi/old/newton/omp,
           NPairHalfSizeMultiOldNewtonOmp,
           NP_HALF | NP_SIZE | NP_MULTI_OLD | NP_OMP | NP_NEWTON | NP_ORTHO);

typedef NPairMultiOldOmp<1, 1, 1, 1> NPairHalfSizeMultiOldNewtonTriOmp;
NPairStyle(half/size/multi/old/newton/tri/omp,
           NPairHalfSizeMultiOldNewtonTriOmp,
           NP_HALF | NP_SIZE | NP_MULTI_OLD | NP_OMP | NP_NEWTON | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_MULTI_OLD_OMP_H
#define LMP_NPAIR_MULTI_OLD_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int HALF, int NEWTON, int TRI, int SIZE>
class NPairMultiOldOmp : public NPair {
 public:
  NPairMultiOldOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
