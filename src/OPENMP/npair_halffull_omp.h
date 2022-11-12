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
typedef NPairHalffullOmp<0, 0> NPairHalffullOmpNewtoffOmp;
NPairStyle(halffull/newtoff/omp,
           NPairHalffullOmpNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_OMP);

typedef NPairHalffullOmp<0, 0> NPairHalffullOmpNewtoffOmp;
NPairStyle(halffull/newtoff/skip/omp,
           NPairHalffullOmpNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_OMP);

typedef NPairHalffullOmp<0, 0> NPairHalffullOmpNewtoffOmp;
NPairStyle(halffull/newtoff/ghost/omp,
           NPairHalffullOmpNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_OMP);

typedef NPairHalffullOmp<0, 0> NPairHalffullOmpNewtoffOmp;
NPairStyle(halffull/newtoff/skip/ghost/omp,
           NPairHalffullOmpNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST | NP_OMP);

typedef NPairHalffullOmp<1, 0> NPairHalffullOmpNewtonOmp;
NPairStyle(halffull/newton/omp,
           NPairHalffullOmpNewtonOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_OMP);

typedef NPairHalffullOmp<1, 0> NPairHalffullOmpNewtonOmp;
NPairStyle(halffull/newton/skip/omp,
           NPairHalffullOmpNewtonOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_OMP);

typedef NPairHalffullOmp<0, 1> NPairHalffullOmpNewtoffTrimOmp;
NPairStyle(halffull/newtoff/trim/omp,
           NPairHalffullOmpNewtoffTrimOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<0, 1> NPairHalffullOmpNewtoffTrimOmp;
NPairStyle(halffull/newtoff/skip/trim/omp,
           NPairHalffullOmpNewtoffTrimOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<0, 1> NPairHalffullOmpNewtoffTrimOmp;
NPairStyle(halffull/newtoff/ghost/trim/omp,
           NPairHalffullOmpNewtoffTrimOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<0, 1> NPairHalffullOmpNewtoffTrimOmp;
NPairStyle(halffull/newtoff/skip/ghost/trim/omp,
           NPairHalffullOmpNewtoffTrimOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<1, 1> NPairHalffullOmpNewtonTrimOmp;
NPairStyle(halffull/newton/trim/omp,
           NPairHalffullOmpNewtonTrimOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<1, 1> NPairHalffullOmpNewtonTrimOmp;
NPairStyle(halffull/newton/skip/trim/omp,
           NPairHalffullOmpNewtonTrimOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM | NP_OMP);
// clang-format on
#else

#ifndef LMP_NPAIR_HALFFULL_OMP_H
#define LMP_NPAIR_HALFFULL_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int NEWTON, int TRIM>
class NPairHalffullOmp : public NPair {
 public:
  NPairHalffullOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
