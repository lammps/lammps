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
typedef NPairHalffullOmp<0, 0, 0> NPairHalffullNewtoffOmp;
NPairStyle(halffull/newtoff/omp,
           NPairHalffullNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_OMP);

typedef NPairHalffullOmp<0, 0, 0> NPairHalffullNewtoffOmp;
NPairStyle(halffull/newtoff/skip/omp,
           NPairHalffullNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_OMP);

typedef NPairHalffullOmp<0, 0, 0> NPairHalffullNewtoffOmp;
NPairStyle(halffull/newtoff/ghost/omp,
           NPairHalffullNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_OMP);

typedef NPairHalffullOmp<0, 0, 0> NPairHalffullNewtoffOmp;
NPairStyle(halffull/newtoff/skip/ghost/omp,
           NPairHalffullNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST | NP_OMP);

typedef NPairHalffullOmp<1, 0, 0> NPairHalffullNewtonOmp;
NPairStyle(halffull/newton/omp,
           NPairHalffullNewtonOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_OMP);

typedef NPairHalffullOmp<1, 1, 0> NPairHalffullNewtonTriOmp;
NPairStyle(halffull/newton/tri/omp,
           NPairHalffullNewtonTriOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI | NP_OMP);

typedef NPairHalffullOmp<1, 0, 0> NPairHalffullNewtonOmp;
NPairStyle(halffull/newton/skip/omp,
           NPairHalffullNewtonOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_SKIP | NP_OMP);

typedef NPairHalffullOmp<1, 1, 0> NPairHalffullNewtonTriOmp;
NPairStyle(halffull/newton/tri/skip/omp,
           NPairHalffullNewtonTriOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI | NP_SKIP | NP_OMP);

typedef NPairHalffullOmp<0, 0, 1> NPairHalffullTrimNewtoffOmp;
NPairStyle(halffull/trim/newtoff/omp,
           NPairHalffullTrimNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<0, 0, 1> NPairHalffullTrimNewtoffOmp;
NPairStyle(halffull/trim/newtoff/skip/omp,
           NPairHalffullTrimNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<0, 0, 1> NPairHalffullTrimNewtoffOmp;
NPairStyle(halffull/trim/newtoff/ghost/omp,
           NPairHalffullTrimNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<0, 0, 1> NPairHalffullTrimNewtoffOmp;
NPairStyle(halffull/trim/newtoff/skip/ghost/omp,
           NPairHalffullTrimNewtoffOmp,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<1, 0, 1> NPairHalffullTrimNewtonOmp;
NPairStyle(halffull/trim/newton/omp,
           NPairHalffullTrimNewtonOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<1, 1, 1> NPairHalffullTrimNewtonTriOmp;
NPairStyle(halffull/trim/newton/tri/omp,
           NPairHalffullTrimNewtonTriOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<1, 0, 1> NPairHalffullTrimNewtonOmp;
NPairStyle(halffull/trim/newton/skip/omp,
           NPairHalffullTrimNewtonOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_SKIP | NP_TRIM | NP_OMP);

typedef NPairHalffullOmp<1, 1, 1> NPairHalffullTrimNewtonTriOmp;
NPairStyle(halffull/trim/newton/tri/skip/omp,
           NPairHalffullTrimNewtonTriOmp,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI | NP_SKIP | NP_TRIM | NP_OMP);
// clang-format on
#else

#ifndef LMP_NPAIR_HALFFULL_OMP_H
#define LMP_NPAIR_HALFFULL_OMP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int NEWTON, int TRI, int TRIM>
class NPairHalffullOmp : public NPair {
 public:
  NPairHalffullOmp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
