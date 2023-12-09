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
typedef NPairHalffull<0, 0, 0> NPairHalffullNewtoff;
NPairStyle(halffull/newtoff,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI);

typedef NPairHalffull<0, 0, 0> NPairHalffullNewtoff;
NPairStyle(halffull/newtoff/skip,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP);

typedef NPairHalffull<0, 0, 0> NPairHalffullNewtoff;
NPairStyle(halffull/newtoff/ghost,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST);

typedef NPairHalffull<0, 0, 0> NPairHalffullNewtoff;
NPairStyle(halffull/newtoff/skip/ghost,
           NPairHalffullNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST);

typedef NPairHalffull<1, 0, 0> NPairHalffullNewton;
NPairStyle(halffull/newton,
           NPairHalffullNewton,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO);

typedef NPairHalffull<1, 1, 0> NPairHalffullNewtonTri;
NPairStyle(halffull/newton/tri,
           NPairHalffullNewtonTri,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI);

typedef NPairHalffull<1, 0, 0> NPairHalffullNewton;
NPairStyle(halffull/newton/skip,
           NPairHalffullNewton,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_SKIP);

typedef NPairHalffull<1, 1, 0> NPairHalffullNewtonTri;
NPairStyle(halffull/newton/skip/tri,
           NPairHalffullNewtonTri,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI | NP_SKIP);

typedef NPairHalffull<0, 0, 1> NPairHalffullTrimNewtoff;
NPairStyle(halffull/trim/newtoff,
           NPairHalffullTrimNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_TRIM);

typedef NPairHalffull<0, 0, 1> NPairHalffullTrimNewtoff;
NPairStyle(halffull/trim/newtoff/skip,
           NPairHalffullTrimNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM);

typedef NPairHalffull<0, 0, 1> NPairHalffullTrimNewtoff;
NPairStyle(halffull/trim/newtoff/ghost,
           NPairHalffullTrimNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM);

typedef NPairHalffull<0, 0, 1> NPairHalffullTrimNewtoff;
NPairStyle(halffull/trim/newtoff/skip/ghost,
           NPairHalffullTrimNewtoff,
           NP_HALF_FULL | NP_NEWTOFF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_HALF |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_GHOST | NP_TRIM);

typedef NPairHalffull<1, 0, 1> NPairHalffullTrimNewton;
NPairStyle(halffull/trim/newton,
           NPairHalffullTrimNewton,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRIM);

typedef NPairHalffull<1, 1, 1> NPairHalffullTrimNewtonTri;
NPairStyle(halffull/trim/newton/tri,
           NPairHalffullTrimNewtonTri,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI | NP_TRIM);

typedef NPairHalffull<1, 0, 1> NPairHalffullTrimNewton;
NPairStyle(halffull/trim/newton/skip,
           NPairHalffullTrimNewton,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_SKIP | NP_TRIM);

typedef NPairHalffull<1, 1, 1> NPairHalffullTrimNewtonTri;
NPairStyle(halffull/trim/newton/tri/skip,
           NPairHalffullTrimNewtonTri,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_TRI | NP_SKIP | NP_TRIM);
// clang-format on
#else

#ifndef LMP_NPAIR_HALFFULL_H
#define LMP_NPAIR_HALFFULL_H

#include "npair.h"

namespace LAMMPS_NS {

template<int NEWTON, int TRI, int TRIM>
class NPairHalffull : public NPair {
 public:
  NPairHalffull(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
