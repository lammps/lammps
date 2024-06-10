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
typedef NPairSkipTemp<0> NPairSkip;
NPairStyle(skip,
           NPairSkip,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairSkipTemp<0> NPairSkip;
NPairStyle(skip/ghost,
           NPairSkip,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_GHOST);

typedef NPairSkipTemp<0> NPairSkipSize;
NPairStyle(skip/half/size,
           NPairSkipSize,
           NP_SKIP | NP_SIZE | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI);

typedef NPairSkipTemp<1> NPairSkipTrim;
NPairStyle(skip/trim,
           NPairSkipTrim,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_TRIM);

typedef NPairSkipTemp<1> NPairSkipTrim;
NPairStyle(skip/ghost/trim,
           NPairSkipTrim,
           NP_SKIP | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_GHOST | NP_TRIM);

typedef NPairSkipTemp<1> NPairSkipTrimSize;
NPairStyle(skip/trim/half/size,
           NPairSkipTrimSize,
           NP_SKIP | NP_SIZE | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_TRIM);

// clang-format on
#else

#ifndef LMP_NPAIR_SKIP_H
#define LMP_NPAIR_SKIP_H

#include "npair.h"

namespace LAMMPS_NS {

template<int TRIM>
class NPairSkipTemp : public NPair {
 public:
  NPairSkipTemp(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
