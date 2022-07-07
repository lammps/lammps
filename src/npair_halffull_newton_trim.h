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
NPairStyle(halffull/newton/trim,
           NPairHalffullNewtonTrim,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_TRIM);

NPairStyle(halffull/newton/skip/trim,
           NPairHalffullNewtonTrim,
           NP_HALF_FULL | NP_NEWTON | NP_HALF | NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_ORTHO | NP_TRI | NP_SKIP | NP_TRIM);
// clang-format on
#else

#ifndef LMP_NPAIR_HALFFULL_NEWTON_TRIM_H
#define LMP_NPAIR_HALFFULL_NEWTON_TRIM_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairHalffullNewtonTrim : public NPair {
 public:
  NPairHalffullNewtonTrim(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
