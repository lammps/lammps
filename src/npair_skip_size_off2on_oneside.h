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
NPairStyle(skip/size/off2on/oneside,
           NPairSkipSizeOff2onOneside,
           NP_SKIP | NP_SIZE | NP_OFF2ON | NP_ONESIDE | NP_HALF |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD | NP_NEWTON | NP_NEWTOFF |
           NP_ORTHO | NP_TRI);
// clang-format on
#else

#ifndef LMP_NPAIR_SKIP_SIZE_OFF2ON_ONESIDE_H
#define LMP_NPAIR_SKIP_SIZE_OFF2ON_ONESIDE_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairSkipSizeOff2onOneside : public NPair {
 public:
  NPairSkipSizeOff2onOneside(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
