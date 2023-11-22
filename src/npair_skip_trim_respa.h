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
NPairStyle(skip/trim/half/respa,
           NPairSkipTrimRespa,
           NP_SKIP | NP_RESPA | NP_HALF | NP_FULL |
           NP_NSQ | NP_BIN | NP_MULTI | NP_MULTI_OLD |
           NP_NEWTON | NP_NEWTOFF | NP_ORTHO | NP_TRI | NP_TRIM);
// clang-format on
#else

#ifndef LMP_NPAIR_SKIP_TRIM_RESPA_H
#define LMP_NPAIR_SKIP_TRIM_RESPA_H

#include "npair.h"

namespace LAMMPS_NS {

class NPairSkipTrimRespa : public NPair {
 public:
  NPairSkipTrimRespa(class LAMMPS *);
  void build(class NeighList *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
