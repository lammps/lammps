/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(meam/ms,PairMEAMMS);
// clang-format on
#else

#ifndef LMP_PAIR_MEAM_MS_H
#define LMP_PAIR_MEAM_MS_H

#include "pair_meam.h"

namespace LAMMPS_NS {

class PairMEAMMS : public PairMEAM {
 public:
  PairMEAMMS(class LAMMPS *);
};
}    // namespace LAMMPS_NS
#endif
#endif
