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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(hybrid/overlay,PairHybridOverlay);
// clang-format on
#else

#ifndef LMP_PAIR_HYBRID_OVERLAY_H
#define LMP_PAIR_HYBRID_OVERLAY_H

#include "pair_hybrid.h"

namespace LAMMPS_NS {

class PairHybridOverlay : public PairHybrid {
 public:
  PairHybridOverlay(class LAMMPS *);

  void coeff(int, char **) override;

  void init_svector() override;
  void copy_svector(int, int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
