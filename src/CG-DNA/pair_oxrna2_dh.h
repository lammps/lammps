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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxrna2/dh,PairOxrna2Dh);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_DH_H
#define LMP_PAIR_OXRNA2_DH_H

#include "pair_oxdna2_dh.h"

namespace LAMMPS_NS {

class PairOxrna2Dh : public PairOxdna2Dh {
 public:
  PairOxrna2Dh(class LAMMPS *lmp) : PairOxdna2Dh(lmp) {}

  void compute_interaction_sites(double *, double *, double *, double *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
