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
PairStyle(oxrna2/excv,PairOxrna2Excv);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_EXCV_H
#define LMP_PAIR_OXRNA2_EXCV_H

#include "pair_oxdna_excv.h"

namespace LAMMPS_NS {

class PairOxrna2Excv : public PairOxdnaExcv {
 public:
  PairOxrna2Excv(class LAMMPS *lmp) : PairOxdnaExcv(lmp) {}

  void compute_interaction_sites(double *, double *, double *, double *, double *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
