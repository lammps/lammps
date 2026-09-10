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

#include "nucleotide_oxdna.h"
#include "pair_oxdna_excv.h"

namespace LAMMPS_NS {

class PairOxrna2Excv : public PairOxdnaExcv {
 public:
  PairOxrna2Excv(class LAMMPS *lmp) : PairOxdnaExcv(lmp) {}
  // inline below has to be here in the header file, otherwise KOKKOS
  // compilation fails due to undefined vtable symbols.
  void compute_backbone_site(double e1[3], double /*e2*/[3], double e3[3],
                             double rbk[3]) const override
  {
    NucleotideOxrna2 oxrna2;
    oxrna2.backbone_site(e1, nullptr, e3, rbk);
  };
};

}    // namespace LAMMPS_NS

#endif
#endif
