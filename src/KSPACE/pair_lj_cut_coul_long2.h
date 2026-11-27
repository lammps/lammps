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

/* ----------------------------------------------------------------------
   lj/cut/coul/long2 -- CPU base for the experimental lj/cut/coul/long2/kk
   accelerator style.  Physics is identical to lj/cut/coul/long; this class
   exists so the KOKKOS style-registration machinery has a base style to
   attach the /kk variant to.  See src/KOKKOS/pair_lj_cut_coul_long2_kokkos.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/long2,PairLJCutCoulLong2);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_LONG2_H
#define LMP_PAIR_LJ_CUT_COUL_LONG2_H

#include "pair_lj_cut_coul_long.h"

namespace LAMMPS_NS {

class PairLJCutCoulLong2 : public PairLJCutCoulLong {
 public:
  PairLJCutCoulLong2(class LAMMPS *lmp) : PairLJCutCoulLong(lmp) {}
};

}    // namespace LAMMPS_NS

#endif
#endif
