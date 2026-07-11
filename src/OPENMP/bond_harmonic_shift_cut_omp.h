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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(harmonic/shift/cut/omp,BondHarmonicShiftCutOMP);
// clang-format on
#else

#ifndef LMP_BOND_HARMONIC_SHIFT_CUT_OMP_H
#define LMP_BOND_HARMONIC_SHIFT_CUT_OMP_H

#include "bond_harmonic_shift_cut.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class BondHarmonicShiftCutOMP : public BondHarmonicShiftCut, public ThrOMP {

 public:
  BondHarmonicShiftCutOMP(class LAMMPS *lmp);
  void compute(int, int) override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
