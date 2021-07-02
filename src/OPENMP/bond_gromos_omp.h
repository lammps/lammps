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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(gromos/omp,BondGromosOMP);
// clang-format on
#else

#ifndef LMP_BOND_GROMOS_OMP_H
#define LMP_BOND_GROMOS_OMP_H

#include "bond_gromos.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class BondGromosOMP : public BondGromos, public ThrOMP {

 public:
  BondGromosOMP(class LAMMPS *lmp);
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
