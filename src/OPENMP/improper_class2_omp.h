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

#ifdef IMPROPER_CLASS
// clang-format off
ImproperStyle(class2/omp,ImproperClass2OMP);
// clang-format on
#else

#ifndef LMP_IMPROPER_CLASS2_OMP_H
#define LMP_IMPROPER_CLASS2_OMP_H

#include "improper_class2.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class ImproperClass2OMP : public ImproperClass2, public ThrOMP {

 public:
  ImproperClass2OMP(class LAMMPS *lmp);
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval(int ifrom, int ito, ThrData *const thr);

  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void angleangle_thr(int, int, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
