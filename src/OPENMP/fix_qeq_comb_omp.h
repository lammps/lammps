/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/comb/omp,FixQEQCombOMP);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_COMB_OMP_H
#define LMP_FIX_QEQ_COMB_OMP_H

#include "fix_qeq_comb.h"

namespace LAMMPS_NS {

class FixQEQCombOMP : public FixQEQComb {
 public:
  FixQEQCombOMP(class LAMMPS *, int, char **);
  void init() override;
  void post_force(int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
