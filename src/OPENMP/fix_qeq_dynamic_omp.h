/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/dynamic/omp,FixQEqDynamicOMP);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_DYNAMIC_OMP_H
#define LMP_FIX_QEQ_DYNAMIC_OMP_H

#include "fix_qeq_dynamic.h"

namespace LAMMPS_NS {

class FixQEqDynamicOMP : public FixQEqDynamic {
 public:
  FixQEqDynamicOMP(class LAMMPS *, int, char **);
  ~FixQEqDynamicOMP() override;

 protected:
  double compute_eneg() override;

  double **b_temp;
  int nmax_btmp;
};

}    // namespace LAMMPS_NS

#endif
#endif
