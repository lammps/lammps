/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/ctip/omp,FixQEqCTIPOMP);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_CTIP_OMP_H
#define LMP_FIX_QEQ_CTIP_OMP_H

#include "fix_qeq_ctip.h"

namespace LAMMPS_NS {

class FixQEqCTIPOMP : public FixQEqCTIP {
 public:
  FixQEqCTIPOMP(class LAMMPS *, int, char **);
  ~FixQEqCTIPOMP() override;
  void pre_force(int) override;

 protected:
  void init_matvec_thr();
  void compute_H_thr();
  void sparse_matvec(sparse_matrix *, double *, double *) override;

  double **b_temp;
  int nmax_btmp;
};

}    // namespace LAMMPS_NS

#endif
#endif
