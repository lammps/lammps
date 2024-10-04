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

#ifdef FIX_CLASS
// clang-format off
FixStyle(qeq/ctip,FixQEqCTIP);
// clang-format on
#else

#ifndef LMP_FIX_QEQ_CTIP_H
#define LMP_FIX_QEQ_CTIP_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqCTIP : public FixQEq {
 public:
  FixQEqCTIP(class LAMMPS *, int, char **);
  ~FixQEqCTIP() override;

  void init() override;
  void pre_force(int) override;

 protected:
  void init_matvec();
  void sparse_matvec(sparse_matrix *, double *, double *) override;
  void compute_H();
  void extract_ctip();
  int calculate_check_Q();
  double *reff, *reffsq, *reff4, *reff7, *s2d_self;
  double **shield, **shieldcu, **reffc, **reffcsq, **reffc4, **reffc7;
  double **s2d_shift, **f_shift, **e_shift;

  double cdamp;
  int maxrepeat;
  int nout;
};
}    // namespace LAMMPS_NS
#endif
#endif
