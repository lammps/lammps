/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/disp/dielectric,PPPMDispDielectric);
// clang-format on
#else

#ifndef LMP_PPPM_DISP_DIELECTRIC_H
#define LMP_PPPM_DISP_DIELECTRIC_H

#include "pppm_disp.h"

namespace LAMMPS_NS {

class PPPMDispDielectric : public PPPMDisp {
 public:
  PPPMDispDielectric(class LAMMPS *);
  ~PPPMDispDielectric() override;
  double memory_usage() override;
  void compute(int, int) override;
  void slabcorr(int) override;

  double **efield;
  double *phi;
  int potflag;    // 1/0 if per-atom electrostatic potential phi is needed

 protected:
  void make_rho_c() override;
  void fieldforce_c_ik() override;
  void fieldforce_c_ad() override;
  void fieldforce_c_peratom() override;
  void qsum_qsq(int warning_flag = 1) override;

  class AtomVecDielectric *avec;
  bool use_qscaled;

  void compute_ave_epsilon();
  double epsilon_ave;
  int mu_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif
