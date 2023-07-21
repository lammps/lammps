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
KSpaceStyle(pppm/dielectric,PPPMDielectric);
// clang-format on
#else

#ifndef LMP_PPPM_DIELECTRIC_H
#define LMP_PPPM_DIELECTRIC_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMDielectric : public PPPM {
 public:
  PPPMDielectric(class LAMMPS *);
  ~PPPMDielectric() override;
  void compute(int, int) override;

  double **efield;
  double *phi;
  int potflag;    // 1/0 if per-atom electrostatic potential phi is needed

 protected:
  void slabcorr() override;
  void make_rho() override;
  void fieldforce_ik() override;
  void fieldforce_ad() override;
  void qsum_qsq(int warning_flag = 1) override;

  class AtomVecDielectric *avec;
  bool use_qscaled;

  void compute_ave_epsilon();
  double epsilon_ave;
};

}    // namespace LAMMPS_NS

#endif
#endif
