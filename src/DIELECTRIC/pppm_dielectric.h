/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

  void fieldforce_ik() override;
  void fieldforce_ad() override;

  class AtomVecDielectric *avec;
};

}    // namespace LAMMPS_NS

#endif
#endif
