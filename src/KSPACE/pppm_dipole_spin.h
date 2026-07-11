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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/dipole/spin,PPPMDipoleSpin);
// clang-format on
#else

#ifndef LMP_PPPM_DIPOLE_SPIN_H
#define LMP_PPPM_DIPOLE_SPIN_H

#include "pppm_dipole.h"

namespace LAMMPS_NS {

class PPPMDipoleSpin : public PPPMDipole {
 public:
  PPPMDipoleSpin(class LAMMPS *);
  ~PPPMDipoleSpin() override;
  void init() override;
  void compute(int, int) override;

 protected:
  double hbar;            // reduced Planck's constant
  double mub;             // Bohr's magneton
  double mu_0;            // vacuum permeability
  double mub2mu0;         // prefactor for mech force
  double mub2mu0hbinv;    // prefactor for mag force

  void slabcorr() override;

  // spin

  void make_rho_spin();
  void fieldforce_ik_spin();
  void fieldforce_peratom_spin();
  void spsum_spsq();
};

}    // namespace LAMMPS_NS

#endif
#endif
