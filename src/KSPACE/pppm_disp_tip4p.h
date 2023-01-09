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
KSpaceStyle(pppm/disp/tip4p,PPPMDispTIP4P);
// clang-format on
#else

#ifndef LMP_PPPM_DISP_TIP4P_H
#define LMP_PPPM_DISP_TIP4P_H

#include "pppm_disp.h"

namespace LAMMPS_NS {

class PPPMDispTIP4P : public PPPMDisp {
 public:
  PPPMDispTIP4P(class LAMMPS *);

  void init() override;

 protected:
  void particle_map_c(double, double, double, double, int **, int, int, int, int, int, int, int,
                      int) override;
  void make_rho_c() override;
  void fieldforce_c_ik() override;
  void fieldforce_c_ad() override;
  void fieldforce_c_peratom() override;
  void slabcorr(int) override;

 private:
  void find_M(int, int &, int &, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
