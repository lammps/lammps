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
KSpaceStyle(pppm/tip4p,PPPMTIP4P);
// clang-format on
#else

#ifndef LMP_PPPM_TIP4P_H
#define LMP_PPPM_TIP4P_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMTIP4P : public PPPM {
 public:
  PPPMTIP4P(class LAMMPS *);
  void init() override;

 protected:
  void particle_map() override;
  void make_rho() override;
  void fieldforce_ik() override;
  void fieldforce_ad() override;
  void fieldforce_peratom() override;
  void slabcorr() override;

 private:
  void find_M(int, int &, int &, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
