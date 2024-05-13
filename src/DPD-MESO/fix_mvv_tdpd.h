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
FixStyle(mvv/tdpd,FixMvvTDPD);
// clang-format on
#else

#ifndef LMP_FIX_MVV_TDPD_H
#define LMP_FIX_MVV_TDPD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMvvTDPD : public Fix {
 public:
  FixMvvTDPD(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;

 protected:
  double dtv, dtf;
  double verlet;
  int cc_species;
};

}    // namespace LAMMPS_NS

#endif
#endif
