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
FixStyle(smd,FixSMD);
// clang-format on
#else

#ifndef LMP_FIX_SMD_H
#define LMP_FIX_SMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSMD : public Fix {
 public:
  FixSMD(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  double compute_vector(int) override;

  void write_restart(FILE *) override;
  void restart(char *) override;

 private:
  double xc, yc, zc, xn, yn, zn, r0;
  double k_smd, f_smd, v_smd;
  int xflag, yflag, zflag;
  int styleflag;
  double r_old, r_now, pmf;

  int igroup2, group2bit;
  double masstotal, masstotal2;
  int ilevel_respa;
  double ftotal[3], ftotal_all[7];
  int force_flag;

  void smd_tether();
  void smd_couple();
};

}    // namespace LAMMPS_NS

#endif
#endif
