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

/* ----------------------------------------------------------------------
Contributing author: Gabriel Alkuino (Syracuse University) - gsalkuin@syr.edu
Modified from fix_efield
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(efield/lepton,FixEfieldLepton);
// clang-format on
#else

#ifndef LMP_FIX_EFIELD_LEPTON_H
#define LMP_FIX_EFIELD_LEPTON_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEfieldLepton : public Fix {

 public:
  FixEfieldLepton(class LAMMPS *, int, char **);
  ~FixEfieldLepton() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 protected:
  char *idregion;
  class Region *region;
  int ilevel_respa;
  std::string expr;

  int force_flag;
  double h = -1.0;
  double qe2f, mue2e;
  double fsum[4], fsum_all[4];
};
}    // namespace LAMMPS_NS
#endif
#endif
