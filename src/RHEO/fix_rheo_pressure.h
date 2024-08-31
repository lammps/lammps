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
FixStyle(rheo/pressure,FixRHEOPressure)
// clang-format on
#else

#ifndef LMP_FIX_RHEO_PRESSURE_H
#define LMP_FIX_RHEO_PRESSURE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRHEOPressure : public Fix {
 public:
  FixRHEOPressure(class LAMMPS *, int, char **);
  ~FixRHEOPressure() override;
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double calc_pressure(double, int);
  double calc_rho(double, int);

 private:
  double *c_cubic, *csq, *csqinv, *rho0, *rho0inv, *tpower, *pbackground;
  int *pressure_style;

  class FixRHEO *fix_rheo;
};

}    // namespace LAMMPS_NS

#endif
#endif
