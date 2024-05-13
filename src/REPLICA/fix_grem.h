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
FixStyle(grem,FixGrem);
// clang-format on
#else

#ifndef LMP_FIX_GREM_H
#define LMP_FIX_GREM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGrem : public Fix {
 public:
  FixGrem(class LAMMPS *, int, char **);
  ~FixGrem() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void *extract(const char *, int &) override;
  double compute_scalar() override;
  double scale_grem, lambda, eta, h0;
  int pressflag;

 private:
  double tbath, pressref;

 protected:
  char *id_temp, *id_press, *id_ke, *id_pe, *id_nh;
  class Compute *temperature, *pressure, *ke, *pe;
};

}    // namespace LAMMPS_NS

#endif
#endif
