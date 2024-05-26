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
FixStyle(heat/flow,FixHeatFlow);
// clang-format on
#else

#ifndef LMP_FIX_HEAT_FLOW_H
#define LMP_FIX_HEAT_FLOW_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHeatFlow : public Fix {
 public:
  FixHeatFlow(class LAMMPS *, int, char **);

  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void final_integrate() override;
  void final_integrate_respa(int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  void reset_dt() override;

 protected:
  double dt;
  double cp, *cp_type;
  int cp_style;
  int first_flag;

  double calc_cp(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
