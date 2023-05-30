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
FixStyle(nphug,FixNPHug);
// clang-format on
#else

#ifndef LMP_FIX_NPHUG_H
#define LMP_FIX_NPHUG_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNPHug : public FixNH {
 public:
  FixNPHug(class LAMMPS *, int, char **);
  ~FixNPHug() override;
  void init() override;
  void setup(int) override;
  int modify_param(int, char **) override;
  int pack_restart_data(double *) override;    // pack restart data
  void restart(char *) override;

 private:
  class Compute *pe;    // PE compute pointer

  void compute_temp_target() override;
  double compute_vector(int) override;
  double compute_etotal();
  double compute_vol();
  double compute_hugoniot();
  double compute_us();
  double compute_up();

  char *id_pe;
  int peflag;
  int v0_set, p0_set, e0_set;
  double v0, p0, e0, rho0;
  int idir;
  int uniaxial;

  int size_restart_global() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
