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
FixStyle(nvt/sllod/omp,FixNVTSllodOMP);
// clang-format on
#else

#ifndef LMP_FIX_NVT_SLLOD_OMP_H
#define LMP_FIX_NVT_SLLOD_OMP_H

#include "fix_nh_omp.h"

namespace LAMMPS_NS {

class FixNVTSllodOMP : public FixNHOMP {
 public:
  FixNVTSllodOMP(class LAMMPS *, int, char **);
  void init() override;

 private:
  int nondeformbias;
  int psllod_flag;     // 0 for SLLOD, 1 for p-SLLOD
  int peculiar_flag;   // 0 for lab frame, 1 for peculiar
  int kick_flag;       // 0 for no initial velocity kick, 1 for kick
  enum {REVERSIBLE, LEGACY} integrator;

  void nh_v_temp() override;
  void nve_x() override;
  int size_restart_global() override;
  int pack_restart_data(double *list) override;
  void restart(char *buf) override;
  int modify_param(int narg, char **arg) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
