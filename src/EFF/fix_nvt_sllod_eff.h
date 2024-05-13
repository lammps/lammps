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
FixStyle(nvt/sllod/eff,FixNVTSllodEff);
// clang-format on
#else

#ifndef LMP_FIX_NVT_SLODD_EFF_H
#define LMP_FIX_NVT_SLODD_EFF_H

#include "fix_nh_eff.h"

namespace LAMMPS_NS {

class FixNVTSllodEff : public FixNHEff {
 public:
  FixNVTSllodEff(class LAMMPS *, int, char **);
  void init() override;

 private:
  int nondeformbias;
  int psllod_flag;

  void nh_v_temp() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
