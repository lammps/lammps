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

  void nh_v_temp() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
