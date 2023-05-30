/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nvt/uef,FixNVTUef);
// clang-format on
#else

#ifndef LMP_FIX_NVT_UEF_H
#define LMP_FIX_NVT_UEF_H

#include "fix_nh_uef.h"

namespace LAMMPS_NS {

class FixNVTUef : public FixNHUef {
 public:
  FixNVTUef(class LAMMPS *, int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
