/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Akhlak Mahmood

   Contact:
     Department of Materials Science and Engineering,
     North Carolina State University,
     Raleigh, NC, USA

     amahmoo3@ncsu.edu; mahmoodakhlak@gmail.com
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(mspin/nvt, FixMspinNVT)
// clang-format on
#else

#ifndef LMP_FIX_MSPIN_NVT_H
#define LMP_FIX_MSPIN_NVT_H

#include "fix_mspin_nh.h"

namespace LAMMPS_NS {
class FixMspinNVT : public FixMspinNH {
 public:
  FixMspinNVT(class LAMMPS *, int, char **);
  ~FixMspinNVT() {}
};
}    // namespace LAMMPS_NS

#endif
#endif
