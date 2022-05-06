/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nvt/body,FixNVTBody);
// clang-format on
#else

#ifndef LMP_FIX_NVT_BODY_H
#define LMP_FIX_NVT_BODY_H

#include "fix_nh_body.h"

namespace LAMMPS_NS {

class FixNVTBody : public FixNHBody {
 public:
  FixNVTBody(class LAMMPS *, int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
