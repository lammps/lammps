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
FixStyle(edpd/source,FixEDPDSource);
// clang-format on
#else

#ifndef LMP_FIX_EDPDSOURCE_H
#define LMP_FIX_EDPDSOURCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEDPDSource : public Fix {
 public:
  FixEDPDSource(class LAMMPS *, int, char **);
  ~FixEDPDSource();
  int setmask();
  void init();
  void post_force(int);

 protected:
  int option;
  double center[3], radius, dLx, dLy, dLz;
  double value;
};

}    // namespace LAMMPS_NS

#endif
#endif
