/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nvt/body,FixNVTBody)

#else

#ifndef LMP_FIX_NVT_BODY_H
#define LMP_FIX_NVT_BODY_H

#include "fix_nh_body.h"

namespace LAMMPS_NS {

class FixNVTBody : public FixNHBody {
 public:
  FixNVTBody(class LAMMPS *, int, char **);
  ~FixNVTBody() {}
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control must be used with fix nvt/body

Self-explanatory.

E: Pressure control can not be used with fix nvt/body

Self-explanatory.

*/
