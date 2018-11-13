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

FixStyle(npt/cauchy,FixNPTCauchy)

#else

#ifndef LMP_FIX_NPT_CAUCHY_H
#define LMP_FIX_NPT_CAUCHY_H

#include "fix_cauchy.h"

namespace LAMMPS_NS {

class FixNPTCauchy : public FixCauchy {
 public:
  FixNPTCauchy(class LAMMPS *, int, char **);
  ~FixNPTCauchy() {}
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control must be used with fix npt

Self-explanatory.

E: Pressure control must be used with fix npt

Self-explanatory.

*/
