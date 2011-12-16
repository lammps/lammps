/* ----------------------------------------------------------------------
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

FixStyle(npt,FixNPT)

#else

#ifndef LMP_FIX_NPT_H
#define LMP_FIX_NPT_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNPT : public FixNH {
 public:
  FixNPT(class LAMMPS *, int, char **);
  ~FixNPT() {}
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
