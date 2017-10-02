/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS

DumpStyle(cfg/uef,DumpCFGUef)

#else

#ifndef LMP_DUMP_CFG_UEF_H
#define LMP_DUMP_CFG_UEF_H

#include "dump_cfg.h"

namespace LAMMPS_NS {

class DumpCFGUef : public DumpCFG {
 public:
  DumpCFGUef(LAMMPS *lmp, int narg, char **arg) :
    DumpCFG(lmp, narg, arg){}
  void init_style();
  void write_header(bigint);

 protected:
  int ifix_uef;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Can't use dump cfg/uef without defining a fix nvt/npt/uef

Self-explanatory.

*/
