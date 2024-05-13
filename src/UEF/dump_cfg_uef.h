/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(cfg/uef,DumpCFGUef);
// clang-format on
#else

#ifndef LMP_DUMP_CFG_UEF_H
#define LMP_DUMP_CFG_UEF_H

#include "dump_cfg.h"

namespace LAMMPS_NS {

class DumpCFGUef : public DumpCFG {
 public:
  DumpCFGUef(LAMMPS *lmp, int narg, char **arg) : DumpCFG(lmp, narg, arg) {}
  void init_style() override;
  void write_header(bigint) override;

 protected:
  int ifix_uef;
};

}    // namespace LAMMPS_NS

#endif
#endif
