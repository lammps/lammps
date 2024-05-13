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
FixStyle(print,FixPrint);
// clang-format on
#else

#ifndef LMP_FIX_PRINT_H
#define LMP_FIX_PRINT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPrint : public Fix {
 public:
  FixPrint(class LAMMPS *, int, char **);
  ~FixPrint() override;
  void init() override;
  void setup(int) override;
  int setmask() override;
  void end_of_step() override;

 private:
  int screenflag;
  FILE *fp;
  char *text, *copy, *work;
  int maxcopy, maxwork;
  char *var_print;
  int ivar_print;
  bigint next_print;
};

}    // namespace LAMMPS_NS

#endif
#endif
