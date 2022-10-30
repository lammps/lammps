/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Pair zero is a dummy pair interaction useful for requiring a
   force cutoff distance in the absence of pair-interactions or
   with hybrid/overlay if a larger force cutoff distance is required.

   This can be used in conjunction with bond/create to create bonds
   that are longer than the cutoff of a given force field, or to
   calculate radial distribution functions for models without
   pair interactions.

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(python/move,FixPythonMove);
// clang-format on
#else

#ifndef LMP_FIX_PYTHON_MOVE_H
#define LMP_FIX_PYTHON_MOVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPythonMove : public Fix {
 public:
  FixPythonMove(LAMMPS *lmp, int narg, char **arg);
  ~FixPythonMove() override;

  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
  void reset_dt() override;

 protected:
  void *py_move;
};

}    // namespace LAMMPS_NS

#endif
#endif
