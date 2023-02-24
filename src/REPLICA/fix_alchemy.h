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
FixStyle(alchemy,FixAlchemy);
// clang-format on
#else

#ifndef LMP_FIX_ALCHEMY_H
#define LMP_FIX_ALCHEMY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAlchemy : public Fix {
 public:
  FixAlchemy(class LAMMPS *, int, char **);
  ~FixAlchemy() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_integrate() override;
  void post_force(int) override;
  double compute_vector(int) override;

 protected:
  MPI_Comm samerank;
  MPI_Comm rankzero;
  double *coordbuf;
  double lambda;     // changes from 0 to 1 during run
  double epot[2];    // last (unscaled) potential energy from each replica
};
}    // namespace LAMMPS_NS

#endif
#endif
