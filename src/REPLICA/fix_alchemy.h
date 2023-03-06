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
  double compute_scalar() override;
  double compute_vector(int) override;
  void *extract(const char *, int &) override;

 protected:
  MPI_Comm samerank;
  double *commbuf;
  class Compute *pe, *temp, *press;
  std::string id_pe, id_temp, id_press, id_lambda;
  double lambda;         // scaling prefactor for combining the partitions
  double epot[3];        // last (unscaled) potential energy from each replica and combined energy
  double pressure[6];    // joined pressure
  int progress;          // for progress indicator
  int sync_box;          // 1 of box dimensions need to be synchronized
  int nmax;
  int ivar;

  void check_consistency_atoms();
};
}    // namespace LAMMPS_NS

#endif
#endif
