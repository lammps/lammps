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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(hyper,Hyper);
// clang-format on
#else

#ifndef LMP_HYPER_H
#define LMP_HYPER_H

#include "command.h"

namespace LAMMPS_NS {

class Hyper : public Command {
 public:
  Hyper(class LAMMPS *);

  void command(int, char **) override;

 private:
  int me, nprocs;
  int t_event;
  double etol, ftol;
  int maxiter, maxeval;
  int stepmode, dumpflag, rebond;
  std::vector<class Dump *> dumplist;

  int neigh_every, neigh_delay, neigh_dist_check;
  int quench_reneighbor;
  bigint nbuild, ndanger;

  double time_dynamics, time_quench;
  double time_start;

  class FixHyper *fix_hyper;
  class FixEventHyper *fix_event;
  class ComputeEventDisplace *compute_event;
  class Finish *finish;

  void dynamics(int, double &);
  void quench(int flag);
  void options(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
