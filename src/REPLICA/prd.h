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
CommandStyle(prd,PRD);
// clang-format on
#else

#ifndef LMP_PRD_H
#define LMP_PRD_H

#include "command.h"

namespace LAMMPS_NS {

class PRD : public Command {
 public:
  PRD(class LAMMPS *);

  void command(int, char **) override;

 private:
  int me, nprocs;
  int t_event, n_dephase, t_dephase, t_corr;
  double etol, ftol, temp_dephase;
  int maxiter, maxeval, temp_flag, stepmode, cmode;
  char *loop_setting, *dist_setting;

  int equal_size_replicas, natoms;
  int neigh_every, neigh_delay, neigh_dist_check;
  int quench_reneighbor;
  bigint nbuild, ndanger;

  double time_dephase, time_dynamics, time_quench, time_comm, time_output;
  double time_start;

  MPI_Comm comm_replica;
  int *counts, *displacements;
  tagint *tagall;
  double **xall;
  imageint *imageall;

  int ncoincident;

  class RanPark *random_select, *random_clock;
  class RanMars *random_dephase;
  class Compute *compute_event;
  class FixEventPRD *fix_event;
  class Velocity *velocity;
  class Compute *temperature;
  class Finish *finish;

  void dephase();
  void dynamics(int, double &);
  void quench();
  int check_event(int replica = -1);
  void share_event(int, int, int);
  void log_event();
  void replicate(int);
  void options(int, char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
