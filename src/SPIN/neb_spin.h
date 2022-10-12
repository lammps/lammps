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
CommandStyle(neb/spin,NEBSpin);
// clang-format on
#else

#ifndef LMP_NEB_SPIN_H
#define LMP_NEB_SPIN_H

#include "command.h"

namespace LAMMPS_NS {

class NEBSpin : public Command {
 public:
  NEBSpin(class LAMMPS *);
  ~NEBSpin() override;
  void command(int, char **) override;    // process neb/spin command
  void run();                             // run NEBSpin

  double ebf, ebr;    // forward and reverse energy barriers

 private:
  int me, me_universe;    // my proc ID in world and universe
  int ireplica, nreplica;
  bool verbose;
  MPI_Comm uworld;
  MPI_Comm roots;    // MPI comm with 1 root proc from each world
  FILE *fp;
  int compressed;
  double etol;             // energy tolerance convergence criterion
  double ttol;             // torque tolerance convergence criterion
  int n1steps, n2steps;    // number of steps in stage 1 and 2
  int nevery;              // output interval
  char *inpfile;           // name of file containing final state

  class FixNEBSpin *fneb;
  int numall;                // per-replica dimension of array all
  double **all;              // PE,plen,nlen,gradvnorm from each replica
  double *rdist;             // normalize reaction distance, 0 to 1
  double *freplica;          // force on an image
  double *fmaxatomInRepl;    // force on an image

  void readfile(char *, int);
  int initial_rotation(double *, double *, double);
  void open(char *);
  void print_status();
};

}    // namespace LAMMPS_NS

#endif
#endif
