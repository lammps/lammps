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

#ifdef COMMAND_CLASS

CommandStyle(neb,NEB)

#else

#ifndef LMP_NEB_H
#define LMP_NEB_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class NEB : protected Pointers {
 public:
  NEB(class LAMMPS *);
  NEB(class LAMMPS *, double, double, int, int, int, double *, double *);
  ~NEB();
  void command(int, char **);  // process neb command
  void run();                  // run NEB

  double ebf,ebr;              // forward and reverse energy barriers

 private:
  int me,me_universe;          // my proc ID in world and universe
  int ireplica,nreplica;
  MPI_Comm uworld;
  MPI_Comm roots;              // MPI comm with 1 root proc from each world
  FILE *fp;
  int compressed;
  double etol;                 // energy tolerance convergence criterion
  double ftol;                 // force tolerance convergence criterion
  int n1steps, n2steps;        // number of steps in stage 1 and 2
  int nevery;                  // output interval
  char *infile;                // name of file containing final state

  class FixNEB *fneb;
  int nall;                    // per-replica dimension of array all
  double **all;                // PE,plen,nlen,gradvnorm from each replica
  double *rdist;               // normalize reaction distance, 0 to 1

  void readfile(char *);
  void open(char *);
  void print_status();
};

}

#endif
#endif
