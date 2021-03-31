/* -*- c++ -*- ----------------------------------------------------------
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

CommandStyle(neb/spin,NEBSpin)

#else

#ifndef LMP_NEB_SPIN_H
#define LMP_NEB_SPIN_H

#include "pointers.h"

namespace LAMMPS_NS {

class NEBSpin : protected Pointers {
 public:
  NEBSpin(class LAMMPS *);
  ~NEBSpin();
  void command(int, char **);  // process neb/spin command
  void run();                  // run NEBSpin

  double ebf,ebr;              // forward and reverse energy barriers

 private:
  int me,me_universe;          // my proc ID in world and universe
  int ireplica,nreplica;
  bool verbose;
  MPI_Comm uworld;
  MPI_Comm roots;              // MPI comm with 1 root proc from each world
  FILE *fp;
  int compressed;
  double etol;                 // energy tolerance convergence criterion
  double ttol;                 // torque tolerance convergence criterion
  int n1steps, n2steps;        // number of steps in stage 1 and 2
  int nevery;                  // output interval
  char *inpfile;                // name of file containing final state

  class FixNEBSpin *fneb;
  int numall;                  // per-replica dimension of array all
  double **all;                // PE,plen,nlen,gradvnorm from each replica
  double *rdist;               // normalize reaction distance, 0 to 1
  double *freplica;            // force on an image
  double *fmaxatomInRepl;      // force on an image

  void readfile(char *, int);
  int initial_rotation(double *, double *, double);
  void open(char *);
  void print_status();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: NEBSpin command before simulation box is defined

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use NEBSpin with a single replica

Self-explanatory.

E: Cannot use NEBSpin unless atom map exists

Use the atom_modify command to create an atom map.

E: NEBSpin requires use of fix neb

Self-explanatory.

E: NEBSpin requires damped dynamics minimizer

Use a different minimization style.

E: Too many timesteps for NEBSpin

You must use a number of timesteps that fit in a 32-bit integer
for NEBSpin.

E: Too many timesteps

The cumulative timesteps must fit in a 64-bit integer.

E: Unexpected end of neb/spin file

A read operation from the file failed.

E: Incorrect atom format in neb/spin file

The number of fields per line is not what expected.

E: Invalid atom IDs in neb/spin file

An ID in the file was not found in the system.

E: Cannot open gzipped file

LAMMPS was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DLAMMPS_GZIP.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

U: Can only use NEBSpin with 1-processor replicas

This is current restriction for NEBSpin as implemented in LAMMPS.

U: Cannot use NEBSpin with atom_modify sort enabled

This is current restriction for NEBSpin implemented in LAMMPS.

*/
