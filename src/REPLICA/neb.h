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
  ~NEB();
  void command(int, char **);

 private:
  int me,me_universe;          // my proc ID in world and universe
  int ireplica,nreplica;
  MPI_Comm uworld;
  MPI_Comm roots;              // MPI comm with 1 root proc from each world
  FILE *fp;
  int compressed;

  class FixNEB *fneb;
  double **all;                // PE,plen,nlen from each replica
  double *rdist;               // normalize reaction distance, 0 to 1

  void readfile(char *);
  void open(char *);
  void print_status();
};

}

#endif
#endif
