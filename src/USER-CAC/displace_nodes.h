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

CommandStyle(displace_nodes,DisplaceNodes)

#else

#ifndef LMP_DISPLACE_NODES_H
#define LMP_DISPLACE_NODES_H

#include "command.h"

namespace LAMMPS_NS {

class DisplaceNodes : public Command {
 public:
  DisplaceNodes(class LAMMPS *);
  ~DisplaceNodes();
  void command(int, char **);

 private:
  int igroup,groupbit;
  int scaleflag;
  double *mvec;

  void move(int, char *, double);
  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: displace_nodes command before simulation box is defined

The displace_nodes command cannot be used before a read_data,
read_restart, or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot displace_nodes after reading restart file with per-atom info

This is because the restart file info cannot be migrated with the
atoms.  You can get around this by performing a 0-timestep run which
will assign the restart file info to actual atoms.

E: Could not find displace_nodes group ID

Group ID used in the displace_nodes command does not exist.

W: Attempting to displace atoms in rigid bodies

UNDOCUMENTED

E: Invalid displace_nodes rotate axis for 2d

Axis must be in z direction.

E: Zero length rotation vector with displace_nodes

Self-explanatory.

W: Lost atoms via displace_nodes: original %ld current %ld

The command options you have used caused atoms to be lost.

E: Variable name for displace_nodes does not exist

Self-explanatory.

E: Variable for displace_nodes is invalid style

It must be an equal-style or atom-style variable.

*/