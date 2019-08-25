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

CommandStyle(replicate,Replicate)

#else

#ifndef LMP_REPLICATE_H
#define LMP_REPLICATE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Replicate : protected Pointers {
 public:
  Replicate(class LAMMPS *);
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Replicate command before simulation box is defined

The replicate command cannot be used before a read_data, read_restart,
or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot replicate 2d simulation in z dimension

The replicate command cannot replicate a 2d simulation in the z
dimension.

W: Replicating in a non-periodic dimension

The parameters for a replicate command will cause a non-periodic
dimension to be replicated; this may cause unwanted behavior.

E: Cannot replicate with fixes that store atom quantities

Either fixes are defined that create and store atom-based vectors or a
restart file was read which included atom-based vectors for fixes.
The replicate command cannot duplicate that information for new atoms.
You should use the replicate command before fixes are applied to the
system.

E: Replicated system atom IDs are too big

See the setting for tagint in the src/lmptype.h file.

E: Replicated system is too big

See the setting for bigint in the src/lmptype.h file.

E: Replicate did not assign all atoms correctly

Atoms replicated by the replicate command were not assigned correctly
to processors.  This is likely due to some atom coordinates being
outside a non-periodic simulation box.

*/
