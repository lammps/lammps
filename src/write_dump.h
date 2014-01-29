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

CommandStyle(write_dump,WriteDump)

#else

#ifndef LMP_WRITE_DUMP_H
#define LMP_WRITE_DUMP_H

#include "pointers.h"

namespace LAMMPS_NS {

class WriteDump : protected Pointers {
 public:
  WriteDump(class LAMMPS *lmp) : Pointers(lmp) {};
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Invalid dump style

UNDOCUMENTED

U: Too many total bits for bitmapped lookup table

Table size specified via pair_modify command is too large.  Note that
a value of N generates a 2^N size table.

U: Cannot have both pair_modify shift and tail set to yes

These 2 options are contradictory.

U: Cannot use pair tail corrections with 2d simulations

The correction factors are only currently defined for 3d systems.

U: Using pair tail corrections with nonperiodic system

This is probably a bogus thing to do, since tail corrections are
computed by integrating the density of a periodic system out to
infinity.

U: Using a manybody potential with bonds/angles/dihedrals and special_bond exclusions

UNDOCUMENTED

U: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

U: Pair style requres a KSpace style

UNDOCUMENTED

U: Pair style does not support pair_write

The pair style does not have a single() function, so it can
not be invoked by pair write.

U: Invalid atom types in pair_write command

Atom types must range from 1 to Ntypes inclusive.

U: Invalid style in pair_write command

Self-explanatory.  Check the input script.

U: Invalid cutoffs in pair_write command

Inner cutoff must be larger than 0.0 and less than outer cutoff.

U: Cannot open pair_write file

The specified output file for pair energies and forces cannot be
opened.  Check that the path and name are correct.

U: Bitmapped lookup tables require int/float be same size

Cannot use pair tables on this machine, because of word sizes.  Use
the pair_modify command with table 0 instead.

U: Table inner cutoff >= outer cutoff

You specified an inner cutoff for a Coulombic table that is longer
than the global cutoff.  Probably not what you wanted.

U: Too many exponent bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

U: Too many mantissa bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

U: Too few bits for lookup table

Table size specified via pair_modify command does not work with your
machine's floating point representation.

*/
