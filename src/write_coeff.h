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

CommandStyle(write_coeff,WriteCoeff)

#else

#ifndef LMP_WRITE_COEFF_H
#define LMP_WRITE_COEFF_H

#include <cstdio>
#include "pointers.h"

namespace LAMMPS_NS {

class WriteCoeff : protected Pointers {
 public:
  WriteCoeff(class LAMMPS *lmp) : Pointers(lmp) {};
  void command(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Write_coeff command before simulation box is defined

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open coeff file %s

The specified file cannot be opened.  Check that the path and name are
correct.

U: write_coeff command before simulation box is defined

Self-explanatory.

U: Atom count is inconsistent, cannot write data file

The sum of atoms across processors does not equal the global number
of atoms.  Probably some atoms have been lost.

*/
