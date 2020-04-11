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

CommandStyle(pair_species_coeff,PairSpeciesCoeff)

#else

#ifndef LMP_PAIR_SPECIES_COEFF_H
#define LMP_PAIR_SPECIES_COEFF_H

#include "pointers.h"

namespace LAMMPS_NS {

class PairSpeciesCoeff : protected Pointers {
 public:
  PairSpeciesCoeff(class LAMMPS *lmp): Pointers(lmp) {};
  void command(int, char **);

 private:
};

}

#endif
#endif

//FIXFIXFIX below @@@@@@@@@@@@@@@@@@@@@@@

/* ERROR/WARNING messages:

E: Displace_atoms command before simulation box is defined

The displace_atoms command cannot be used before a read_data,
read_restart, or create_box command.

*/
