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

#ifndef DELETE_ATOMS_H
#define DELETE_ATOMS_H

#include "pointers.h"

namespace LAMMPS_NS {

class DeleteAtoms : protected Pointers {
 public:
  DeleteAtoms(class LAMMPS *);
  void command(int, char **);

 private:
  int *dlist;

  void delete_group(int, char **);
  void delete_region(int, char **);
  void delete_overlap(int, char **);
  void delete_porosity(int, char **);
};

}

#endif
