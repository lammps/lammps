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

#ifdef FIX_CLASS
// clang-format off
FixStyle(evaporate,FixEvaporate);
// clang-format on
#else

#ifndef LMP_FIX_EVAPORATE_H
#define LMP_FIX_EVAPORATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEvaporate : public Fix {
 public:
  FixEvaporate(class LAMMPS *, int, char **);
  ~FixEvaporate();
  int setmask();
  void init();
  void pre_exchange();
  double compute_scalar();
  double memory_usage();

 private:
  int nevery, nflux, iregion;
  int molflag;
  int ndeleted;
  char *idregion;

  int nmax;
  int *list, *mark;

  class RanPark *random;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix evaporate does not exist

Self-explanatory.

E: Cannot evaporate atoms in atom_modify first group

This is a restriction due to the way atoms are organized in
a list to enable the atom_modify first command.

W: Fix evaporate may delete atom with non-zero molecule ID

This is probably an error, since you should not delete only one atom
of a molecule.

E: Fix evaporate molecule requires atom attribute molecule

The atom style being used does not define a molecule ID.

*/
