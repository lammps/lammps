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

#ifdef FIX_CLASS

FixStyle(property/atom/kk,FixPropertyAtomKokkos)

#else

#ifndef LMP_FIX_PROPERTY_ATOM_KOKKOS_H
#define LMP_FIX_PROPERTY_ATOM_KOKKOS_H

#include "fix_property_atom.h"

namespace LAMMPS_NS {

class FixPropertyAtomKokkos : public FixPropertyAtom {
 public:
  FixPropertyAtomKokkos(class LAMMPS *, int, char **);
  virtual ~FixPropertyAtomKokkos() {}

  void grow_arrays(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix property/atom mol when atom_style already has molecule attribute

Self-explanatory.

E: Fix property/atom cannot specify mol twice

Self-explanatory.

E: Fix property/atom q when atom_style already has charge attribute

Self-explanatory.

E: Fix property/atom cannot specify q twice

Self-explanatory.

E: Fix property/atom vector name already exists

The name for an integer or floating-point vector must be unique.

W: Fix property/atom mol or charge w/out ghost communication

A model typically needs these properties defined for ghost atoms.

E: Atom style was redefined after using fix property/atom

This is not allowed.

E: Incorrect %s format in data file

A section of the data file being read by fix property/atom does
not have the correct number of values per line.

E: Too few lines in %s section of data file

Self-explanatory.

E: Invalid atom ID in %s section of data file

An atom in a section of the data file being read by fix property/atom
has an invalid atom ID that is <= 0 or > the maximum existing atom ID.

*/
