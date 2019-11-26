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

#ifdef ATOM_CLASS

AtomStyle(full,AtomVecFull)

#else

#ifndef LMP_ATOM_VEC_FULL_H
#define LMP_ATOM_VEC_FULL_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecFull : public AtomVec {
 public:
  AtomVecFull(class LAMMPS *);
  ~AtomVecFull();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void data_atom(double *, imageint, char **);

 private:
  int bond_per_atom,angle_per_atom,dihedral_per_atom,improper_per_atom;
  int *bond_negative,*angle_negative,*dihedral_negative,*improper_negative;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
