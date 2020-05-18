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

  void grow_pointers();
  void pack_restart_pre(int);
  void pack_restart_post(int);
  void unpack_restart_init(int);
  void data_atom_post(int);

 private:
  int *num_bond,*num_angle,*num_dihedral,*num_improper;
  int **bond_type,**angle_type,**dihedral_type,**improper_type;
  int **nspecial;

  int any_bond_negative,any_angle_negative,
    any_dihedral_negative,any_improper_negative;
  int bond_per_atom,angle_per_atom,dihedral_per_atom,improper_per_atom;
  int *bond_negative,*angle_negative,*dihedral_negative,*improper_negative;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
