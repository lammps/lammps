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

#ifndef ATOM_VEC_MOLECULAR_H
#define ATOM_VEC_MOLECULAR_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecMolecular : public AtomVec {
 public:
  AtomVecMolecular(class LAMMPS *, int, char **);
  void grow(int);
  void reset_special();
  void copy(int, int);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_one(int, double *);
  void unpack_border(int, int, double *);
  int unpack_border_one(int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, int, char **);
  int data_atom_hybrid(int, char **);
  double memory_usage();

 private:
  int *tag,*type,*mask,*image;
  double **x,**v,**f;
  int *molecule;
  int **nspecial,**special;
  int *num_bond;
  int **bond_type,**bond_atom;
  int *num_angle;
  int **angle_type;
  int **angle_atom1,**angle_atom2,**angle_atom3;
  int *num_dihedral;
  int **dihedral_type;
  int **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;
  int *num_improper;
  int **improper_type;
  int **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;
};

}

#endif
