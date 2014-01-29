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
  virtual ~AtomVecFull() {}
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual int pack_comm_vel(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  virtual void unpack_border(int, int, double *);
  virtual void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  int data_atom_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  bigint memory_usage();

 protected:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;
  double *q;
  tagint *molecule;
  int **nspecial;
  tagint **special;
  int *num_bond;
  int **bond_type;
  tagint **bond_atom;
  int *num_angle;
  int **angle_type;
  tagint **angle_atom1,**angle_atom2,**angle_atom3;
  int *num_dihedral;
  int **dihedral_type;
  tagint **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;
  int *num_improper;
  int **improper_type;
  tagint **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

U: Invalid atom ID in Atoms section of data file

Atom IDs must be positive integers.

*/
