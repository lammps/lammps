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

AtomStyle(template,AtomVecTemplate)

#else

#ifndef LMP_ATOM_VEC_TEMPLATE_H
#define LMP_ATOM_VEC_TEMPLATE_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecTemplate : public AtomVec {
 public:
  AtomVecTemplate(class LAMMPS *);
  virtual ~AtomVecTemplate() {}
  void process_args(int, char **);
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
  void data_atom(double *, tagint, char **);
  int data_atom_hybrid(int, char **);
  void pack_data(double **);
  int pack_data_hybrid(int, double *);
  void write_data(FILE *, int, double **);
  int write_data_hybrid(FILE *, double *);
  bigint memory_usage();

 protected:
  tagint *tag;
  int *type,*mask;
  tagint *image;
  double **x,**v,**f;
  tagint *molecule;
  int *molindex,*molatom;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Molecule template ID for atom_style template does not exist

Self-explanatory.

E: Atom style template molecule must have atom types

The defined molecule(s) does not specify atom types.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom ID in Atoms section of data file

Atom IDs must be positive integers.

E: Invalid template index in Atoms section of data file

The template indices must be between 1 to N, where N is the number of
molecules in the template.

E: Invalid template atom in Atoms section of data file

The atom indices must be between 1 to N, where N is the number of
atoms in the template molecule the atom belongs to.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

*/
