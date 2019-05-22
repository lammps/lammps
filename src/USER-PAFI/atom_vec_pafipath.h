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
/* ------------------------------------------------------------------------
   Contributing authors: Thomas Swinburne (CNRS & CINaM, Marseille, France)

   Please cite the related publication:
   T.D. Swinburne and M.-C. Marinica, Unsupervised calculation of free energy barriers in large crystalline systems, Physical Review Letters 2018
------------------------------------------------------------------------- */


#ifdef ATOM_CLASS

AtomStyle(pafipath,AtomVecPAFIPATH)

#else

#ifndef LMP_ATOM_VEC_PAFIPATH_H
#define LMP_ATOM_VEC_PAFIPATH_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecPAFIPATH : public AtomVec {
 public:
  AtomVecPAFIPATH(class LAMMPS *);
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  int pack_border_hybrid(int, int *, double *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int unpack_border_hybrid(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  virtual void data_atom(double *, imageint, char **);
  virtual int data_atom_hybrid(int, char **);
  virtual void pack_data(double **);
  virtual int pack_data_hybrid(int, double *);
  virtual void write_data(FILE *, int, double **);
  virtual int write_data_hybrid(FILE *, double *);
  bigint memory_usage();

 private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  double **x,**v,**f;		// lattice quantities
  // 0th, 1st and 2nd derivative of reference path w.r.t. to path coordinate r
  double **path,**norm,**dnorm;
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

*/
