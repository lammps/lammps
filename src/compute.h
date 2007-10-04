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

#ifndef COMPUTE_H
#define COMPUTE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Compute : protected Pointers {
 public:
  char *id,*style;
  int igroup,groupbit;

  double scalar;            // computed global scalar
  double *vector;           // computed global vector
  double *scalar_atom;      // computed per-atom scalar
  double **vector_atom;     // computed per-atom vector

  int scalar_flag;          // 0/1 if compute_scalar() function exists
  int vector_flag;          // 0/1 if compute_vector() function exists
  int peratom_flag;         // 0/1 if compute_peratom() function exists
  int size_vector;          // N = size of global vector
  int size_peratom;         // 0 = just scalar_atom, N = size of vector_atom

  int extensive;      // 0/1 if scalar,vector are intensive/extensive values
  int tempflag;       // 1 if Compute can be used as temperature
                      // must have both compute_scalar, compute_vector
  int pressflag;      // 1 if Compute can be used as pressure (uses virial)
                      // must have both compute_scalar, compute_vector
  double dof;         // degrees-of-freedom for temperature

  int npre;           // # of computes to compute before this one
  char **id_pre;      // IDs of Computes to compute before this one

  int comm_forward;     // size of forward communication (0 if none)
  int comm_reverse;     // size of reverse communication (0 if none)

  Compute(class LAMMPS *, int, char **);
  virtual ~Compute();
  void modify_params(int, char **);
  virtual void init() = 0;
  virtual void init_list(int, class NeighList *) {}
  virtual double compute_scalar() {return 0.0;}
  virtual void compute_vector() {}
  virtual void compute_peratom() {}

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual double memory_usage() {return 0.0;}

 protected:
  int extra_dof,dynamic;
};

}

#endif
