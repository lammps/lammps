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

#ifndef ATOM_VEC_GRANULAR_H
#define ATOM_VEC_GRANULAR_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecGranular : public AtomVec {
 public:
  AtomVecGranular(class LAMMPS *, int, char **);
  ~AtomVecGranular() {}
  void grow(int);
  void reset_ptrs();
  void zero_owned(int);
  void zero_ghost(int, int);
  void copy(int, int);
  int pack_comm(int, int *, double *, int, double *);
  int pack_comm_one(int, double *);
  void unpack_comm(int, int, double *);
  int unpack_comm_one(int, double *);
  int pack_reverse(int, int, double *);
  int pack_reverse_one(int, double *);
  void unpack_reverse(int, int *, double *);
  int unpack_reverse_one(int, double *);
  int pack_border(int, int *, double *, int, double *);
  int pack_border_one(int, double *);
  void unpack_border(int, int, double *);
  int unpack_border_one(int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int size_restart_one(int);
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *, int);
  void data_atom(double *, int, char **, int);
  void data_vel(int, char *, int);
  int memory_usage();

 private:
  double PI;
  int *tag,*type,*mask,*image;
  double **x,**v,**f;
  double *radius,*density,*rmass;
  double **phix,**phiv,**phia;
};

}

#endif
