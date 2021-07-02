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

FixStyle(CAC/ALLOCVECTOR,FixCACAllocVector)

#else

#ifndef LMP_FIX_CAC_ALLOCVECTOR_H
#define LMP_FIX_CAC_ALLOCVECTOR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCACAllocVector : public Fix {

 public:
  FixCACAllocVector(class LAMMPS *, int, char **);
  ~FixCACAllocVector();
  int setmask();
  void init() {}

  double memory_usage();
  void grow_arrays(int);
  void shrink_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  virtual void clear_arrays(int, size_t);

  void add_vector(int);
  void allocate_element(int, int, int);
  double ****request_vector(int);

 private:
  int nvector, alloc_counter, callback;
  int *peratom;
  double *****nodal_vectors;
};

}

#endif
#endif
/* ERROR/WARNING messages:

*/
