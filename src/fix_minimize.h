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

FixStyle(MINIMIZE,FixMinimize)

#else

#ifndef LMP_FIX_MINIMIZE_H
#define LMP_FIX_MINIMIZE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMinimize : public Fix {
  friend class MinLineSearch;

 public:
  FixMinimize(class LAMMPS *, int, char **);
  ~FixMinimize();
  int setmask();
  void init() {}

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void add_vector(int);
  double *request_vector(int);
  void store_box();
  void reset_coords();

 private:
  int nvector;
  int *peratom;
  double **vectors;
  double boxlo[3],boxhi[3];

  void box_swap();
};

}

#endif
#endif
/* ERROR/WARNING messages:

*/
