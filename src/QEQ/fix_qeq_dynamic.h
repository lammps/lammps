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

FixStyle(qeq/dynamic,FixQEqDynamic)

#else

#ifndef LMP_FIX_QEQ_DYNAMIC_H
#define LMP_FIX_QEQ_DYNAMIC_H

#include "fix_qeq.h"

namespace LAMMPS_NS {

class FixQEqDynamic : public FixQEq {
 public:
  FixQEqDynamic(class LAMMPS *, int, char **);
  ~FixQEqDynamic() {}
  void init();
  void pre_force(int);

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 private:
  double compute_eneg();
};

}

#endif
#endif
