// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MIN_LSRCH_KOKKOS_H
#define LMP_MIN_LSRCH_KOKKOS_H

#include "min_kokkos.h"

namespace LAMMPS_NS {

  struct s_double2 {
    double d0, d1;
    KOKKOS_INLINE_FUNCTION
    s_double2() {
      d0 = d1 = 0.0;
    }
    KOKKOS_INLINE_FUNCTION
    s_double2& operator+=(const s_double2 &rhs) {
      d0 += rhs.d0;
      d1 += rhs.d1;
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    void operator+=(const volatile s_double2 &rhs) volatile {
      d0 += rhs.d0;
      d1 += rhs.d1;
    }
  };
  //typedef s_double2 double2;

class MinLineSearchKokkos : public MinKokkos {
 public:
  MinLineSearchKokkos(class LAMMPS *);
  ~MinLineSearchKokkos();
  void init();
  void setup_style();
  void reset_vectors();

 //protected: // won't compile with CUDA
  // vectors needed by linesearch minimizers
  // allocated and stored by fix_minimize
  // x,f are stored by parent or Atom class or Pair class

  DAT::t_ffloat_1d x0;   // coords at start of linesearch
  DAT::t_ffloat_1d g;    // old gradient vector
  DAT::t_ffloat_1d h;    // search direction vector

  double *gextra;             // g,h for extra global dof, x0 is stored by fix
  double *hextra;

  typedef int (MinLineSearchKokkos::*FnPtr)(double, double &);
  FnPtr linemin;
  int linemin_quadratic(double, double &);

  double alpha_step(double, int);
  double compute_dir_deriv(double &);
};

}

#endif
