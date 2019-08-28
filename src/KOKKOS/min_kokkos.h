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

#ifndef LMP_MIN_KOKKOS_H
#define LMP_MIN_KOKKOS_H

#include "min.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class MinKokkos : public Min {
 public:
  MinKokkos(class LAMMPS *);
  virtual ~MinKokkos();
  void init();
  void setup(int flag=1);
  void setup_minimal(int);
  void run(int);
  void cleanup();
  double fnorm_sqr();
  double fnorm_inf();

  virtual void init_style() {}
  virtual void setup_style() = 0;
  virtual void reset_vectors() = 0;
  virtual int iterate(int) = 0;

  // possible return values of iterate() method
  enum{MAXITER,MAXEVAL,ETOL,FTOL,DOWNHILL,ZEROALPHA,ZEROFORCE,
       ZEROQUAD,TRSMALL,INTERROR,TIMEOUT};

 //protected: // won't compile with CUDA
  class FixMinimizeKokkos *fix_minimize_kk;  // fix that stores auxiliary data

  DAT::t_ffloat_1d xvec;            // variables for atomic dof, as 1d vector
  DAT::t_ffloat_1d fvec;            // force vector for atomic dof, as 1d vector

  double energy_force(int);
  void force_clear();
};

}

#endif

/* ERROR/WARNING messages:

*/
