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
/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods (SELMs) Package

  Paul J. Atzberger
  http://atzberger.org/

  Please cite the follow paper when referencing this package

  "Fluctuating Hydrodynamics Methods for Dynamic Coarse-Grained Implicit-Solvent Simulations in LAMMPS,"
  Wang, Y. and Sigurdsson, J. K. and Atzberger, P. J., SIAM Journal on Scientific Computing, 38(5), 2016.

  @article{atz_selm_lammps_fluct_hydro,
    title = {Fluctuating Hydrodynamics Methods for Dynamic
    Coarse-Grained Implicit-Solvent Simulations in LAMMPS},
    author = {Wang, Y. and Sigurdsson, J. K. and Atzberger, P. J.},
    journal = {SIAM Journal on Scientific Computing},
    volume = {38},
    number = {5},
    pages = {S62-S77},
    year = {2016},
    doi = {10.1137/15M1026390},
    URL = {https://doi.org/10.1137/15M1026390},
  }

  For latest version of the codes, examples, and additional information see
  http://mango-selm.org/

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(selm,FixSELM);
// clang-format on
#else

#ifndef LMP_FIX_SELM_H
#define LMP_FIX_SELM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSELM : public Fix {

 public:

  FixSELM(class LAMMPS *, int, char **);
  virtual ~FixSELM();

  int setmask();
  virtual void init();
  void setup(int vflag);
  virtual void initial_integrate(int);
  virtual void final_integrate();
  void reset_dt();
  void post_force(int vflag);

  void pre_exchange();
  void end_of_step();
 protected:
  class DriverSELM *driver_selm;
  static int instances;
};
}    // namespace LAMMPS_NS
#endif
#endif
