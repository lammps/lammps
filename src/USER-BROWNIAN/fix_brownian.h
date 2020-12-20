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

FixStyle(brownian,FixBrownian)

#else

#ifndef LMP_FIX_BROWNIAN_H
#define LMP_FIX_BROWNIAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBrownian : public Fix {
 public:
  FixBrownian(class LAMMPS *, int, char **);
  virtual ~FixBrownian();
  void init();
  void initial_integrate(int);
  void setup(int);
  void post_force(int);
  int setmask();
  void reset_dt();

 private:
  int seed;               // RNG seed
  double dt, sqrtdt;      // time step interval and its sqrt


  double gamma_t;          // translational damping param
  double diff_t;           // translational diffusion coeff

  double g1,g2;            // prefactors in time stepping
  int noise_flag;          // 0/1 for noise off/on
  int gaussian_noise_flag; // 0/1 for uniform/gaussian noise
  
protected:
  class RanMars *random;
  typedef double (RanMars::*rng_member)();
  rng_member rng_func;    // placeholder for RNG function

};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal fix brownian command.

Wrong number/type of input arguments.

E: Fix brownian viscous drag coefficient must be > 0.

Self-explanatory.

E: Fix brownian diffusion coefficient must be > 0.

Self-explanatory.

E: Fix brownian seed must be > 0.

Self-explanatory.

*/
