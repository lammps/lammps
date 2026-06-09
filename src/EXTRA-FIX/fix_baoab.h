/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(baoab,FixBAOAB)
// clang-format on
#else

#ifndef LMP_FIX_BAOAB_H
#define LMP_FIX_BAOAB_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBAOAB : public Fix {
 public:
  FixBAOAB(class LAMMPS *, int, char **);
  ~FixBAOAB() override;

  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  double compute_scalar() override;

 private:
  double t_start, t_stop;    // temperature ramp endpoints
  double t_target;           // current target temperature
  double gamma;              // friction coefficient = 1/damp
  double damp;               // damping time (user-specified)
  int seed;                  // RNG seed

  // Precomputed timestep-dependent quantities
  double dtf;      // 0.5*dt*ftm2v  (force->velocity kick)
  double dtby2;    // 0.5*dt        (half-step drift)
  double c1;       // exp(-gamma*dt) for O step

  // Energy accounting
  double energy;            // cumulative thermostat energy
  double energy_onestep;    // thermostat energy this step

  // Options
  int zeroflag;    // zero total random momentum?

  class RanMars *random;

  void compute_target();         // update t_target from ramp
  void recompute_dt_coeffs();    // recompute dtf, dtby2, c1
};

}    // namespace LAMMPS_NS
#endif
#endif
