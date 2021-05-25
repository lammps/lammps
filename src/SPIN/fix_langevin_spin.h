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

#ifdef FIX_CLASS
// clang-format off
FixStyle(langevin/spin,FixLangevinSpin);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_SPIN_H
#define LMP_FIX_LANGEVIN_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevinSpin : public Fix {
 public:
  int tdamp_flag, temp_flag;    // damping and temperature flags

  FixLangevinSpin(class LAMMPS *, int, char **);
  virtual ~FixLangevinSpin();
  int setmask();
  void init();
  void setup(int);
  void add_tdamping(double *, double *);    // add transverse damping
  void add_temperature(double[3]);
  void compute_single_langevin(int, double *, double *);

 protected:
  double alpha_t;       // transverse mag. damping
  double dts;           // magnetic timestep
  double temp;          // spin bath temperature
  double D, sigma;      // bath intensity var.
  double gil_factor;    // gilbert's prefactor

  int nlevels_respa;
  class RanMars *random;
  int seed;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal langevin/spin command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix langevin period must be > 0.0

The time window for temperature relaxation must be > 0

W: Energy tally does not account for 'zero yes'

The energy removed by using the 'zero yes' flag is not accounted
for in the energy tally and thus energy conservation cannot be
monitored in this case.


E: Variable for fix langevin is invalid style

It must be an equal-style variable.


E: Cannot zero Langevin force of 0 atoms

The group has zero atoms, so you cannot request its force
be zeroed.

E: Fix langevin variable returned negative temperature

Self-explanatory.

E: Could not find fix_modify temperature ID

The compute ID for computing temperature does not exist.

E: Fix_modify temperature ID does not compute temperature

The compute ID assigned to the fix must compute temperature.

W: Group for fix_modify temp != fix group

The fix_modify command is specifying a temperature computation that
computes a temperature on a different group of atoms than the fix
itself operates on.  This is probably not what you want to do.

*/
