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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(PRESSURE/GREM,ComputePressureGrem);
// clang-format on
#else

#ifndef LMP_COMPUTE_PRESSURE_GREM_H
#define LMP_COMPUTE_PRESSURE_GREM_H

#include "compute_pressure.h"

namespace LAMMPS_NS {

class ComputePressureGrem : public ComputePressure {
 public:
  ComputePressureGrem(class LAMMPS *, int, char **);
  virtual ~ComputePressureGrem();
  virtual void init();
  virtual double compute_scalar();
  virtual void compute_vector();

 protected:
  // Access to gREM fix scale factor
  char *fix_grem;
  double *scale_grem;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute pressure must use group all

Virial contributions computed by potentials (pair, bond, etc) are
computed on all atoms.

E: Could not find compute pressure temperature ID

The compute ID for calculating temperature does not exist.

E: Compute pressure temperature ID does not compute temperature

The compute ID assigned to a pressure computation must compute
temperature.

E: Compute pressure requires temperature ID to include kinetic energy

The keflag cannot be used unless a temperature compute is provided.

E: Virial was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied the virial, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

E: Must use 'kspace_modify pressure/scalar no' for tensor components with kspace_style msm

Otherwise MSM will compute only a scalar pressure.  See the kspace_modify
command for details on this setting.

E: Fix grem ID for compute PRESSURE/GREM does not exist

Compute PRESSURE/GREM was passed an invalid fix id

E: Cannot extract gREM scale factor from fix grem

The fix id passed to compute PRESSURE/GREM refers to an incompatible fix

*/
