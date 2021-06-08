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

/*------------------------------------------------------------------------
  Contributing Authors : Romain Vermorel (LFCR), Laurent Joly (ULyon)
  --------------------------------------------------------------------------*/

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(stress/mop,ComputeStressMop);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRESS_MOP_H
#define LMP_COMPUTE_STRESS_MOP_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressMop : public Compute {
 public:
  ComputeStressMop(class LAMMPS *, int, char **);
  virtual ~ComputeStressMop();
  void init();
  void init_list(int, class NeighList *);
  void compute_vector();

 private:
  void compute_pairs();

  int me, nvalues, dir;
  int *which;

  double *values_local, *values_global;
  double pos, pos1, dt, nktv2p, ftm2v;
  double area;
  class NeighList *list;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

 E: Illegal ... command

 Self-explanatory.  Check the input script syntax and compare to the
 documentation for the command.  You can use -echo screen as a
 command-line option when running LAMMPS to see the offending line.

 E: Compute stress/mop incompatible with simulation dimension

 Compute stress/mop only works with 3D simulations.

 E: Compute stress/mop incompatible with triclinic simulation box

 Self-explanatory.

 E: Compute stress/mop requires a fixed simulation box

 Compute stress/mop is not compatible with any change of volume or shape
 or boundary conditions of the simulation box.

 E: No pair style is defined for compute stress/mop

 Self-explanatory. Compute stress/mop requires the definition of a pair style.

 E: Pair style does not support compute stress/mop

 The pair style does not have a single() function, so it can
 not be invoked by compute stress/mop.

 W: compute stress/mop does not account for bond potentials

 W: compute stress/mop does not account for angle potentials

 W: compute stress/mop does not account for dihedral potentials

 W: compute stress/mop does not account for improper potentials

 W: compute stress/mop does not account for kspace contributions

 Compute stress/mop only accounts for pairwise additive interactions for
 the computation of local stress tensor components.

 */
