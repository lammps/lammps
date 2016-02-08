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

FixStyle(deform/kk,FixDeformKokkos)

#else

#ifndef LMP_FIX_DEFORM_KOKKOS_H
#define LMP_FIX_DEFORM_KOKKOS_H

#include "fix_deform.h"

namespace LAMMPS_NS {

class FixDeformKokkos : public FixDeform {
 public:

  FixDeformKokkos(class LAMMPS *, int, char **);
  virtual ~FixDeformKokkos() {}
  void pre_exchange();
  void end_of_step();

 private:
  class DomainKokkos *domainKK;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot (yet) use rigid bodies with fix deform and Kokkos

UNDOCUMENTED

U: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

U: Fix deform tilt factors require triclinic box

Cannot deform the tilt factors of a simulation box unless it
is a triclinic (non-orthogonal) box.

U: Cannot use fix deform on a shrink-wrapped boundary

The x, y, z options cannot be applied to shrink-wrapped
dimensions.

U: Cannot use fix deform tilt on a shrink-wrapped 2nd dim

This is because the shrink-wrapping will change the value
of the strain implied by the tilt factor.

U: Fix deform volume setting is invalid

Cannot use volume style unless other dimensions are being controlled.

U: More than one fix deform

Only one fix deform can be defined at a time.

U: Variable name for fix deform does not exist

Self-explantory.

U: Variable for fix deform is invalid style

The variable must be an equal-style variable.

U: Final box dimension due to fix deform is < 0.0

Self-explanatory.

U: Cannot use fix deform trate on a box with zero tilt

The trate style alters the current strain.

U: Fix deform cannot use yz variable with xy

The yz setting cannot be a variable if xy deformation is also
specified.  This is because LAMMPS cannot determine if the yz setting
will induce a box flip which would be invalid if xy is also changing.

U: Fix deform is changing yz too much with xy

When both yz and xy are changing, it induces changes in xz if the
box must flip from one tilt extreme to another.  Thus it is not
allowed for yz to grow so much that a flip is induced.

*/
