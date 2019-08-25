/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------

   This file is a part of the USER-MANIFOLD package.

   Copyright (2013-2015) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterised by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   This fix subtracts force components that point out of the manifold,
   useful for minimisation on surfaces.


------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(manifoldforce,FixManifoldForce)

#else

#ifndef LMP_FIX_MANIFOLDFORCE_H
#define LMP_FIX_MANIFOLDFORCE_H

#include "fix.h"
#include "manifold.h"

namespace LAMMPS_NS {

class FixManifoldForce : public Fix {
 public:
  FixManifoldForce(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);


 private:
  user_manifold::manifold *ptr_m;

  // Stuff to store the parameters in.
  int nvars; // # of args after manifold name.
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
