/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   -----------------------------------------------------------------------

   This file is a part of the MANIFOLD package.

   Copyright (2013-2014) Stefan Paquay, Eindhoven University of Technology.
   License: GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   This file is part of the user-manifold package written by
   Stefan Paquay at the Eindhoven University of Technology.
   This module makes it possible to do MD with particles constrained
   to pretty arbitrary manifolds characterized by some constraint function
   g(x,y,z) = 0 and its normal grad(g). The number of manifolds available
   right now is limited but can be extended straightforwardly by making
   a new class that inherits from manifold and implements all pure virtual
   methods.

   Thanks to Remy Kusters for beta-testing!

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nvt/manifold/rattle,FixNVTManifoldRattle);
// clang-format on
#else

#ifndef LMP_FIX_NVT_MANIFOLD_RATTLE_H
#define LMP_FIX_NVT_MANIFOLD_RATTLE_H

#include "fix_nve_manifold_rattle.h"

/*
  FixNVTManifoldRattle works by wrapping some Nose-Hoover thermostat routines
  around the time integration functions of FixNVEManifoldRattle.
*/
namespace LAMMPS_NS {

// namespace user_manifold {

class FixNVTManifoldRattle : public FixNVEManifoldRattle {
 public:
  FixNVTManifoldRattle(LAMMPS *, int, char **, int = 1);
  virtual ~FixNVTManifoldRattle();

  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void init();
  virtual void reset_dt();
  virtual int setmask();
  virtual void setup(int);    // Not needed for fixNVE but is for fixNVT
  virtual double memory_usage();

 protected:
  void compute_temp_target();
  void nhc_temp_integrate();
  void nh_v_temp();

  double dthalf, dt4, dt8;

  char *id_temp;
  class Compute *temperature;
  double t_start, t_stop, t_period;
  double t_current, t_target, ke_target;
  double t_freq, drag, tdrag_factor;
  double boltz, nktv2p, tdof;
  double *eta, *eta_dot;
  double *eta_dotdot;
  double *eta_mass;
  int mtchain;
  double factor_eta;
  int which, got_temp;

  const char *fix_id;
};
}    // namespace LAMMPS_NS

#endif    // LMP_FIX_NVE_MANIFOLD_RATTLE_H
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: There is no manifold named ...

Self-explanatory.  You requested a manifold whose name was not
registered at the factory.

E: Manifold pointer was nullptr for some reason!

This indicates a bug.  The factory was unable to properly create
the requested manifold even though it was registered. Send the
maintainer an e-mail.

E: Manifold ... needs at least ... argument(s)!

Self-explanatory.  Provide enough arguments for the proper
creating of the requested manifold.

E: Parameter pointer was nullptr!

This indicates a bug.  The array that contains the parameters
could not be allocated. Send the maintainer an e-mail.

E: Could not allocate space for arg!

One of the arguments provided was too long (it contained
too many characters)

E: Option ... needs ... argument(s):

Self-explanatory.  Read the documentation of this fix properly.


E: Illegal fix nve/manifold/rattle command! Option ... not recognized!

Self-explanatory.  The given option(s) do not exist.

E: Variable name for fix nve/manifold/rattle does not exist

Self-explanatory.

E: Variable for fix nve/manifold/rattle is invalid style

fix nve/manifold/rattle only supports equal style variables.



*/
