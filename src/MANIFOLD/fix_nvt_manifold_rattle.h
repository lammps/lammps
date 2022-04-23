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
  ~FixNVTManifoldRattle() override;

  void initial_integrate(int) override;
  void final_integrate() override;
  void init() override;
  void reset_dt() override;
  int setmask() override;
  void setup(int) override;    // Not needed for fixNVE but is for fixNVT
  double memory_usage() override;

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
