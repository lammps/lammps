// clang-format off
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

#ifndef LMP_MODIFY_KOKKOS_H
#define LMP_MODIFY_KOKKOS_H

#include "modify.h"

namespace LAMMPS_NS {

class ModifyKokkos : public Modify {
 public:
  ModifyKokkos(class LAMMPS *);

  void setup(int) override;
  void setup_pre_exchange() override;
  void setup_pre_neighbor() override;
  void setup_post_neighbor() override;
  void setup_pre_force(int) override;
  void setup_pre_reverse(int, int) override;
  void initial_integrate(int) override;
  void post_integrate() override;
  void pre_decide();
  void pre_exchange() override;
  void pre_neighbor() override;
  void post_neighbor() override;
  void pre_force(int) override;
  void pre_reverse(int,int) override;
  void post_force(int) override;
  void final_integrate() override;
  void end_of_step() override;
  double energy_couple() override;
  double energy_global() override;
  void energy_atom(int, double *) override;
  void post_run() override;

  void setup_pre_force_respa(int, int) override;
  void initial_integrate_respa(int, int, int) override;
  void post_integrate_respa(int, int) override;
  void pre_force_respa(int, int, int) override;
  void post_force_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;

  void min_pre_exchange() override;
  void min_pre_neighbor() override;
  void min_post_neighbor() override;
  void min_pre_force(int) override;
  void min_pre_reverse(int,int) override;
  void min_post_force(int) override;

  double min_energy(double *) override;
  void min_store() override;
  void min_step(double, double *) override;
  void min_clearstore() override;
  void min_pushstore() override;
  void min_popstore() override;
  double max_alpha(double *) override;
  int min_dof() override;
  int min_reset_ref() override;

 protected:

};

}

#endif

