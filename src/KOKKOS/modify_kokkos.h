/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

  void init() override;
  void setup(int) override;
  void setup_pre_exchange() override;
  void setup_pre_neighbor() override;
  void setup_post_neighbor() override;
  void setup_pre_force(int) override;
  void setup_pre_reverse(int, int) override;
  void initial_integrate(int) override;
  void post_integrate() override;
  void pre_exchange() override;
  void pre_neighbor() override;
  void post_neighbor() override;
  void pre_force(int) override;
  void pre_reverse(int, int) override;
  void post_force(int) override;
  void final_integrate() override;
  void fused_integrate(int) override;
  void end_of_step() override;
  double energy_couple() override;
  double energy_global() override;
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
  void min_pre_reverse(int, int) override;
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

  int check_fuse_integrate();

 protected:
  int fuse_integrate_ok;    // cached result of check_fuse_integrate(), set in init()

  // invoke func(obj) on a fix or compute with the per-atom data synchronized
  // to its execution space before and marked as modified afterwards

  template <typename T, typename Func> void sync_call(T *, Func &&);

  // same as Modify::for_each_fix() but with the data synchronization

  template <typename Func> void for_each_fix_kk(int, Func &&);
};

}    // namespace LAMMPS_NS

#endif
