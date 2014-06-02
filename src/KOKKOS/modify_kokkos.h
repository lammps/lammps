/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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
  ~ModifyKokkos() {}
  void setup(int);
  void setup_pre_exchange();
  void setup_pre_neighbor();
  void setup_pre_force(int);
  void initial_integrate(int);
  void post_integrate();
  void pre_decide();
  void pre_exchange();
  void pre_neighbor();
  void pre_force(int);
  void post_force(int);
  void final_integrate();
  void end_of_step();
  double thermo_energy();
  void post_run();

  void setup_pre_force_respa(int, int);
  void initial_integrate_respa(int, int, int);
  void post_integrate_respa(int, int);
  void pre_force_respa(int, int, int);
  void post_force_respa(int, int, int);
  void final_integrate_respa(int, int);

  void min_pre_exchange();
  void min_pre_neighbor();
  void min_pre_force(int);
  void min_post_force(int);

  double min_energy(double *);
  void min_store();
  void min_step(double, double *);
  void min_clearstore();
  void min_pushstore();
  void min_popstore();
  double max_alpha(double *);
  int min_dof();
  int min_reset_ref();

 protected:
  class AtomKokkos *atomKK;
};

}

#endif

/* ERROR/WARNING messages:

*/
