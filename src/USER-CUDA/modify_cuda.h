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

#ifndef LMP_MODIFY_CUDA_H
#define LMP_MODIFY_CUDA_H

#include <cstdio>
#include "modify.h"

namespace LAMMPS_NS {

class ModifyCuda : public Modify {
 public:

  int n_initial_integrate_cuda;
  int n_post_integrate_cuda;
  int n_pre_exchange_cuda;
  int n_pre_neighbor_cuda;
  int n_pre_force_cuda;
  int n_post_force_cuda;
  int n_final_integrate_cuda;
  int n_end_of_step_cuda;
  int n_thermo_energy_cuda;

  int n_initial_integrate_host;
  int n_post_integrate_host;
  int n_pre_exchange_host;
  int n_pre_neighbor_host;
  int n_pre_force_host;
  int n_post_force_host;
  int n_final_integrate_host;
  int n_end_of_step_host;
  int n_thermo_energy_host;

  ModifyCuda(class LAMMPS *);
  ~ModifyCuda();
  void init();
  void initial_integrate(int);
  void post_integrate();
  //void pre_decide();
  void pre_exchange();
  void pre_neighbor();
  void setup_pre_force(int);
  void pre_force(int);
  void post_force(int);
  void final_integrate();
  void end_of_step();
  double thermo_energy();


 protected:
  class Cuda *cuda;

  // lists of fixes to apply at different stages of timestep

  // list of cuda fixes
  int *list_initial_integrate_cuda;
  int *list_post_integrate_cuda;
  int *list_pre_exchange_cuda;
  int *list_pre_neighbor_cuda;
  int *list_pre_force_cuda;
  int *list_post_force_cuda;
  int *list_final_integrate_cuda;
  int *list_end_of_step_cuda;
  int *list_thermo_energy_cuda;
  int *end_of_step_every_cuda;

  void list_init_end_of_step_cuda(int, int &, int *&);
};

}

#endif
