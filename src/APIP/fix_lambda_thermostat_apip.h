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
/* ----------------------------------------------------------------------
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(lambda_thermostat/apip,FixLambdaThermostatAPIP);
// clang-format on
#else

#ifndef LMP_FIX_LAMBDA_THERMOSTAT_APIP_H
#define LMP_FIX_LAMBDA_THERMOSTAT_APIP_H

#include "fix.h"
#include <random>

namespace LAMMPS_NS {

class FixLambdaThermostatAPIP : public Fix {
 public:
  FixLambdaThermostatAPIP(class LAMMPS *, int, char **);
  ~FixLambdaThermostatAPIP() override;

  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void post_force(int) override;
  void end_of_step() override;
  double memory_usage() override;
  void reset_dt() override;
  double compute_vector(int) override;

 protected:
  void apply_thermostat();
  void calculate_energy_change();
  double calculate_kinetic_energy(int);
  void local_neighbour_list();
  void init_peratom_stats();

  double dtf;    // constant for time integration

  class NeighList *list;
  double *energy_change_atom;    // energy violation compared to const lambda case
  double **peratom_stats;        // peratom output
  int nmax_energy;               // number of atoms for which energy_change_atom is allocated
  int nmax_stats;                // number of atoms for which peratom_stats is allocated

  // own neighbour list
  int pgsize;                // size of neighbor page
  int oneatom;               // max # of neighbors for one atom
  int nmax_list;             // size of numneigh, firstneigh arrays
  int *local_numneigh;       // # of pair neighbors for each atom
  int **local_firstneigh;    // ptr to 1st neighbor of each atom
  MyPage<int> *ipage;        // neighbor list pages
  int *jlist_copy;           // jlist for one atom

  int reduceflag;              // 1/0 calculation of compute_vector required/not required
  double outvec[6];            // vector returned by compute_vector
  double energy_change_pot;    // energy conservation violation of all atoms
  double energy_change_kin;    // energy conservation violation of all atoms
  int n_energy_differences;    // number of atoms whose energy has changed compared to the constant lambda case
  double sum_energy_change;       // absolute value of all energy changes due to rescaling
  double sum_energy_violation;    // energy that could not be compensated, accumulated over time
  int n_energy_violation;    // number of atoms whose energy could not be compensated, accumulated over time

  int rescaling_N_neighbours;    // requested neighbour list size used for rescaling

  std::mt19937 random_mt;    // mersenne twister for shuffle

  bool update_stats;    // true(false) peratom output needs(not) to be calculated
};

}    // namespace LAMMPS_NS

#endif
#endif
