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

/* ----------------------------------------------------------------------
   Contributing authors: Shen,Yuan, Qi,Tingting, and Reed,Evan
   Implementation of the Multi-Scale Shock Method with quantum nuclear effects
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qbmsst,FixQBMSST);
// clang-format on
#else

#ifndef FIX_QBMSST_H
#define FIX_QBMSST_H

#include "fix.h"

namespace LAMMPS_NS {

class FixQBMSST : public Fix {
 public:
  FixQBMSST(class LAMMPS *, int, char **);
  ~FixQBMSST();
  int setmask();
  void init();
  void setup(int);
  void initial_integrate(int);
  void final_integrate();
  double compute_scalar();
  double compute_vector(int);
  void write_restart(FILE *);
  void restart(char *);
  int modify_param(int, char **);
  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

 private:
  // msst parameters
  int direction;      // Direction of shock
  double velocity;    // Velocity of the shock
  double qmass;       // Effective cell mass
  double mu;          // Effective cell viscosity
  double e0;          // Initial energy
  double v0;          // Initial volume
  double p0;          // Initial pressure
  int p0_set;         // Is pressure set
  int v0_set;         // Is volume set
  int e0_set;         // Is energy set
  double tscale;      // Converts thermal energy to compressive strain ke at simulation start
  char *id_temp, *id_press, *id_pe;    // Strings with identifiers of created computes.
  int tflag, pflag, vsflag, peflag;    // Flags to keep track of computes that were created.
  double dtv;                          // update->dt
  double dtf;                          // Full step size
  double dthalf;                       // half step size
  bigint ntotal;                       // atom->natoms
  double boltz, nktv2p, mvv2e;         // Boltzmann factor and unit conversions
  class Compute *temperature;          // Computes created to evaluate
  class Compute *pressure;             // thermodynamic quantities.
  class Compute *pe;
  double **old_velocity;    // Saved velocities.
  int atoms_allocated;      // size of old_velocity
  double dilation[3];
  double omega[3];               // Time derivative of the volume.
  double total_mass;             // Mass of the computational cell
  int kspace_flag;               // 1 if KSpace invoked, 0 if not
  int nrigid;                    // number of rigid fixes
  int *rfix;                     // indices of rigid fixes
  double p_current[3];           // pressure
  double velocity_sum;           // Sum of the velocities squared.
  double lagrangian_position;    // Lagrangian location of computational cell

  // qbmsst parameters
  double t_period, fric_coef;    // friction coefficient
  int seed;                      // seed for the random number generator
  double f_max;                  // frequency cutoff
  int N_f;                       // number of frequency grid
  double eta;                    // coupling coefficient between shock and the qtb
  int beta;                      // average beta steps before updating the qtb temperature
  double t_init;                 // initial qtb temperature
  int qtb_set;                   // 1 if its a restarting qbmsst, 0 if not
  int counter_l, counter_mu;     // counter l and mu
  double t_current;              // qtb temperature
  double h_timestep;             // time step to update the random forces
  int alpha;                     // number of time steps to update the random forces
  class RanMars *random;         // random number generator
  double *gfactor;               // factors of random forces
  double *omega_H, *time_H;      // H gives the desired power spectrum
  double **random_array_0, **random_array_1,
      **random_array_2;    // random number arrays give independence between atoms and directions
  double **fran;           // random forces
  double old_eavg;         // average previous energies

  // functions
  void couple();
  void remap(int);
  void check_alloc(int n);
  double compute_etotal();
  double compute_egrand();
  double compute_vol();
  double compute_vsum();
  double compute_hugoniot();
  double compute_rayleigh();
  double compute_lagrangian_speed();
  double compute_lagrangian_position();
};

}    // namespace LAMMPS_NS

#endif
#endif
