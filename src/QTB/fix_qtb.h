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
   Contributing authors: Shen,Yuan, Qi,Tingting, and Reed,Evan
   Implementation of the colored thermostat for quantum nuclear effects
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(qtb,FixQTB);
// clang-format on
#else

#ifndef LMP_FIX_QTB_H
#define LMP_FIX_QTB_H

#include "fix.h"

namespace LAMMPS_NS {

class FixQTB : public Fix {
 public:
  FixQTB(class LAMMPS *, int, char **);
  ~FixQTB() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  int modify_param(int, char **) override;
  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 private:
  // qtb parameters
  int counter_mu;                // counter l and mu
  double t_period, fric_coef;    // friction coefficient
  int seed;                      // seed for the random number generator
  double f_max;                  // frequency cutoff
  int N_f;                       // number of frequency grid
  double t_target;               // target qtb temperature
  char *id_temp;
  class Compute *temperature;
  double h_timestep;              // time step to update the random forces
  int alpha;                      // number of time steps to update the random forces
  class RanMars *random;          // random number generator
  double *gfactor1, *gfactor3;    // factors of frictions and random forces
  double *omega_H, *time_H;       // H gives the desired power spectrum
  // random number arrays give independence between atoms and directions
  double **random_array_0, **random_array_1, **random_array_2;
  int nlevels_respa;
  double **fran, fsum[3], fsumall[3];    // random forces and their sums
};

}    // namespace LAMMPS_NS

#endif
#endif
