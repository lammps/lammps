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

#ifndef FIX_TTM_H
#define FIX_TTM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTTM : public Fix {
 public:
  FixTTM(class LAMMPS *, int, char **);
  ~FixTTM();
  int setmask();
  void init();
  void read_initial_electron_temperatures();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void reset_dt();
  void update_electron_temperatures();

 private:
  int me;
  int nfileevery;
  int nlevels_respa;
  class RanMars *random;
  FILE *fp, *fpr;
  int nxnodes,nynodes,nznodes,total_nnodes;
  int ***nsum, ***nsum_prime, ***nsum_all, ***nsum_prime_all, ***T_initial_set;
  double *gfactor1,*gfactor2,*ratio;
  double ***T_electron, ***T_a, ***T_a_prime, ***g_p, ***g_s;
  double ***sum_vsq, ***sum_vsq_prime, ***sum_mass_vsq, ***sum_mass_vsq_prime;
  double ***sum_vsq_all, ***sum_vsq_prime_all, ***sum_mass_vsq_all, ***sum_mass_vsq_prime_all;
  double ***T_electron_old;
  double electronic_specific_heat,electronic_density,electronic_thermal_conductivity;
  double gamma_p,gamma_s,v_0,v_0_sq;
};

}

#endif
