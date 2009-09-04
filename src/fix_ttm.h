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
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void post_force_setup(int);
  void post_force_respa_setup(int, int, int);
  void end_of_step();
  void reset_dt();
  void write_restart(FILE *);
  void restart(char *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  double memory_usage();
  void grow_arrays(int); 
  double compute_vector(int);

 private:
  int me;
  int nfileevery;
  int nlevels_respa;
  int seed;
  class RanMars *random;
  FILE *fp,*fpr;
  int nxnodes,nynodes,nznodes,total_nnodes;
  int ***nsum;
  int ***nsum_all,***T_initial_set;
  double *gfactor1,*gfactor2,*ratio;
  double **flangevin; 
  double ***T_electron,***T_electron_old;
  double ***sum_vsq,***sum_mass_vsq;
  double ***sum_vsq_all,***sum_mass_vsq_all;
  double ***net_energy_transfer,***net_energy_transfer_all;
  double electronic_specific_heat,electronic_density;
  double electronic_thermal_conductivity;
  double gamma_p,gamma_s,v_0,v_0_sq;

  void read_initial_electron_temperatures();
};

}

#endif
