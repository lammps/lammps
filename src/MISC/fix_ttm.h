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

#ifdef FIX_CLASS

FixStyle(ttm,FixTTM)

#else

#ifndef LMP_FIX_TTM_H
#define LMP_FIX_TTM_H

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
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Cannot open fix ttm file %s

The output file for the fix ttm command cannot be opened.  Check that
the path and name are correct.

E: Invalid random number seed in fix ttm command

Random number seed must be > 0.

E: Fix ttm electronic_specific_heat must be > 0.0

Self-explanatory.

E: Fix ttm electronic_density must be > 0.0

Self-explanatory.

E: Fix ttm electronic_thermal_conductivity must be >= 0.0

Self-explanatory.

E: Fix ttm gamma_p must be > 0.0

Self-explanatory.

E: Fix ttm gamma_s must be >= 0.0

Self-explanatory.

E: Fix ttm v_0 must be >= 0.0

Self-explanatory.

E: Fix ttm number of nodes must be > 0

Self-explanatory.

E: Cannot use fix ttm with 2d simulation

This is a current restriction of this fix due to the grid it creates.

E: Cannot use nonperiodic boundares with fix ttm

This fix requires a fully periodic simulation box.

E: Cannot use fix ttm with triclinic box

This is a current restriction of this fix due to the grid it creates.

E: Electronic temperature dropped below zero

Something has gone wrong with the fix ttm electron temperature model.

E: Fix ttm electron temperatures must be > 0.0

Self-explanatory.

E: Initial temperatures not all set in fix ttm

Self-explantory.

W: Too many inner timesteps in fix ttm

Self-explanatory.

*/
