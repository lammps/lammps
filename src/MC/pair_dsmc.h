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

#ifdef PAIR_CLASS

PairStyle(dsmc,PairDSMC)

#else

#ifndef LMP_PAIR_DSMC_H
#define LMP_PAIR_DSMC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairDSMC : public Pair {
 public:
  PairDSMC(class LAMMPS *);
  virtual ~PairDSMC();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

 private:
  double cut_global;
  double **cut;
  double **sigma;

  double cellx;
  double celly;
  double cellz;
  int ncellsx;
  int ncellsy;
  int ncellsz;
  int total_ncells;
  int total_number_of_collisions;
  int recompute_vsigmamax_stride;
  int vsigmamax_samples;
  double T_ref;
  double kT_ref;
  double two_pi;
  double max_cell_size;

  int seed;
  int number_of_A;
  int number_of_B;
  int max_particle_list;

  class RanMars *random;

  int **particle_list;
  int **first;
  int **number;

  double **V_sigma_max;

  int max_particles;
  int *next_particle;

  int itype;
  int jtype;

  double imass;
  double jmass;
  double total_mass;
  double reduced_mass;
  double imass_tmass;
  double jmass_tmass;
  double vol;
  double weighting;

  void allocate();
  void recompute_V_sigma_max(int);
  double V_sigma(int, int);
  void scatter_random(int, int, int);
  int convert_double_to_equivalent_int(double);

  inline void subtract3d(const double *v1, const double *v2, double *v3) {
    v3[0] = v2[0] - v1[0];
    v3[1] = v2[1] - v1[1];
    v3[2] = v2[2] - v1[2];
  }
  
  inline double dot3d(const double *v1, const double *v2) {
    return( v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] ); 
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: Pair dsmc: num_of_collisions > number_of_A

Collision model in DSMC is breaking down.

W: Pair dsmc: num_of_collisions > number_of_B

Collision model in DSMC is breaking down.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Tried to convert a double to int, but input_double > INT_MAX

Self-explanatory.

*/
