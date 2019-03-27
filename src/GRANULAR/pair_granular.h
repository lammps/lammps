/* ----------------------------------------------------------
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

PairStyle(granular,PairGranular)

#else

#ifndef LMP_PAIR_GRANULAR_H
#define LMP_PAIR_GRANULAR_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGranular : public Pair {
 public:
  PairGranular(class LAMMPS *);
  ~PairGranular();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void reset_dt();
  double single(int, int, int, int, double, double, double, double &);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 protected:
  double dt;
  int freeze_group_bit;
  int use_history;

  int neighprev;
  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;
  double **cut;

  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  void allocate();
  void transfer_history(double*, double*);

 private:
  int size_history;
  int *history_transfer_factors;

  // model choices
  int **normal_model, **damping_model;
  int **tangential_model, **roll_model, **twist_model;

  // history flags
  int normal_history, tangential_history, roll_history, twist_history;

  // indices of history entries
  int normal_history_index;
  int tangential_history_index;
  int roll_history_index;
  int twist_history_index;

  // per-type material coefficients
  double **Emod, **poiss, **Gmod;

  // per-type coefficients, set in pair coeff command
  double ***normal_coeffs;
  double ***tangential_coeffs;
  double ***roll_coeffs;
  double ***twist_coeffs;

  // optional user-specified global cutoff, per-type user-specified cutoffs
  double **cutoff_type;
  double cutoff_global;

  double mix_stiffnessE(double, double, double, double);
  double mix_stiffnessG(double, double, double, double);
  double mix_geom(double, double);
  double pulloff_distance(double, double, int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

 */
