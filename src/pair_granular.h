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
  virtual ~PairGranular();

  void compute(int, int);
  // comment next line to turn off templating
#define TEMPLATED_PAIR_GRANULAR
#ifdef TEMPLATED_PAIR_GRANULAR
  template < int Tp_normal, int Tp_damping, int Tp_tangential,
             int Tp_roll, int Tp_twist>
  void compute_templated(int, int);
#else
  void compute_untemplated(int, int, int, int, int,
      int, int);
#endif

  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void reset_dt();
  virtual double single(int, int, int, int, double, double, double, double &);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  double memory_usage();

 protected:
  double cut_global;
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

  virtual void allocate();

private:
  int size_history;

  //Models
  int normal, damping, tangential, roll, twist;

  //History flags
  int tangential_history, roll_history, twist_history;

  //Indices of history entries
  int tangential_history_index, roll_history_index, twist_history_index;

  //Flags for whether model choices have been set
  int normal_set, tangential_set, damping_set, roll_set, twist_set;

  //Per-type coefficients, set in pair coeff command
  double ***normal_coeffs;
  double ***tangential_coeffs;
  double ***roll_coeffs;
  double ***twist_coeffs;

  //Optional user-specified global cutoff
  double cutoff_global;

  double mix_stiffnessE(double Eii, double Ejj, double Gii, double Gjj);
  double mix_stiffnessG(double Eii, double Ejj, double Gii, double Gjj);
  double mix_geom(double valii, double valjj);
  double pulloff_distance(double radius, int itype);
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
