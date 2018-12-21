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
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
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

  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  virtual void allocate();
  int beyond_contact;


  // comment next line to turn off templating
/*#define TEMPLATED_PAIR_GRANULAR
#ifdef TEMPLATED_PAIR_GRANULAR
  template < int Tp_coeff_types,
  int Tp_normal, int Tp_damping, int Tp_tangential,
  int Tp_rolling, int Tp_twisting >
  void compute_templated(int eflag, int vflag);
#else
*/
  void compute_untemplated(
      int,
      int, int, int,
      int, int,
      int, int);
//#endif

 private:
  int coeff_types;
  int size_history;

  //Per-type models
  int **normal, **damping, **tangential, **rolling, **twisting;

  int normal_global, damping_global;
  int tangential_global, rolling_global, twisting_global;

  int tangential_history, rolling_history, twisting_history;
  int tangential_history_index;
  int rolling_history_index;
  int twisting_history_index;

  double *normal_coeffs_global;
  double *tangential_coeffs_global;
  double *rolling_coeffs_global;
  double *twisting_coeffs_global;

  double ***normal_coeffs;
  double ***tangential_coeffs;
  double ***rolling_coeffs;
  double ***twisting_coeffs;

  double mix_stiffnessE(double Eii, double Ejj, double Gii, double Gjj);
  double mix_stiffnessG(double Eii, double Ejj, double Gii, double Gjj);
  double mix_geom(double valii, double valjj);
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
