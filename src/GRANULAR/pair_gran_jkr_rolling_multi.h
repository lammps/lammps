/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(gran/jkr/rolling/multi,PairGranJKRRollingMulti)

#else

#ifndef LMP_PAIR_GRAN_JKR_ROLLING_MULTI_H
#define LMP_PAIR_GRAN_JKR_ROLLING_MULTI_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGranJKRRollingMulti : public Pair {
public:
  PairGranJKRRollingMulti(class LAMMPS *);
  virtual ~PairGranJKRRollingMulti();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **); // Made Virtual by IS Oct 7 2017
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
  double **E,**G,**alpha,**gamman,**muS,**Ecoh,**kR,**muR,**etaR,**cut;
  int **normaldamp, **rollingdamp;
  double dt;
  int freeze_group_bit;
  int history;

  int neighprev;
  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;

  class FixNeighHistory *fix_history;

  // storage of rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  virtual void allocate(); // Made Virtual by IS Oct 7 2017

private:
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
