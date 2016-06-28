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

PairStyle(lj/cut/hars/at,PairLJCutHARSAT)

#else

#ifndef LMP_PAIR_LJ_CUT_HARS_AT_H
#define LMP_PAIR_LJ_CUT_HARS_AT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutHARSAT : public Pair {
 public:
  PairLJCutHARSAT(class LAMMPS *);
  virtual ~PairLJCutHARSAT();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);

  void AT_Print_Compensation_Energy();
  void AT_Update_Compensation_Energy();

 protected:
  double cut_global;
  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double *cut_respa;
  class FixLambdaHCalc *lambda_H_fix;

  int AllAtomistic;

  virtual void allocate();

//private:
  int me;

  int H_AdResS_allocated;

  int **Comp_Energy_Num_H,**Comp_Energy_Num_all_H, Comp_Counter_H;

  double AT_lambda_Increment;
  int AT_Bin_Num, AT_Update_Frequency, AT_Update_Time_End, AT_Update_Time_Begin;
  int AT_Pressure_Compensation_Run;
  int AT_Pressure_Comp_Flag;

  double *AT_center_box,AT_x0lo;
  int AT_Hybrid_Style;
  double **Int_Mean_Energy_H, **Comp_Energy_H, **Comp_Energy_all_H, **Mean_Energy_H, **Mean_Comp_Energy_H;

   int nmolecules;
   tagint idlo,idhi;

   double *AT_massproc_H,*AT_masstotal_H;

   double **AT_mol_f_H, **AT_mol_f_all_H;

   int *AT_molmap_H;                 // convert molecule ID to local index
   int Load_File_Flag;

   int molecules_in_group(tagint &, tagint &);
   void Load_Compensation_Pressure();
   void H_AdResS_Allocation();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
