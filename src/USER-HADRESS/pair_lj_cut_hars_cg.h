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

PairStyle(lj/cut/hars/cg,PairLJCutHARSCG)

#else

#ifndef LMP_PAIR_LJ_CUT_HARS_CG_H
#define LMP_PAIR_LJ_CUT_HARS_CG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutHARSCG : public Pair {
 public:
  PairLJCutHARSCG(class LAMMPS *);
  virtual ~PairLJCutHARSCG();
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

  void CG_Print_Compensation_Energy();
  void CG_Update_Compensation_Energy();

 protected:
  double cut_global;
  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double *cut_respa;
  class FixLambdaHCalc *lambda_H_fix;
  int AllCoarseGrained;

  double *CG_Ave_Mean_Density_H,CG_grad_Ave_Mean_Density_H;
 //private:

  int me;
  virtual void allocate();

  int H_AdResS_allocated;
  int **Comp_Energy_Num_H,**Comp_Energy_Num_all_H, Comp_Counter_H;

  double CG_lambda_Increment;
  int CG_Bin_Num, CG_Update_Frequency, CG_Update_Time_End, CG_Update_Time_Begin;
  int CG_Pressure_Compensation_Run, Density_Compensation_Run;

  int CG_Pressure_Comp_Flag, CG_Density_Comp_Flag;

  double CG_Density_Bin_Size,**CG_Mean_grad_Comp_Density_Conv_H,CG_x0BoxSize;
  int CG_Density_Bin_Num,CG_Density_Update_Frequency,CG_Density_Update_Time_Begin,CG_Density_Update_Time_End;
  double **Int_Mean_Energy_H, **Comp_Energy_H, **Comp_Energy_all_H, **Mean_Energy_H, **Mean_Comp_Energy_H;

  double *CG_center_box, CG_x0lo;
  int CG_Hybrid_Style;
    int nmolecules;
    tagint idlo,idhi;

    double *massproc_H,*masstotal_H;

    double **mol_f_H, **mol_f_all_H;

   int CG_Restart_Time_Step;
    //double **drift_f_H, **drift_f_all_H;
    //int nbin_H;
    int *molmap_H;                 // convert molecule ID to local index
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
