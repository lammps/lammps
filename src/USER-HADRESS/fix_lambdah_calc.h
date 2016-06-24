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

#ifdef FIX_CLASS

FixStyle(lambdah/calc,FixLambdaHCalc)

#else

#ifndef LMP_FIX_LAMBDAH_CALC_H
#define LMP_FIX_LAMBDAH_CALC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLambdaHCalc : public Fix {
 public:
  FixLambdaHCalc(class LAMMPS *, int, char **);
  ~FixLambdaHCalc();
  int setmask();
  void init();
  void setup(int);

  void post_integrate();
  double memory_usage();


 double Length_Hyb,Length_ATRegion,Length_CGRegion;
 int Hybrid_Style;
 int Pressure_Comp_Flag,Pressure_Bin_Num,Pressure_Update_Frequency,Pressure_Update_Time_End, Pressure_Update_Time_Begin;
 double Pressure_lambda_Increment;
 int Density_Comp_Flag,Density_Compensation_Run,Density_Bin_Num,Density_Update_Frequency,Density_Update_Time_End,Density_Update_Time_Begin;
 double Density_Bin_Size;
 int Comp_Counter_H,Density_Counter_H;
 double **Mean_grad_Comp_Density_Conv_H;

 double Density_Ref;
 double center_box[3];
 int me;
 double x0lo,x0hi,x1lo,x1hi,x2lo,x2hi,x0BoxSize;
 double Density_Sigma_Gauss;
 double Comp_Density_Scaling_factor_H;
 int Density_Gauss_Int_Range;

 double xmin_AT,xmax_AT;
 private:
 int nmolecules;
 tagint idlo,idhi;
 double *massproc_H,*masstotal_H;
 double **com_H,**comall_H;
 double Density_Fluctuation;
 int **Comp_Density_Num_H, **Comp_Density_Num_all_H;
 double **Int_Mean_Density_H, **Mean_Density_H, **grad_Comp_Density_Conv_H, **Mean_Comp_Density_Conv_H;
 int Load_File_Flag;
 double R_Start_Hybrid_1,R_Start_Hybrid_2,S_Start_Hybrid,r,cosr,sinr,sinx,cosx;
 double *lambdaCM, **gradlambdaCM;
 int *molmap_H;                 // convert molecule ID to local index

 int molecules_in_group(tagint &, tagint &);

 void Print_Compensation_Density();
 void Load_Compensation_Density();
 void Clear_File_Compensation_Density();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Variable name for fix indent does not exist

Self-explanatory.

E: Variable for fix indent is invalid style

Only equal-style variables can be used.

E: Variable for fix indent is not equal style

Only equal-style variables can be used.

*/
