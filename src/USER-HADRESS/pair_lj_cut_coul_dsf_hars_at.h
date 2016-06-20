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

PairStyle(lj/cut/coul/dsf/hars/at,PairLJCutCoulDSFHARSAT)

#else

#ifndef LMP_PAIR_LJ_CUT_COUL_DSF_HARS_AT_H
#define LMP_PAIR_LJ_CUT_COUL_DSF_HARS_AT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulDSFHARSAT : public Pair {
 public:
  PairLJCutCoulDSFHARSAT(class LAMMPS *);
  ~PairLJCutCoulDSFHARSAT();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

  void AT_Print_Compensation_Energy();
  void AT_Update_Compensation_Energy();

 protected:
  double cut_lj_global;
  double **cut_lj,**cut_ljsq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  
  double cut_coul,cut_coulsq;
  double alpha;
  double f_shift,e_shift;
  
  void allocate();
  class FixLambdaHCalc *lambda_H_fix;

  int AllAtomistic;

  //private:
    int me;

    int H_AdResS_allocated;

    int **Comp_Energy_Num_H,**Comp_Energy_Num_all_H, Comp_Counter_H;

    double AT_lambda_Increment;
    int AT_Bin_Num, AT_Update_Frequency, AT_Update_Time_End, AT_Update_Time_Begin;
    int AT_Pressure_Compensation_Run;
    int AT_Pressure_Comp_Flag;

    int AT_Restart_Time_Step;
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

E: Pair style lj/cut/coul/dsf requires atom attribute q

The atom style defined does not have these attributes.

*/
