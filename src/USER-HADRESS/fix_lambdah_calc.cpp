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

/* ----------------------------------------------------------------------
   Contributing author: Maziar Heidari (Max Planck Institute for Polymer Research)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_lambdah_calc.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "output.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "comm.h"
#include "citeme.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define BIG MAXTAGINT

static const char cite_HAdResS[] =
  "@Article{Heidari et al.2016\n"
  " author = {M. Heidari, R. Cortes-Huerto, D. Donadio and R. Potestio},\n"
  " title = {Accurate and general treatment of electrostatic interaction in Hamiltonian adaptive resolution simulations},\n"
  " journal = {Eur. Phys. J. Special Topics},\n"
  " year =    2016,\n"
  " volume =  Submitted,\n"
  " pages =   {}\n"
  "}\n\n";


/* ---------------------------------------------------------------------- */

FixLambdaHCalc::FixLambdaHCalc(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{


    if (lmp->citeme) lmp->citeme->add(cite_HAdResS);

    massproc_H = NULL;
    masstotal_H = NULL;
    com_H = NULL;
    comall_H = NULL;
    molmap_H = NULL;
    lambdaCM = NULL;
    gradlambdaCM = NULL;
    Comp_Density_Num_H = NULL;
    Comp_Density_Num_all_H = NULL;
    Mean_Density_H = NULL;
    Int_Mean_Density_H = NULL;
    Mean_Comp_Density_Conv_H = NULL;
    grad_Comp_Density_Conv_H = NULL;
    Mean_grad_Comp_Density_Conv_H = NULL;



  me = comm->me;

  if (narg < 8) error->all(FLERR,"Illegal fix lambdah/calc command");

  atom->nmoltypesH  = force->numeric(FLERR,arg[3]);
  Length_Hyb  = force->numeric(FLERR,arg[4]);
  Length_ATRegion  = force->numeric(FLERR,arg[5]);
  Pressure_Comp_Flag  = force->numeric(FLERR,arg[6]);
  Pressure_lambda_Increment  = force->numeric(FLERR,arg[7]);
  Pressure_Update_Frequency  = force->numeric(FLERR,arg[8]);
  Pressure_Update_Time_Begin  = force->numeric(FLERR,arg[9]);
  Pressure_Update_Time_End  = force->numeric(FLERR,arg[10]);

  if (strcmp(arg[11],"slab") == 0) Hybrid_Style = 0;
  else if (strcmp(arg[11],"sphere") == 0) Hybrid_Style = 1;
  else if (strcmp(arg[11],"cylinder") == 0) Hybrid_Style = 2;
  else error->all(FLERR,"Illegal fix lambdah/calc command");

  Density_Comp_Flag  = force->numeric(FLERR,arg[12]);
  Density_Bin_Size  = force->numeric(FLERR,arg[13]);
  Density_Update_Frequency  = force->numeric(FLERR,arg[14]);
  Density_Update_Time_Begin  = force->numeric(FLERR,arg[15]);
  Density_Update_Time_End  = force->numeric(FLERR,arg[16]);
  Density_Sigma_Gauss = force->numeric(FLERR,arg[17]);
  Density_Gauss_Int_Range = force->numeric(FLERR,arg[18]);
  Density_Ref = force->numeric(FLERR,arg[19]);
  Comp_Density_Scaling_factor_H =  force->numeric(FLERR,arg[20]);
  Load_File_Flag =  force->numeric(FLERR,arg[21]);

  x0lo = domain->boxlo[0];
  x0hi = domain->boxhi[0];
  x1lo = domain->boxlo[1];
  x1hi = domain->boxhi[1];
  x2lo = domain->boxlo[2];
  x2hi = domain->boxhi[2];

  center_box[0] = (x0hi + x0lo)/2.0;
  center_box[1] = (x1hi + x1lo)/2.0;
  center_box[2] = (x2hi + x2lo)/2.0;


  x0BoxSize = x0hi - x0lo;
  Length_CGRegion = x0BoxSize - 2*Length_Hyb - Length_ATRegion;
  R_Start_Hybrid_1 = x0lo + Length_CGRegion/2.0;
  R_Start_Hybrid_2 = x0lo + x0BoxSize - (Length_CGRegion/2.0 + Length_Hyb);
  S_Start_Hybrid = Length_ATRegion;
  Pressure_Bin_Num = 1.0 / Pressure_lambda_Increment;
  xmin_AT = R_Start_Hybrid_1 + Length_Hyb;
  xmax_AT= R_Start_Hybrid_1 + Length_Hyb + Length_ATRegion;
  if(Hybrid_Style==0)Density_Bin_Num = floor(x0BoxSize / Density_Bin_Size);
  else if(Hybrid_Style==1)Density_Bin_Num = floor(sqrt(pow(0.5*(x0hi-x0lo),2.0)+pow(0.5*(x1hi-x1lo),2.0)+pow(0.5*(x2hi-x2lo),2.0)) / Density_Bin_Size);
  else if(Hybrid_Style==2)Density_Bin_Num = floor(sqrt(pow(0.5*(x0hi-x0lo),2.0)+pow(0.5*(x1hi-x1lo),2.0)) / Density_Bin_Size);


  Comp_Counter_H = 0;
  Density_Counter_H = 0;
  Density_Compensation_Run = 0;

  memory->create(Comp_Density_Num_H,Density_Bin_Num,atom->nmoltypesH+1,"lambdaH/calc:Comp_Density_Num_H");
  memory->create(Comp_Density_Num_all_H,Density_Bin_Num,atom->nmoltypesH+1,"lambdaH/calc:Comp_Density_Num_all_H");

  if(me==0){
    if (screen){
      fprintf(screen,"nmoltypes= %d\n",atom->nmoltypesH);
      fprintf(screen,"Length_Hyb= %f\n",Length_Hyb);
      fprintf(screen,"Length_ATRegion= %f\n",Length_ATRegion);
      fprintf(screen,"Pressure_Comp_Flag= %d\n",Pressure_Comp_Flag);
      fprintf(screen,"Pressure_lambda_Increment= %f\n",Pressure_lambda_Increment);
      fprintf(screen,"Pressure_Update_Frequency= %d\n",Pressure_Update_Frequency);
      fprintf(screen,"Pressure_Update_Time_Begin= %d\n",Pressure_Update_Time_Begin);
      fprintf(screen,"Pressure_Update_Time_End= %d\n",Pressure_Update_Time_End);
      fprintf(screen,"Density_Comp_Flag= %d\n",Density_Comp_Flag);
      fprintf(screen,"Density_Bin_Size= %f\n",Density_Bin_Size);
      fprintf(screen,"Density_Update_Frequency= %d\n",Density_Update_Frequency);
      fprintf(screen,"Density_Update_Time_Begin= %d\n",Density_Update_Time_Begin);
      fprintf(screen,"Density_Update_Time_End= %d\n",Density_Update_Time_End);
      fprintf(screen,"Density_Sigma_Gauss= %f\n",Density_Sigma_Gauss);
      fprintf(screen,"Density_Gauss_Int_Range= %d\n",Density_Gauss_Int_Range);
      fprintf(screen,"Density_Ref= %f\n",Density_Ref);
      fprintf(screen,"Comp_Density_Scaling_factor_H = %f\n",Comp_Density_Scaling_factor_H);
      fprintf(screen,"Load_File_Flag = %d\n",Load_File_Flag);
      fprintf(screen,"Center_box = %f %f %f\n",center_box[0],center_box[1],center_box[2]);
      fprintf(screen,"Density_Bin_Size = %f\n",Density_Bin_Size);
      fprintf(screen,"Density_Bin_Num = %d\n",Density_Bin_Num);
      fprintf(screen,"x0lo= %f\n",x0lo);
      fprintf(screen,"x0hi= %f\n",x0hi);
      fprintf(screen,"x0BoxSize= %f\n",x0BoxSize);
      fprintf(screen,"d1= %f\n",R_Start_Hybrid_1);
      fprintf(screen,"d2= %f\n",R_Start_Hybrid_2);
      fprintf(screen,"moltype%d\n",atom->nmoltypesH);

    }
    if (logfile){
        fprintf(logfile,"nmoltypes= %d\n",atom->nmoltypesH);
        fprintf(logfile,"Length_Hyb= %f\n",Length_Hyb);
        fprintf(logfile,"Length_ATRegion= %f\n",Length_ATRegion);
        fprintf(logfile,"Pressure_Comp_Flag= %d\n",Pressure_Comp_Flag);
        fprintf(logfile,"Pressure_lambda_Increment= %f\n",Pressure_lambda_Increment);
        fprintf(logfile,"Pressure_Update_Frequency= %d\n",Pressure_Update_Frequency);
        fprintf(logfile,"Pressure_Update_Time_Begin= %d\n",Pressure_Update_Time_Begin);
        fprintf(logfile,"Pressure_Update_Time_End= %d\n",Pressure_Update_Time_End);
        fprintf(logfile,"Density_Comp_Flag= %d\n",Density_Comp_Flag);
        fprintf(logfile,"Density_Bin_Size= %f\n",Density_Bin_Size);
        fprintf(logfile,"Density_Update_Frequency= %d\n",Density_Update_Frequency);
        fprintf(logfile,"Density_Update_Time_Begin= %d\n",Density_Update_Time_Begin);
        fprintf(logfile,"Density_Update_Time_End= %d\n",Density_Update_Time_End);
        fprintf(logfile,"Density_Sigma_Gauss= %f\n",Density_Sigma_Gauss);
        fprintf(logfile,"Density_Gauss_Int_Range= %d\n",Density_Gauss_Int_Range);
        fprintf(logfile,"Density_Ref= %f\n",Density_Ref);
        fprintf(logfile,"Comp_Density_Scaling_factor_H = %f\n",Comp_Density_Scaling_factor_H);
        fprintf(logfile,"Load_File_Flag = %d\n",Load_File_Flag);
        fprintf(logfile,"Center_box = %f %f %f\n",center_box[0],center_box[1],center_box[2]);
        fprintf(logfile,"Density_Bin_Size = %f\n",Density_Bin_Size);
        fprintf(logfile,"Density_Bin_Num = %d\n",Density_Bin_Num);
        fprintf(logfile,"x0lo= %f\n",x0lo);
        fprintf(logfile,"x0hi= %f\n",x0hi);
        fprintf(logfile,"x0BoxSize= %f\n",x0BoxSize);
        fprintf(logfile,"d1= %f\n",R_Start_Hybrid_1);
        fprintf(logfile,"d2= %f\n",R_Start_Hybrid_2);
        fprintf(logfile,"moltype%d\n",atom->nmoltypesH);

    }
  }

  memory->create(Int_Mean_Density_H,Density_Bin_Num,atom->nmoltypesH+1,"lambdaH/calc:Int_Mean_Density_H");
  memory->create(Mean_Density_H,Density_Bin_Num,atom->nmoltypesH+1,"lambdaH/calc:Mean_Density_H");
  memory->create(Mean_Comp_Density_Conv_H,Density_Bin_Num,atom->nmoltypesH+1,"lambdaH/calc:Mean_Comp_Density_Conv_H");
  memory->create(grad_Comp_Density_Conv_H,Density_Bin_Num,atom->nmoltypesH+1,"lambdaH/calc:grad_Comp_Density_Conv_H");
  memory->create(Mean_grad_Comp_Density_Conv_H,Density_Bin_Num,atom->nmoltypesH+1,"lambdaH/calc:Mean_grad_Comp_Density_Conv_H");


  int This_Step = update->ntimestep;

  for(int i = 0;i < Density_Bin_Num; i++){
      for(int j = 0; j < atom->nmoltypesH; j++){
          Comp_Density_Num_H[i][j] = 0;
          Comp_Density_Num_all_H[i][j] = 0;
          Int_Mean_Density_H[i][j] = 0;
          Mean_Density_H[i][j] = 0;
          Mean_Comp_Density_Conv_H[i][j] = 0;
          Mean_grad_Comp_Density_Conv_H[i][j] = 0;
          grad_Comp_Density_Conv_H[i][j] = 0;
      }
  }

  if((This_Step >= Density_Update_Time_End || Load_File_Flag) && Density_Comp_Flag != 0)Load_Compensation_Density();

  if (atom->molecular == 0)
    error->all(FLERR,"Compute com/molecule requires molecular atom style");

  array_flag = 1;
  size_array_cols = 3;
  extarray = 0;
  // setup molecule-based data
  nmolecules = molecules_in_group(idlo,idhi);
  size_array_rows = nmolecules;

  //printf("Nmolecules= %d\n",nmolecules);

  memory->create(massproc_H,nmolecules,"lambdaH/calc:massproc_H");
  memory->create(masstotal_H,nmolecules,"lambdaH/calc:masstotal_H");
  memory->create(com_H,nmolecules,3,"lambdaH/calc:com_H");
  memory->create(comall_H,nmolecules,3,"lambdaH/calc:comall_H");
  memory->create(lambdaCM,nmolecules,"lambdaH/calc:lambdaCM");
  memory->create(gradlambdaCM,nmolecules,3,"lambdaH/calc:gradlambdaCM");

  // compute masstotal for each molecule

  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  tagint imol;
  double massone;

  for (int i = 0; i < nmolecules; i++) massproc_H[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap_H) imol = molmap_H[imol-idlo];
      else imol--;
      massproc_H[imol] += massone;
    }
}

  MPI_Allreduce(massproc_H,masstotal_H,nmolecules,MPI_DOUBLE,MPI_SUM,world);
  int *replambdaH = atom->replambdaH;
  int pass = 0;
  for (int i = 0; i < nlocal; i++) if(replambdaH[i] != 0)pass=1;
  if(pass==0)      error->all(FLERR,"No representative atom (replambdaH) has been defined!");

}

/* ---------------------------------------------------------------------- */

FixLambdaHCalc::~FixLambdaHCalc()
{

    memory->destroy(massproc_H);
    memory->destroy(masstotal_H);
    memory->destroy(com_H);
    memory->destroy(comall_H);
    memory->destroy(molmap_H);
    memory->destroy(lambdaCM);
    memory->destroy(gradlambdaCM);
    memory->destroy(Comp_Density_Num_H);
    memory->destroy(Comp_Density_Num_all_H);
    memory->destroy(Mean_Density_H);
    memory->destroy(Int_Mean_Density_H);
    memory->destroy(Mean_Comp_Density_Conv_H);
    memory->destroy(grad_Comp_Density_Conv_H);
    memory->destroy(Mean_grad_Comp_Density_Conv_H);
}

/* ---------------------------------------------------------------------- */

int FixLambdaHCalc::setmask()
{
  int mask = 0;
  //mask |= PRE_FORCE;
//  mask |= PRE_NEIGHBOR;
  mask |= POST_INTEGRATE;
  mask |= THERMO_ENERGY;
//  mask |= POST_FORCE_RESPA;
  mask |= MIN_PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLambdaHCalc::init()
{
    int ntmp = molecules_in_group(idlo,idhi);
    if (ntmp != nmolecules)
      error->all(FLERR,"Molecule count changed in compute com/molecule");



}

/* ---------------------------------------------------------------------- */

void FixLambdaHCalc::setup(int vflag)
{

//  if (strstr(update->integrate_style,"verlet"))
    //pre_force(vflag);
    post_integrate();
//  else {
//    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
//    post_force_respa(vflag,nlevels_respa-1,0);
 //   ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
 // }
}

/* ---------------------------------------------------------------------- */


void FixLambdaHCalc::post_integrate()
{

  tagint imol;
  double massone;
  double unwrap[3];


  //invoked_array = update->ntimestep;
  for (int i = 0; i < nmolecules; i++)
    com_H[i][0] = com_H[i][1] = com_H[i][2] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  imageint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double *lambdaH = atom->lambdaH;
  double **gradlambdaH = atom->gradlambdaH;
  double **comH = atom->comH;

  int This_Step = update->ntimestep;
  if(This_Step >= Density_Update_Time_Begin && This_Step < Density_Update_Time_End && Density_Comp_Flag == 1){
      Density_Compensation_Run = 1;
      if(me==0 && This_Step == Density_Update_Time_Begin){
          if(screen)fprintf(screen,"\nStart of constant-density route\n");
          if(logfile)fprintf(logfile,"\nStart of constant-density route\n");
          Clear_File_Compensation_Density();
      }
  }


  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        imol = molecule[i];
        if (molmap_H) imol = molmap_H[imol-idlo];
        else imol--;
        domain->unmap(x[i],image[i],unwrap);
        massone /= masstotal_H[imol];

        com_H[imol][0] += unwrap[0] * massone;
        com_H[imol][1] += unwrap[1] * massone;
        com_H[imol][2] += unwrap[2] * massone;
    }


  MPI_Allreduce(&com_H[0][0],&comall_H[0][0],3*nmolecules,MPI_DOUBLE,MPI_SUM,world);


  double xtmp,ytmp,ztmp;
  double delx,dely,delz;
  double Grad_Factor;
  double rdiff;

  for (int i = 0; i < nmolecules; i++){


      domain->remap(comall_H[i]);

      if(Hybrid_Style==0){
          xtmp = comall_H[i][0];

         if(xtmp < R_Start_Hybrid_1){
             lambdaCM[i] = 0.0;
             gradlambdaCM[i][0] = 0.0;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;

         }
         else if(xtmp >= R_Start_Hybrid_1 && xtmp < R_Start_Hybrid_1+Length_Hyb){
             cosx=cos(MY_PI*(xtmp-R_Start_Hybrid_1+Length_Hyb)/(2.0*Length_Hyb));
             sinx=sin(MY_PI*(xtmp-R_Start_Hybrid_1+Length_Hyb)/(2.0*Length_Hyb));
             lambdaCM[i] = cosx * cosx;
             gradlambdaCM[i][0] = -(MY_PI/Length_Hyb) * sinx * cosx;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;

         }
         else if(xtmp >=R_Start_Hybrid_1+Length_Hyb && xtmp < R_Start_Hybrid_2){
             lambdaCM[i] = 1.0;
             gradlambdaCM[i][0] = 0.0;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;

         }
         else if(xtmp >= R_Start_Hybrid_2 && xtmp < R_Start_Hybrid_2+Length_Hyb){
             cosx = cos(MY_PI*(xtmp-R_Start_Hybrid_2)/(2.0*Length_Hyb));
             sinx = sin(MY_PI*(xtmp-R_Start_Hybrid_2)/(2.0*Length_Hyb));
             lambdaCM[i] = cosx * cosx;
             gradlambdaCM[i][0] = -(MY_PI/Length_Hyb) * sinx * cosx;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;

         }
         else if(xtmp >=R_Start_Hybrid_2+Length_Hyb){
             lambdaCM[i] = 0.0;
             gradlambdaCM[i][0] = 0.0;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;

         }
    }

      else if(Hybrid_Style==1){
          xtmp = comall_H[i][0];
          ytmp = comall_H[i][1];
          ztmp = comall_H[i][2];

          delx = (xtmp-center_box[0]);
          dely = (ytmp-center_box[1]);
          delz = (ztmp-center_box[2]);

          r = sqrt(delx*delx+dely*dely+delz*delz);

         if(r < S_Start_Hybrid){
             lambdaCM[i] = 1.0;
             gradlambdaCM[i][0] = 0.0;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;
         }
         else if(r >= S_Start_Hybrid && r < S_Start_Hybrid+Length_Hyb) {
             rdiff = MY_PI*(r-S_Start_Hybrid)/(2.0*Length_Hyb);
             cosr=cos(rdiff);
             sinr=sin(rdiff);
             lambdaCM[i] = cosr * cosr;
             Grad_Factor = -MY_PI * sinr * cosr / (Length_Hyb * r);
             gradlambdaCM[i][0] = Grad_Factor * delx;
             gradlambdaCM[i][1] = Grad_Factor * dely;
             gradlambdaCM[i][2] = Grad_Factor * delz;
         }
         else if(r >= S_Start_Hybrid+Length_Hyb){
             lambdaCM[i] = 0.0;
             gradlambdaCM[i][0] = 0.0;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;
         }
    }

      else if(Hybrid_Style==2){

          xtmp = comall_H[i][0];
          ytmp = comall_H[i][1];

          delx = (xtmp-center_box[0]);
          dely = (ytmp-center_box[1]);

          r = sqrt(delx*delx+dely*dely);

         if(r < S_Start_Hybrid){
             lambdaCM[i] = 1.0;
             gradlambdaCM[i][0] = 0.0;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;
         }
         else if(r >= S_Start_Hybrid && r < S_Start_Hybrid+Length_Hyb) {
             rdiff = MY_PI*(r-S_Start_Hybrid)/(2.0*Length_Hyb);
             cosr=cos(rdiff);
             sinr=sin(rdiff);
             lambdaCM[i] = cosr * cosr;
             Grad_Factor = -MY_PI * sinr * cosr / (Length_Hyb * r);

             gradlambdaCM[i][0] = Grad_Factor * delx;
             gradlambdaCM[i][1] = Grad_Factor * dely;
             gradlambdaCM[i][2] = 0.0;
         }
         else if(r >= S_Start_Hybrid+Length_Hyb){
             lambdaCM[i] = 0.0;
             gradlambdaCM[i][0] = 0.0;
             gradlambdaCM[i][1] = 0.0;
             gradlambdaCM[i][2] = 0.0;
         }
    }


}

  int ibin,imoltypeH;
  int *moltypeH = atom->moltypeH;

  for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit)
        {
            imol = molecule[i];
            imoltypeH = moltypeH[i];
            if (molmap_H) imol = molmap_H[imol-idlo];
            else imol--;
            lambdaH[i] = lambdaCM[imol];
            gradlambdaH[i][0] = gradlambdaCM[imol][0];
            gradlambdaH[i][1] = gradlambdaCM[imol][1];
            gradlambdaH[i][2] = gradlambdaCM[imol][2];
//            if(gradlambdaH[i][2]!=0) printf("gradlambdaH[i][2] = %f\n",gradlambdaH[i][2]);

//            gradlambdaH[i][2] = 0.0;
            //gradlambdaH[i] = 0;
//              if(replambdaH[i] == 1){
                    comH[i][0] = comall_H[imol][0];
                    comH[i][1] = comall_H[imol][1];
                    comH[i][2] = comall_H[imol][2];
//              }

        if(Density_Comp_Flag != 0 && Density_Compensation_Run != 0){

            if(Hybrid_Style==0){
//                xtmp = comall_H[imol][0]-x0lo;
                xtmp = comH[i][0]-x0lo;
//                if(xtmp <= 0.5*x0BoxSize)ibin = floor(xtmp/Density_Bin_Size);
//                else ibin = floor((x0BoxSize-xtmp)/Density_Bin_Size);

                ibin = floor(xtmp/Density_Bin_Size);

            }

            else if(Hybrid_Style==1){
                xtmp = comall_H[imol][0];
                ytmp = comall_H[imol][1];
                ztmp = comall_H[imol][2];

                delx = (xtmp-center_box[0]);
                dely = (ytmp-center_box[1]);
                delz = (ztmp-center_box[2]);

                r = sqrt(delx*delx+dely*dely+delz*delz);
                ibin = floor(r/Density_Bin_Size);
            }

            else if(Hybrid_Style==2){

                xtmp = comall_H[imol][0];
                ytmp = comall_H[imol][1];

                delx = (xtmp-center_box[0]);
                dely = (ytmp-center_box[1]);

                r = sqrt(delx*delx+dely*dely);
                ibin = floor(r/Density_Bin_Size);
             }

            if(ibin>=Density_Bin_Num)ibin = Density_Bin_Num - 1;
            Comp_Density_Num_H[ibin][imoltypeH-1]++;

        }

        }


//  int Middle_Bin_Num = (Density_Bin_Num-1) / 2;

//  if(Hybrid_Style == 0)for(int i = Middle_Bin_Num+1;i < Density_Bin_Num; i++)Comp_Density_Num_H[i] = Comp_Density_Num_H[Density_Bin_Num-1-i];
  if(Density_Compensation_Run)Density_Counter_H++;

        double Density_Bin_Vol, exponent;
        double Normalization;
        int jmin, jmax,jj;
        if(This_Step % Density_Update_Frequency == 0 && This_Step > Density_Update_Time_Begin && Density_Comp_Flag != 0 && Density_Compensation_Run != 0){

      //   if(atom->nmoltypesH ==1)   MPI_Allreduce(&Comp_Density_Num_H[0][0],&Comp_Density_Num_all_H[0][0],Density_Bin_Num,MPI_INT,MPI_SUM,world);
        MPI_Allreduce(&Comp_Density_Num_H[0][0],&Comp_Density_Num_all_H[0][0],Density_Bin_Num*(atom->nmoltypesH+1),MPI_INT,MPI_SUM,world);

        if(Hybrid_Style == 0)Density_Bin_Vol = 1.0*Density_Bin_Size * (x1hi-x1lo) * (x2hi-x2lo);

        for(int i = 0; i < Density_Bin_Num; i++){
            if(Hybrid_Style == 1)Density_Bin_Vol = (4.0/3.0) * MY_PI * (pow(((i+1) * Density_Bin_Size),3.0) - pow((i * Density_Bin_Size),3.0));
            if(Hybrid_Style == 2)Density_Bin_Vol = 4.0 * MY_PI * (x2hi-x2lo) * (pow(((i+1) * Density_Bin_Size),2.0) - pow((i * Density_Bin_Size),2.0));
            for(int j = 0; j < atom->nmoltypesH; j++){
                Mean_Density_H[i][j] = Comp_Density_Num_all_H[i][j] / (Density_Counter_H * Density_Bin_Vol);
            }
        }

        Density_Counter_H = 0;

        for(int k = 0; k < atom->nmoltypesH;k++){
            for(int i = 0;i < Density_Bin_Num; i++){
                Normalization = 0;
                jmin = i - Density_Gauss_Int_Range*floor(Density_Sigma_Gauss / Density_Bin_Size);
                jmax = i + Density_Gauss_Int_Range*floor(Density_Sigma_Gauss / Density_Bin_Size);
                Mean_Comp_Density_Conv_H[i][k] = 0;
                    if(Hybrid_Style != 0){
                        if(jmin<0)jmin=0;
                        if(jmax>=Density_Bin_Num)jmax=Density_Bin_Num-1;
                    }

                for(int j = jmin;j <= jmax; j++){
                    jj = j;
                    if(Hybrid_Style == 0){
                            if(j < 0)jj = Density_Bin_Num + j;
                            if(j >= Density_Bin_Num)jj = j - Density_Bin_Num;
                    }
                    exponent = (i-j)*Density_Bin_Size/Density_Sigma_Gauss;
                    exponent = pow(exponent,2.0);
                    Normalization += exp(-exponent)*Density_Bin_Size;
                    Mean_Comp_Density_Conv_H[i][k] += Mean_Density_H[jj][k] * exp(-exponent)*Density_Bin_Size;
    //                    Mean_Comp_Density_Conv_H[i] += Mean_Density_H[j] * exp(-exponent);
                }
                Mean_Comp_Density_Conv_H[i][k] /= Normalization;

                }
            }

            Density_Fluctuation = 0;
            for(int k = 0; k < atom->nmoltypesH;k++)for(int i = 0;i < Density_Bin_Num; i++)Density_Fluctuation += (Density_Bin_Size/x0BoxSize)*pow((Mean_Density_H[i][k]-Density_Ref)/Density_Ref,2.0);
            for(int k = 0; k < atom->nmoltypesH;k++)
                for(int i = 1;i < Density_Bin_Num-1; i++)
                    grad_Comp_Density_Conv_H[i][k] = Comp_Density_Scaling_factor_H * (Mean_Comp_Density_Conv_H[i+1][k] - Mean_Comp_Density_Conv_H[i-1][k])/(2*Density_Bin_Size*Density_Ref);


            //for(int i = 1;i < Density_Bin_Num-1; i++)Mean_grad_Comp_Density_Conv_H[i] = (Comp_Counter_H * Mean_grad_Comp_Density_Conv_H[i] + grad_Comp_Density_Conv_H[i]) / (Comp_Counter_H + 1);
            for(int k = 0; k < atom->nmoltypesH;k++)
                for(int i = 1;i < Density_Bin_Num-1; i++)
                    Mean_grad_Comp_Density_Conv_H[i][k] +=  grad_Comp_Density_Conv_H[i][k];

//            for(int i = Middle_Bin_Num+1;i < Density_Bin_Num; i++)Comp_Density_Num_H[i] = Comp_Density_Num_H[Density_Bin_Num-1-i];

            for(int k = 0; k < atom->nmoltypesH;k++)
                for(int i = 0;i < Density_Bin_Num; i++)
                    Int_Mean_Density_H[i][k]=0;

            for(int k = 0; k < atom->nmoltypesH;k++)
                for(int i = 0;i < Density_Bin_Num; i++)
                    for(int j = 0;j <= i; j++)
                        Int_Mean_Density_H[i][k] += Mean_grad_Comp_Density_Conv_H[j][k] * Density_Bin_Size;

                        Comp_Counter_H += 1;

            for(int k = 0; k < atom->nmoltypesH;k++)
                for(int i = 0;i < Density_Bin_Num; i++){
                    Comp_Density_Num_H[i][k] = 0;
                    Comp_Density_Num_all_H[i][k] = 0;
                }


          }

        if(This_Step == Density_Update_Time_End && Density_Compensation_Run != 0){
            Density_Compensation_Run = 0;
            if(me==0){
                if(screen)fprintf(screen,"\nEnd of constant-density route\n");
                if(logfile)fprintf(logfile,"\nEnd of constant-density route\n");
            }
        }

        if (me == 0 && This_Step % Density_Update_Frequency == 0 && Density_Compensation_Run != 0) Print_Compensation_Density();

}


double FixLambdaHCalc::memory_usage()
{
  double bytes = (bigint) nmolecules * 2 * sizeof(double);
  if (molmap_H) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += (bigint) nmolecules * 2*3 * sizeof(double);
  bytes += (bigint) nmolecules * 2*3 * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */


int FixLambdaHCalc::molecules_in_group(tagint &idlo, tagint &idhi)
{
  int i;

  memory->destroy(molmap_H);
  molmap_H = NULL;

  // find lo/hi molecule ID for any atom in group
  // warn if atom in group has ID = 0

  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  tagint lo = BIG;
  tagint hi = -BIG;
  int flag = 0;
  for (i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit) {
      if (molecule[i] == 0) flag = 1;
      lo = MIN(lo,molecule[i]);
      hi = MAX(hi,molecule[i]);
    }
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"Atom with molecule ID = 0 included in "
                   "compute molecule group");

  MPI_Allreduce(&lo,&idlo,1,MPI_LMP_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_LMP_TAGINT,MPI_MAX,world);
  if (idlo == BIG) return 0;

  // molmap = vector of length nlen
  // set to 1 for IDs that appear in group across all procs, else 0

  tagint nlen_tag = idhi-idlo+1;
  if (nlen_tag > MAXSMALLINT)
    error->all(FLERR,"Too many molecules for compute");
  int nlen = (int) nlen_tag;

  memory->create(molmap_H,nlen,"lambdaH/calc:molmap");
  for (i = 0; i < nlen; i++) molmap_H[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      molmap_H[molecule[i]-idlo] = 1;

  int *molmapall;
  memory->create(molmapall,nlen,"compute:molmapall");
  MPI_Allreduce(molmap_H,molmapall,nlen,MPI_INT,MPI_MAX,world);

  // nmolecules = # of non-zero IDs in molmap
  // molmap[i] = index of molecule, skipping molecules not in group with -1

  int nmolecules = 0;
  for (i = 0; i < nlen; i++)
    if (molmapall[i]) molmap_H[i] = nmolecules++;
    else molmap_H[i] = -1;
  memory->destroy(molmapall);

  // warn if any molecule has some atoms in group and some not in group

  flag = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) continue;
    if (molecule[i] < idlo || molecule[i] > idhi) continue;
    if (molmap_H[molecule[i]-idlo] >= 0) flag = 1;
  }

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,
                   "One or more compute molecules has atoms not in group");

  // if molmap simply stores 1 to Nmolecules, then free it

  if (idlo == 1 && idhi == nmolecules && nlen == nmolecules) {
    memory->destroy(molmap_H);
    molmap_H = NULL;
  }
  return nmolecules;
}


void FixLambdaHCalc::Print_Compensation_Density(){


    FILE *fp1;
    fp1 = fopen("Mean_Comp_Density.txt","w");

    if (fp1 == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix Mean_Comp_Density file %s","Mean_Comp_Density.txt");
        error->one(FLERR,str);
    }


    for(int i = 0;i < Density_Bin_Num; i++){
        fprintf(fp1,"%d",i+1);
        for(int k = 0; k < atom->nmoltypesH;k++)
            fprintf(fp1,"\t%.10f",Mean_grad_Comp_Density_Conv_H[i][k]);
        fprintf(fp1,"\n");
    }
    fclose(fp1);

    /*
    FILE *fp2;

    char filename2[1000];
    sprintf(filename2, "Mean_Comp_Density%4.1f_%d.txt", Comp_Density_Scaling_factor_H,Density_Update_Frequency);

    fp2 = fopen(filename2,"a");

    if (fp2 == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix Mean_Comp_Density file %s","Mean_Comp_Density.txt");
        error->one(FLERR,str);
    }

    printf("Density_Fluctuation = %f\n",Density_Fluctuation);

    double x0;
    if(Hybrid_Style==0)x0=x0lo;
    else x0=0;

    for(int i = 0;i < Density_Bin_Num; i++){
        fprintf(fp2,"%d\t %.10f\t",i+1,i*Density_Bin_Size+x0);
        for(int k = 0; k < atom->nmoltypesH;k++)
            fprintf(fp2,"\t%.10f\t %.10f\t %.10f\t %.10f",Mean_Comp_Density_Conv_H[i][k],grad_Comp_Density_Conv_H[i][k],Mean_grad_Comp_Density_Conv_H[i][k]);
        fprintf(fp2,"\n");

    }
    fclose(fp2);


    FILE *fp3;
    char filename3[1000];
    sprintf(filename3, "Density_Fluctuation%4.1f_%d.txt", Comp_Density_Scaling_factor_H,Density_Update_Frequency);

    fp3 = fopen(filename3,"a");

    if (fp3 == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix Density_Fluctuation file %s","Density_Fluctuation.txt");
        error->one(FLERR,str);
    }

    int This_Step = update->ntimestep;
    fprintf(fp3,"%d %.10f\n",This_Step,Density_Fluctuation);

    fclose(fp3);

*/
}

void FixLambdaHCalc::Clear_File_Compensation_Density(){

    FILE *fp1;
    fp1 = fopen("Mean_Comp_Density.txt","w");

    if (fp1 == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix Mean_Comp_Density file %s","Mean_Comp_Density.txt");
        error->one(FLERR,str);
    }


    fclose(fp1);

    /*
    FILE *fp2;

    char filename2[1000];
    sprintf(filename2, "Mean_Comp_Density%4.1f_%d.txt", Comp_Density_Scaling_factor_H,Density_Update_Frequency);

    fp2 = fopen(filename2,"w");

    if (fp2 == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix Mean_Comp_Density file %s","Mean_Comp_Density.txt");
        error->one(FLERR,str);
    }


   fclose(fp2);

   FILE *fp3;
   char filename3[1000];
   sprintf(filename3, "Density_Fluctuation%4.1f_%d.txt", Comp_Density_Scaling_factor_H,Density_Update_Frequency);

   fp3 = fopen(filename3,"w");

   if (fp3 == NULL) {
       char str[128];
       sprintf(str,"Cannot open fix Density_Fluctuation file %s","Density_Fluctuation.txt");
       error->one(FLERR,str);
   }

   fclose(fp3);
    */

}

void FixLambdaHCalc::Load_Compensation_Density(){

    if(me == 0){
        FILE *fp1;
        char str[128];

        fp1 = fopen("Mean_Comp_Density.txt","r");
        if (fp1 == NULL) {
            sprintf(str,"Cannot open fix Mean_Comp_Density.txt file %s","Mean_Comp_Density.txt");
            error->one(FLERR,str);
        }

        int i1;
        float i2;

        while (!feof(fp1)){
            fscanf (fp1,"%d",&i1);
            for(int k = 0; k < atom->nmoltypesH;k++){
                fscanf (fp1,"\t%f",&i2);
                Mean_grad_Comp_Density_Conv_H[i1-1][k] = (double) i2;
            }
            fscanf (fp1,"\n");

                if(i1 > Density_Bin_Num){
                    sprintf(str,"Density bin number mismatches %d != %d",Density_Bin_Num,i1);
                    error->one(FLERR,str);
                }

        }

        fclose(fp1);


        if(me==0){
            if(screen)fprintf(screen,"\n\nDensity componsation forces distributed successfully!\n\n");
            if(logfile)fprintf(logfile,"\n\nDensity componsation forces distributed successfully!\n\n");
        }
    }

    MPI_Bcast(Mean_grad_Comp_Density_Conv_H,Density_Bin_Num*(atom->nmoltypesH+1),MPI_DOUBLE,0,world);

}
