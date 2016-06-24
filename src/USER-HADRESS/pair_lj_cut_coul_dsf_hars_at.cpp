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
   Contributing author: Maziar Heidari,
                        Robinson Cortes-Huerto,
                        Davide Donadio and
                        Raffaello Potestio
                        (Max Planck Institute for Polymer Research, 2015-2016)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut_coul_dsf_hars_at.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"
#include "fix_lambdah_calc.h"
#include "update.h"


using namespace LAMMPS_NS;
using namespace MathConst;
#define BIG MAXTAGINT


#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutCoulDSFHARSAT::PairLJCutCoulDSFHARSAT(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

  AT_molmap_H = NULL;
  nmolecules = molecules_in_group(idlo,idhi);

  H_AdResS_allocated = 0;

  AT_Pressure_Compensation_Run = 0;
  Comp_Counter_H = 0;
  memory->create(AT_massproc_H,nmolecules,"pair:AT_massproc_H");
  memory->create(AT_masstotal_H,nmolecules,"pair:AT_masstotal_H");
  memory->create(AT_mol_f_H,nmolecules,3,"pair:AT_mol_f_H");
  memory->create(AT_mol_f_all_H,nmolecules,3,"pair:AT_mol_f_all_H");


  // compute masstotal for each molecule
  MPI_Comm_rank(world, &me);
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  tagint imol;
  double massone;

  for (int i = 0; i < nmolecules; i++) AT_massproc_H[i] = 0.0;


  for (int i = 0; i < nlocal; i++)
  {
//    if (mask[i] & groupbit) {
    if (mask[i]) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (AT_molmap_H) imol = AT_molmap_H[imol-idlo];
      else imol--;
      AT_massproc_H[imol] += massone;

    }
  }

  MPI_Allreduce(AT_massproc_H,AT_masstotal_H,nmolecules,MPI_DOUBLE,MPI_SUM,world);

}

/* ---------------------------------------------------------------------- */

PairLJCutCoulDSFHARSAT::~PairLJCutCoulDSFHARSAT()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj,Vij_Coul,Vij_Lj;;
  double prefactor,erfcc,erfcd,t;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int imoltypeH,jmoltypeH;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *lambdaH = atom->lambdaH;
  double **gradlambdaH = atom->gradlambdaH;

  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  int imol,jmol;
  tagint *molecule = atom->molecule;
  double *mass = atom->mass;
  int *replambdaH = atom->replambdaH;
  int *moltypeH = atom->moltypeH;

  int ibin, jbin;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double xtmpj, iLambda, jLambda, ijLambda;

  int This_Step = update->ntimestep;
  if(This_Step >= AT_Update_Time_Begin && This_Step < AT_Update_Time_End && AT_Pressure_Comp_Flag != 0) AT_Pressure_Compensation_Run = 1;


  for (int i = 0; i < nmolecules; i++) {
      AT_mol_f_H[i][0] = AT_mol_f_H[i][1] = AT_mol_f_H[i][2] = 0;
      AT_mol_f_all_H[i][0] = AT_mol_f_all_H[i][1] = AT_mol_f_all_H[i][2] = 0;
  }

  // loop over neighbors of my atoms
//  if(update->ntimestep<AT_Restart_Time_Step+1)return;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    imoltypeH = moltypeH[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    iLambda = lambdaH[i];

    if (eflag) {
      double e_self = -(e_shift/2.0 + alpha/MY_PIS) * qtmp*qtmp*qqrd2e;
//      e_self = 0;
      ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      jLambda = lambdaH[j];

      if(iLambda==0 && jLambda==0 && AllAtomistic!=1)continue;


      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      jmoltypeH = moltypeH[j];

      imol = molecule[i];
      jmol = molecule[j];
      if (AT_molmap_H) {
          imol = AT_molmap_H[imol-idlo];
          jmol = AT_molmap_H[jmol-idlo];
      }
      else {
          imol--;
          jmol--;
      }

//      if(imol == jmol)continue;
      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        } else forcelj = 0.0;

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          prefactor = factor_coul * qqrd2e*qtmp*q[j]/r;
          erfcd = exp(-alpha*alpha*r*r);
          t = 1.0 / (1.0 + EWALD_P*alpha*r);
          erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
          forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd + 
            r*f_shift) * r;
        } else forcecoul = 0.0;

        if(((iLambda==1 && jLambda==1) || AllAtomistic)){

            fpair = (forcecoul + factor_lj*forcelj) * r2inv;
            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;

            if (newton_pair || j < nlocal) {
                f[j][0] -= delx*fpair;
                f[j][1] -= dely*fpair;
                f[j][2] -= delz*fpair;
            }

            if (eflag) {
            if (rsq < cut_ljsq[itype][jtype]) {
                evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                        offset[itype][jtype];
                evdwl *= factor_lj;
              } else evdwl = 0.0;
          
              if (rsq < cut_coulsq) {
                ecoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
              } else ecoul = 0.0;
            }
            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 evdwl,ecoul,fpair,delx,dely,delz);
        }

        else if(iLambda !=0 || jLambda != 0){

                ijLambda = 0.5 * (iLambda + jLambda);
                fpair = (forcecoul + factor_lj*forcelj) *ijLambda* r2inv;

                if (rsq < cut_coulsq)
                    Vij_Coul = 0.5*prefactor * (erfcc - r*e_shift - rsq*f_shift);
                else Vij_Coul = 0.0;


                if (rsq < cut_ljsq[itype][jtype])
                  Vij_Lj = 0.5*factor_lj*(r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                          offset[itype][jtype]);
                else Vij_Lj = 0;

                f[i][0] += delx*fpair;
                f[i][1] += dely*fpair;
                f[i][2] += delz*fpair;

                ibin = floor(iLambda/AT_lambda_Increment);
                if(ibin==AT_Bin_Num)ibin = AT_Bin_Num - 1;

                if(AT_Pressure_Compensation_Run != 0 && iLambda != 0 && iLambda != 1)Comp_Energy_H[ibin][imoltypeH-1] += Vij_Coul+Vij_Lj;

                AT_mol_f_H[imol][0] += -(Vij_Coul+Vij_Lj)*gradlambdaH[i][0];
                AT_mol_f_H[imol][1] += -(Vij_Coul+Vij_Lj)*gradlambdaH[i][1];
                AT_mol_f_H[imol][2] += -(Vij_Coul+Vij_Lj)*gradlambdaH[i][2];

                if (newton_pair || j < nlocal) {
                    f[j][0] -= delx*fpair;
                    f[j][1] -= dely*fpair;
                    f[j][2] -= delz*fpair;

                    jbin = floor(jLambda/AT_lambda_Increment);
                    if(jbin==AT_Bin_Num)jbin = AT_Bin_Num - 1;

                    if(AT_Pressure_Compensation_Run != 0 && jLambda != 0 && jLambda != 1)Comp_Energy_H[jbin][jmoltypeH-1] += Vij_Coul+Vij_Lj;

                    AT_mol_f_H[jmol][0] += -(Vij_Coul+Vij_Lj)*gradlambdaH[j][0];
                    AT_mol_f_H[jmol][1] += -(Vij_Coul+Vij_Lj)*gradlambdaH[j][1];
                    AT_mol_f_H[jmol][2] += -(Vij_Coul+Vij_Lj)*gradlambdaH[j][2];

                }

                if (eflag) {
                    ecoul = ijLambda*Vij_Coul*2.0;
                    evdwl = ijLambda*Vij_Lj*2.0;

                }
                if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                     evdwl,ecoul,fpair,delx,dely,delz);
            }



      }
    }
  }

    MPI_Allreduce(&AT_mol_f_H[0][0],&AT_mol_f_all_H[0][0],3*nmolecules,MPI_DOUBLE,MPI_SUM,world);

    if(AT_Pressure_Compensation_Run != 0 && AT_Pressure_Comp_Flag != 0){

        for (int i = 0; i < nlocal; i++){
            iLambda = lambdaH[i];
            if(replambdaH[i]!=0 && iLambda != 0 && iLambda != 1){
                ibin = floor(iLambda/AT_lambda_Increment);
                if(ibin==AT_Bin_Num)ibin = AT_Bin_Num - 1;
                imoltypeH = moltypeH[i]-1;
                Comp_Energy_Num_H[ibin][imoltypeH]++;
            }
        }

        if(This_Step % AT_Update_Frequency == 0 && This_Step > AT_Update_Time_Begin)AT_Update_Compensation_Energy();

    }


    if(This_Step == AT_Update_Time_End && AT_Pressure_Compensation_Run != 0)AT_Print_Compensation_Energy();

    double mol_mass,mass_frac;

    if(AllAtomistic != 1){
        if(AT_Pressure_Comp_Flag != 0){
            for (int i = 0; i < nlocal; i++){
                iLambda = lambdaH[i];
                if(iLambda != 0 && iLambda != 1){
                    imol = molecule[i];
                    if (AT_molmap_H)imol = AT_molmap_H[imol-idlo];
                    else imol--;
                    mol_mass = AT_masstotal_H[imol];
                    mass_frac = mass[type[i]] / mol_mass;
                    ibin = floor(iLambda/AT_lambda_Increment);
                    imoltypeH = moltypeH[i] - 1;

                    f[i][0] += mass_frac*(AT_mol_f_all_H[imol][0]+gradlambdaH[i][0]*Mean_Comp_Energy_H[ibin][imoltypeH]);
                    f[i][1] += mass_frac*(AT_mol_f_all_H[imol][1]+gradlambdaH[i][1]*Mean_Comp_Energy_H[ibin][imoltypeH]);
                    f[i][2] += mass_frac*(AT_mol_f_all_H[imol][2]+gradlambdaH[i][2]*Mean_Comp_Energy_H[ibin][imoltypeH]);
                    if (evflag) ev_tally(i,i,nlocal,newton_pair,
                                       -0.5*Int_Mean_Energy_H[ibin][imoltypeH],0.0,0.0,0.0,0.0,0.0);

                }

            }

        }
        else{

            for (int i = 0; i < nlocal; i++){

                iLambda = lambdaH[i];
                if(iLambda != 0 && iLambda != 1){

                    imol = molecule[i];
                    if (AT_molmap_H)imol = AT_molmap_H[imol-idlo];
                    else imol--;
                    mol_mass = AT_masstotal_H[imol];
                    mass_frac = mass[type[i]] / mol_mass;
                    ibin = floor(iLambda/AT_lambda_Increment);

                    f[i][0] += mass_frac*(AT_mol_f_all_H[imol][0]);
                    f[i][1] += mass_frac*(AT_mol_f_all_H[imol][1]);
                    f[i][2] += mass_frac*(AT_mol_f_all_H[imol][2]);

                }

            }

        }

        }



    if(This_Step == AT_Update_Time_End)AT_Pressure_Compensation_Run = 0;

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  alpha = force->numeric(FLERR,arg[0]);
  cut_lj_global = force->numeric(FLERR,arg[1]);
  cut_coul = force->numeric(FLERR,arg[2]);
  
  AllAtomistic = force->numeric(FLERR,arg[3]);
  Load_File_Flag = force->numeric(FLERR,arg[4]);
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) 
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);
  
  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double cut_lj_one = cut_lj_global;
  if (narg == 5) cut_lj_one = force->numeric(FLERR,arg[4]);
    
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;

      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/cut/coul/dsf requires atom attribute q");

  if(me == 0){
      if (screen)fprintf(screen,"AT_H_AdResS_allocated flag = %d\n",H_AdResS_allocated);
      if (logfile)fprintf(logfile,"AT_H_AdResS_allocated flag = %d\n",H_AdResS_allocated);
  }

  if(!H_AdResS_allocated)H_AdResS_Allocation();

  int This_Step = update->ntimestep;
  AT_Restart_Time_Step = This_Step;

  if((This_Step > AT_Update_Time_Begin || Load_File_Flag) && AT_Pressure_Comp_Flag != 0)Load_Compensation_Pressure();

  if(This_Step < AT_Update_Time_End && This_Step >= AT_Update_Time_Begin)Comp_Counter_H = floor((This_Step-AT_Update_Time_Begin)/AT_Update_Frequency);

  if(me==0 && This_Step < AT_Update_Time_End && This_Step > AT_Update_Time_Begin){
      if(screen)fprintf(screen,"AT_Pressure componsation forces are again being updated after previous %d times\n",Comp_Counter_H);
      if(logfile)fprintf(logfile,"AT_Pressure componsation forces are again being updated after previous %d times\n",Comp_Counter_H);
  }

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;
  double erfcc = erfc(alpha*cut_coul); 
  double erfcd = exp(-alpha*alpha*cut_coul*cut_coul);
  f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*alpha*erfcd/cut_coul); 
  e_shift = erfcc/cut_coul - f_shift*cut_coul; 

}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutCoulDSFHARSAT::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i],cut_lj[j][j]);
  }

  double cut = MAX(cut_lj[i][j],cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  
  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
     
  if (offset_flag) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;
  
  cut_ljsq[j][i] = cut_ljsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);
        
    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut_lj[i][j]*cut_lj[i][j]*cut_lj[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] * 
               sig6 * (sig6 - 3.0*rc6) / (9.0*rc9); 
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] * 
               sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9); 
  } 

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
	    }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::write_restart_settings(FILE *fp)
{
  fwrite(&alpha,sizeof(double),1,fp);
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulDSFHARSAT::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&alpha,sizeof(double),1,fp);
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&alpha,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulDSFHARSAT::single(int i, int j, int itype, int jtype, double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,r,erfcc,erfcd,prefactor;
  double forcecoul,forcelj,phicoul,philj;
  
  r2inv = 1.0/rsq;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  } else forcelj = 0.0;

  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = factor_coul * force->qqrd2e * atom->q[i]*atom->q[j]/r;
    erfcc = erfc(alpha*r); 
    erfcd = exp(-alpha*alpha*r*r);
    forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd + 
      r*f_shift) * r;
  } else forcecoul = 0.0;
  
  fforce = (forcecoul + factor_lj*forcelj) * r2inv;
      
  double eng = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
      offset[itype][jtype];
    eng += factor_lj*philj;
  }

  if (rsq < cut_coulsq) { 
    phicoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
    eng += phicoul;
  } 
  
  eng=0;
  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairLJCutCoulDSFHARSAT::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  return NULL;
}


int PairLJCutCoulDSFHARSAT::molecules_in_group(tagint &idlo, tagint &idhi)
{
  int i;

  memory->destroy(AT_molmap_H);
  AT_molmap_H = NULL;

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
//    if (mask[i] & groupbit) {
      if (mask[i]) {
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

  memory->create(AT_molmap_H,nlen,"pair:molmap_H");
  for (i = 0; i < nlen; i++) AT_molmap_H[i] = 0;

  for (i = 0; i < nlocal; i++)
   // if (mask[i] & groupbit)
       if (mask[i])
      AT_molmap_H[molecule[i]-idlo] = 1;

  int *AT_molmapall;
  memory->create(AT_molmapall,nlen,"pair:AT_molmapall");
  MPI_Allreduce(AT_molmap_H,AT_molmapall,nlen,MPI_INT,MPI_MAX,world);

  // nmolecules = # of non-zero IDs in molmap
  // molmap[i] = index of molecule, skipping molecules not in group with -1

  int nmolecules = 0;
  for (i = 0; i < nlen; i++)
    if (AT_molmapall[i]) AT_molmap_H[i] = nmolecules++;
    else AT_molmap_H[i] = -1;
  memory->destroy(AT_molmapall);

  // warn if any molecule has some atoms in group and some not in group

  flag = 0;
  for (i = 0; i < nlocal; i++) {
//    if (mask[i] & groupbit) continue;
      if (mask[i]) continue;
    if (molecule[i] < idlo || molecule[i] > idhi) continue;
    if (AT_molmap_H[molecule[i]-idlo] >= 0) flag = 1;
  }

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,
                   "One or more compute molecules has atoms not in group");

  // if molmap simply stores 1 to Nmolecules, then free it

  if (idlo == 1 && idhi == nmolecules && nlen == nmolecules) {
    memory->destroy(AT_molmap_H);
    AT_molmap_H = NULL;
  }
  return nmolecules;
}


void PairLJCutCoulDSFHARSAT::AT_Print_Compensation_Energy(){

    FILE *fp1;
    fp1 = fopen("Mean_Comp_Energy_AT.txt","w");
    if (fp1 == NULL) {
        char str[128];
        sprintf(str,"Cannot open Mean_Comp_Energy_AT.txt file %s","Mean_Comp_Energy_AT.txt");
        error->one(FLERR,str);
    }

    for(int i = 0;i < AT_Bin_Num; i++){
        fprintf(fp1,"%d",i+1);
            for(int j = 0; j < atom->nmoltypesH; j++){
                fprintf(fp1,"\t%.10f",Mean_Comp_Energy_H[i][j]);
            }
            fprintf(fp1,"\n");
    }
    fclose(fp1);
}



void PairLJCutCoulDSFHARSAT::AT_Update_Compensation_Energy(){
    MPI_Allreduce(&Comp_Energy_H[0][0],&Comp_Energy_all_H[0][0],AT_Bin_Num*(atom->nmoltypesH+1),MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&Comp_Energy_Num_H[0][0],&Comp_Energy_Num_all_H[0][0],AT_Bin_Num*(atom->nmoltypesH+1),MPI_INT,MPI_SUM,world);

    for(int j = 0;j < atom->nmoltypesH; j++){
        for(int i = 0;i < AT_Bin_Num; i++){
            Mean_Energy_H[i][j] = Comp_Energy_all_H[i][j] / Comp_Energy_Num_all_H[i][j];
            Mean_Comp_Energy_H[i][j] = (Comp_Counter_H * Mean_Comp_Energy_H[i][j] + Mean_Energy_H[i][j]) / (Comp_Counter_H + 1);
            Int_Mean_Energy_H[i][j]=0;
        }
    }

    Comp_Counter_H++;

    for(int j = 0;j < atom->nmoltypesH; j++){
        for(int i = 0;i < AT_Bin_Num; i++){
            Comp_Energy_Num_H[i][j] = 0;
            Comp_Energy_Num_all_H[i][j] = 0;
            Comp_Energy_H[i][j] = 0;
            Comp_Energy_all_H[i][j] =0;
        }

    }
    if (me == 0)AT_Print_Compensation_Energy();


}


void PairLJCutCoulDSFHARSAT::Load_Compensation_Pressure(){

    if(me == 0){
        FILE *fp1;
        char str[128];

        fp1 = fopen("Mean_Comp_Energy_AT.txt","r");
        if (fp1 == NULL) {
            sprintf(str,"Cannot open Mean_Comp_Energy_AT.txt file %s","Mean_Comp_Energy_AT.txt");
            error->one(FLERR,str);
        }

        int i1;
        float f1;

        while (!feof(fp1)){
                fscanf (fp1,"%d",&i1);
                for(int j=0;j<atom->nmoltypesH;j++){
                    fscanf (fp1,"\t%f",&f1);
                    Mean_Comp_Energy_H[i1-1][j] = f1;
                }
                if(i1 > AT_Bin_Num){
                    sprintf(str,"At drift force compensation bin number mismatches %d != %d",AT_Bin_Num,i1);
                    error->one(FLERR,str);
                }

        }

        fclose(fp1);


        if(me==0){
            if(screen)fprintf(screen,"AT_Pressure componsation forces distributed successfully!\n");
            if(logfile)fprintf(logfile,"AT_Pressure componsation forces distributed successfully!\n");
        }

    }

    MPI_Bcast(Mean_Comp_Energy_H,AT_Bin_Num,MPI_DOUBLE,0,world);


}

void PairLJCutCoulDSFHARSAT::H_AdResS_Allocation(){
    for (int i = 0; i < modify->nfix; i++){
      if (strcmp(modify->fix[i]->style,"lambdah/calc") == 0){
             lambda_H_fix = (FixLambdaHCalc *) modify->fix[i];
             AT_lambda_Increment = lambda_H_fix->Pressure_lambda_Increment;
             AT_Bin_Num = lambda_H_fix->Pressure_Bin_Num;
             AT_Update_Frequency = lambda_H_fix->Pressure_Update_Frequency;
             AT_Update_Time_Begin = lambda_H_fix->Pressure_Update_Time_Begin;
             AT_Update_Time_End = lambda_H_fix->Pressure_Update_Time_End;
             AT_Pressure_Comp_Flag = lambda_H_fix->Pressure_Comp_Flag;
             AT_center_box = lambda_H_fix->center_box;
             AT_Hybrid_Style = lambda_H_fix->Hybrid_Style;
      }
    }

    if(me == 0){
        if (screen){
            fprintf(screen,"AT_lambda_Increment= %f\n",AT_lambda_Increment);
            fprintf(screen,"AT_Bin_Num= %d\n",AT_Bin_Num);
            fprintf(screen,"AT_Update_Frequency= %d\n",AT_Update_Frequency);
            fprintf(screen,"AT_Update_Time_Begin= %d\n",AT_Update_Time_Begin);
            fprintf(screen,"AT_Update_Time_End= %d\n",AT_Update_Time_End);
            fprintf(screen,"AT_Pressure_Comp_Flag= %d\n",AT_Pressure_Comp_Flag);
        }

        if (logfile){
            fprintf(logfile,"AT_lambda_Increment= %f\n",AT_lambda_Increment);
            fprintf(logfile,"AT_Bin_Num= %d\n",AT_Bin_Num);
            fprintf(logfile,"AT_Update_Frequency= %d\n",AT_Update_Frequency);
            fprintf(logfile,"AT_Update_Time_Begin= %d\n",AT_Update_Time_Begin);
            fprintf(logfile,"AT_Update_Time_End= %d\n",AT_Update_Time_End);
            fprintf(logfile,"AT_Pressure_Comp_Flag= %d\n",AT_Pressure_Comp_Flag);
        }
    }

    memory->create(Comp_Energy_Num_H,AT_Bin_Num,atom->nmoltypesH+1,"pairLJHAT:Comp_Energy_Num_H");
    memory->create(Comp_Energy_Num_all_H,AT_Bin_Num,atom->nmoltypesH+1,"pairLJHAT:Comp_Energy_Num_all_H");

    memory->create(Int_Mean_Energy_H,AT_Bin_Num,atom->nmoltypesH+1,"pairLJHAT:Int_Mean_Energy_H");
    memory->create(Comp_Energy_H,AT_Bin_Num,atom->nmoltypesH+1,"pairLJHAT:Comp_Energy_H");
    memory->create(Comp_Energy_all_H,AT_Bin_Num,atom->nmoltypesH+1,"pairLJHAT:Comp_Energy_all_H");
    memory->create(Mean_Comp_Energy_H,AT_Bin_Num,atom->nmoltypesH+1,"pairLJHAT:Mean_Comp_Energy_H");
    memory->create(Mean_Energy_H,AT_Bin_Num,atom->nmoltypesH+1,"pairLJHAT:Mean_Energy_H");


    for(int j=0;j<atom->nmoltypesH;j++){
        for(int i = 0;i < AT_Bin_Num; i++){
            Int_Mean_Energy_H[i][j]=0;
            Comp_Energy_H[i][j]=0;
            Comp_Energy_all_H[i][j]=0;
            Mean_Comp_Energy_H[i][j]=0;
            Comp_Energy_Num_H[i][j] = 0;
            Comp_Energy_Num_all_H[i][j] = 0;
        }

    }

    H_AdResS_allocated = 1;
}

