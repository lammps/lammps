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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut_hars_at.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "fix.h"
#include "fix_lambdah_calc.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;
#define BIG MAXTAGINT

/* ---------------------------------------------------------------------- */

PairLJCutHARSAT::PairLJCutHARSAT(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;

  AT_massproc_H = NULL;
  AT_masstotal_H = NULL;
  AT_molmap_H = NULL;
  AT_mol_f_H = NULL;
  AT_mol_f_all_H = NULL;
  Comp_Energy_Num_H = NULL;
  Comp_Energy_Num_all_H = NULL;

  Int_Mean_Energy_H = NULL;
  Mean_Energy_H = NULL;
  Comp_Energy_H = NULL;
  Comp_Energy_all_H = NULL;
  Mean_Comp_Energy_H = NULL;

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

PairLJCutHARSAT::~PairLJCutHARSAT()
{
  if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);
      memory->destroy(cut);
      memory->destroy(epsilon);
      memory->destroy(sigma);
      memory->destroy(lj1);
      memory->destroy(lj2);
      memory->destroy(lj3);
      memory->destroy(lj4);
      memory->destroy(offset);
      memory->destroy(AT_massproc_H);
      memory->destroy(AT_masstotal_H);
      memory->destroy(AT_molmap_H);
      memory->destroy(AT_mol_f_H);
      memory->destroy(AT_mol_f_all_H);
      memory->destroy(Comp_Energy_Num_H);
      memory->destroy(Comp_Energy_Num_all_H);
      memory->destroy(Int_Mean_Energy_H);
      memory->destroy(Mean_Energy_H);
      memory->destroy(Comp_Energy_H);
      memory->destroy(Comp_Energy_all_H);
      memory->destroy(Mean_Comp_Energy_H);

  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutHARSAT::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj, factor_lj,Vij;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int imoltypeH,jmoltypeH;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *lambdaH = atom->lambdaH;
  double **gradlambdaH = atom->gradlambdaH;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
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

  //if(update->ntimestep<3)return;

  double xtmpj, iLambda, jLambda, ijLambda;

  int This_Step = update->ntimestep;
  if(This_Step >= AT_Update_Time_Begin && This_Step < AT_Update_Time_End && AT_Pressure_Comp_Flag != 0) AT_Pressure_Compensation_Run = 1;


  for (int i = 0; i < nmolecules; i++) {
      AT_mol_f_H[i][0] = AT_mol_f_H[i][1] = AT_mol_f_H[i][2] = 0;
      AT_mol_f_all_H[i][0] = AT_mol_f_all_H[i][1] = AT_mol_f_all_H[i][2] = 0;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    imoltypeH = moltypeH[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    iLambda = lambdaH[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jLambda = lambdaH[j];


      if(iLambda==0 && jLambda==0 && AllAtomistic!=1)continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      xtmpj = x[j][0];

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


      if (rsq < cutsq[itype][jtype] && lj1[itype][jtype] != 0 && (imol != jmol)) {


          if(((iLambda==1 && jLambda==1) || AllAtomistic)){

                r2inv = 1.0/rsq;
                r6inv = r2inv*r2inv*r2inv;
                forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
                fpair = factor_lj*forcelj*r2inv;

                f[i][0] += delx*fpair;
                f[i][1] += dely*fpair;
                f[i][2] += delz*fpair;

                if (newton_pair || j < nlocal) {
                  f[j][0] -= delx*fpair;
                  f[j][1] -= dely*fpair;
                  f[j][2] -= delz*fpair;
                }

                if (eflag) {
                    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                            offset[itype][jtype];
                    evdwl *= factor_lj;
                }

                if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                     evdwl,0.0,fpair,delx,dely,delz);

            }
            else if(iLambda !=0 || jLambda != 0){


                ijLambda = 0.5 * (iLambda + jLambda);
                r2inv = 1.0/rsq;
                r6inv = r2inv*r2inv*r2inv;
                forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);

                fpair = factor_lj*(forcelj*ijLambda)*r2inv;

                Vij = 0.5*(r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                        offset[itype][jtype]);

                f[i][0] += delx*fpair;
                f[i][1] += dely*fpair;
                f[i][2] += delz*fpair;

                ibin = floor(iLambda/AT_lambda_Increment);
                if(ibin==AT_Bin_Num)ibin = AT_Bin_Num - 1;
                if(AT_Pressure_Compensation_Run != 0 && iLambda != 0 && iLambda != 1)Comp_Energy_H[ibin][imoltypeH-1] += Vij;

                AT_mol_f_H[imol][0] += -Vij*gradlambdaH[i][0];
                AT_mol_f_H[imol][1] += -Vij*gradlambdaH[i][1];
                AT_mol_f_H[imol][2] += -Vij*gradlambdaH[i][2];

                if (newton_pair || j < nlocal) {

                  f[j][0] -= delx*fpair;
                  f[j][1] -= dely*fpair;
                  f[j][2] -= delz*fpair;

                  jbin = floor(jLambda/AT_lambda_Increment);
                  if(jbin==AT_Bin_Num)jbin = AT_Bin_Num - 1;
                  if(AT_Pressure_Compensation_Run != 0 && jLambda != 0 && jLambda != 1)Comp_Energy_H[jbin][jmoltypeH-1] += Vij;

                  AT_mol_f_H[jmol][0] += -Vij*gradlambdaH[j][0];
                  AT_mol_f_H[jmol][1] += -Vij*gradlambdaH[j][1];
                  AT_mol_f_H[jmol][2] += -Vij*gradlambdaH[j][2];

                }

                if (eflag) {

                    evdwl = ijLambda*Vij*2.0;
                    evdwl *= factor_lj;
                }

                if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fpair,delx,dely,delz);

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

/* ---------------------------------------------------------------------- */

void PairLJCutHARSAT::compute_inner()
{
    error->all(FLERR,"Rrespa has not been included!");
}

/* ---------------------------------------------------------------------- */

void PairLJCutHARSAT::compute_middle()
{
    error->all(FLERR,"Rrespa has not been included!");
}

/* ---------------------------------------------------------------------- */

void PairLJCutHARSAT::compute_outer(int eflag, int vflag)
{
    error->all(FLERR,"Rrespa has not been included!");

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutHARSAT::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pairLJHAT:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pairLJHAT:cutsq");

  memory->create(cut,n+1,n+1,"pairLJHAT:cut");
  memory->create(epsilon,n+1,n+1,"pairLJHAT:epsilon");
  memory->create(sigma,n+1,n+1,"pairLJHAT:sigma");
  memory->create(lj1,n+1,n+1,"pairLJHAT:lj1");
  memory->create(lj2,n+1,n+1,"pairLJHAT:lj2");
  memory->create(lj3,n+1,n+1,"pairLJHAT:lj3");
  memory->create(lj4,n+1,n+1,"pairLJHAT:lj4");
  memory->create(offset,n+1,n+1,"pairLJHAT:offset");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutHARSAT::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set
  AllAtomistic = force->numeric(FLERR,arg[1]);

  Load_File_Flag = force->numeric(FLERR,arg[2]);


  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }



}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutHARSAT::coeff(int narg, char **arg)
{


  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutHARSAT::init_style()
{
  // request regular or rRESPA neighbor lists


    if(me == 0){
        if (screen)fprintf(screen,"AT_H_AdResS_allocated flag = %d\n",H_AdResS_allocated);
        if (logfile)fprintf(logfile,"AT_H_AdResS_allocated flag = %d\n",H_AdResS_allocated);
    }

    if(!H_AdResS_allocated)H_AdResS_Allocation();

    int This_Step = update->ntimestep;
    if((This_Step > AT_Update_Time_Begin  || Load_File_Flag) && AT_Pressure_Comp_Flag != 0)Load_Compensation_Pressure();

    if(This_Step < AT_Update_Time_End && This_Step >= AT_Update_Time_Begin)Comp_Counter_H = floor((This_Step-AT_Update_Time_Begin)/AT_Update_Frequency);

    if(me==0 && This_Step < AT_Update_Time_End && This_Step > AT_Update_Time_Begin){
        if(screen)fprintf(screen,"AT_Pressure componsation forces are again being updated after previous %d times\n",Comp_Counter_H);
        if(logfile)fprintf(logfile,"AT_Pressure componsation forces are again being updated after previous %d times\n",Comp_Counter_H);
    }


    int irequest;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this,instance_me);
    else if (respa == 1) {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else {
      irequest = neighbor->request(this,instance_me);
  }
  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairLJCutHARSAT::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutHARSAT::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }


  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);


  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

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
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }


  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutHARSAT::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutHARSAT::read_restart(FILE *fp)
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
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutHARSAT::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutHARSAT::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCutHARSAT::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCutHARSAT::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJCutHARSAT::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r2inv,r6inv,forcelj,philj;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  fforce = factor_lj*forcelj*r2inv;

  philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
    offset[itype][jtype];
  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJCutHARSAT::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return NULL;
}


int PairLJCutHARSAT::molecules_in_group(tagint &idlo, tagint &idhi)
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


void PairLJCutHARSAT::AT_Print_Compensation_Energy(){

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



void PairLJCutHARSAT::AT_Update_Compensation_Energy(){

    MPI_Allreduce(&Comp_Energy_H[0][0],&Comp_Energy_all_H[0][0],AT_Bin_Num*atom->nmoltypesH,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&Comp_Energy_Num_H[0][0],&Comp_Energy_Num_all_H[0][0],AT_Bin_Num*atom->nmoltypesH,MPI_INT,MPI_SUM,world);

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


void PairLJCutHARSAT::Load_Compensation_Pressure(){

    if(me == 0){
        FILE *fp1;
        char str[128];

        fp1 = fopen("Mean_Comp_Energy_AT.txt","r");
        if (fp1 == NULL) {
            sprintf(str,"Cannot open fix Mean_Comp_Energy_AT.txt file %s","Mean_Comp_Energy_AT.txt");
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
                fscanf (fp1,"\n");
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

void PairLJCutHARSAT::H_AdResS_Allocation(){

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

    memory->create(Comp_Energy_Num_H,AT_Bin_Num,atom->nmoltypesH,"pairLJHAT:Comp_Energy_Num_H");
    memory->create(Comp_Energy_Num_all_H,AT_Bin_Num,atom->nmoltypesH,"pairLJHAT:Comp_Energy_Num_all_H");

    memory->create(Int_Mean_Energy_H,AT_Bin_Num,atom->nmoltypesH,"pairLJHAT:Int_Mean_Energy_H");
    memory->create(Comp_Energy_H,AT_Bin_Num,atom->nmoltypesH,"pairLJHAT:Comp_Energy_H");
    memory->create(Comp_Energy_all_H,AT_Bin_Num,atom->nmoltypesH,"pairLJHAT:Comp_Energy_all_H");
    memory->create(Mean_Comp_Energy_H,AT_Bin_Num,atom->nmoltypesH,"pairLJHAT:Mean_Comp_Energy_H");
    memory->create(Mean_Energy_H,AT_Bin_Num,atom->nmoltypesH,"pairLJHAT:Mean_Energy_H");

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

