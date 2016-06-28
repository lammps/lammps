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
#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut_hars_cg.h"
#include "atom.h"
#include "atom_vec.h"
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
#include "iostream"
#include "fix.h"
#include "fix_lambdah_calc.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define BIG MAXTAGINT

/* ---------------------------------------------------------------------- */

PairLJCutHARSCG::PairLJCutHARSCG(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;

  massproc_H = NULL;
  masstotal_H = NULL;
  molmap_H = NULL;
  mol_f_H = NULL;
  mol_f_all_H = NULL;
  Comp_Energy_Num_H = NULL;
  Comp_Energy_Num_all_H = NULL;

  Int_Mean_Energy_H = NULL;
  Mean_Energy_H = NULL;
  Comp_Energy_H = NULL;
  Comp_Energy_all_H = NULL;
  Mean_Comp_Energy_H = NULL;
  CG_Mean_grad_Comp_Density_Conv_H = NULL;

  molmap_H = NULL;
  nmolecules = molecules_in_group(idlo,idhi);

  H_AdResS_allocated = 0;

  CG_Pressure_Compensation_Run = 0;
  Density_Compensation_Run = 0;
  Comp_Counter_H = 0;
  CG_Density_Comp_Flag = 0;
  CG_Pressure_Comp_Flag = 0;


  memory->create(massproc_H,nmolecules,"pair:massproc_H");
  memory->create(masstotal_H,nmolecules,"pair:masstotal_H");
  memory->create(mol_f_H,nmolecules,3,"pair:mol_f_H");
  memory->create(mol_f_all_H,nmolecules,3,"pair:mol_f_all_H");


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

  for (int i = 0; i < nmolecules; i++) massproc_H[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
  {
//    if (mask[i] & groupbit) {
    if (mask[i]) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap_H) imol = molmap_H[imol-idlo];
      else imol--;
      massproc_H[imol] += massone;

    }
  }

  MPI_Allreduce(massproc_H,masstotal_H,nmolecules,MPI_DOUBLE,MPI_SUM,world);

}

/* ---------------------------------------------------------------------- */

PairLJCutHARSCG::~PairLJCutHARSCG()
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
      memory->destroy(massproc_H);
      memory->destroy(masstotal_H);
      memory->destroy(molmap_H);
      memory->destroy(mol_f_H);
      memory->destroy(mol_f_all_H);
      memory->destroy(Comp_Energy_Num_H);
      memory->destroy(Comp_Energy_Num_all_H);
      memory->destroy(Int_Mean_Energy_H);
      memory->destroy(Mean_Energy_H);
      memory->destroy(Comp_Energy_H);
      memory->destroy(Comp_Energy_all_H);
      memory->destroy(Mean_Comp_Energy_H);
      //    memory->destroy(CG_Mean_grad_Comp_Density_Conv_H);

//    delete lambda_H_fix;
}
}

/* ---------------------------------------------------------------------- */

void PairLJCutHARSCG::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj, factor_lj,Vij;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int imoltype,jmoltype;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *lambdaH = atom->lambdaH;
  double **gradlambdaH = atom->gradlambdaH;
  double **comH = atom->comH;
  int *replambdaH = atom->replambdaH;
  tagint *molecule = atom->molecule;
  double *mass = atom->mass;
  int *moltypeH = atom->moltypeH;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  int ibin, jbin;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double iLambda, jLambda, ijLambda;
  int imol,jmol;



  int This_Step = update->ntimestep;
  if(This_Step >= CG_Update_Time_Begin && This_Step < CG_Update_Time_End && CG_Pressure_Comp_Flag != 0){
      CG_Pressure_Compensation_Run = 1;
      if(me==0 && This_Step == CG_Update_Time_Begin){
          if(screen)fprintf(screen,"\nStart of constant-pressure route\n");
          if(logfile)fprintf(logfile,"\nStart of constant-pressure route\n");
      }
  }


  if(update->ntimestep<CG_Restart_Time_Step+1)return;


  for (int i = 0; i < nmolecules; i++) {
      mol_f_H[i][0] = mol_f_H[i][1] = mol_f_H[i][2] = 0;
      mol_f_all_H[i][0] = mol_f_all_H[i][1] = mol_f_all_H[i][2] = 0;
  }


  Density_Compensation_Run = lambda_H_fix->Density_Compensation_Run;

  if(Density_Compensation_Run){
      CG_Mean_grad_Comp_Density_Conv_H = lambda_H_fix->Mean_grad_Comp_Density_Conv_H;
  }

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    if(replambdaH[i] == 0)continue;

    xtmp = comH[i][0];
    ytmp = comH[i][1];
    ztmp = comH[i][2];

//    itype = moltypeH[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    imoltype = itype - 1;

    iLambda = 1 - lambdaH[i];

    for (jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jLambda = 1 - lambdaH[j];

      if(replambdaH[j] == 0)continue;

      if((iLambda==0 && jLambda==0) && AllCoarseGrained != 1)continue;

      delx = xtmp - comH[j][0];
      dely = ytmp - comH[j][1];
      delz = ztmp - comH[j][2];

      domain->minimum_image(delx,dely,delz);

      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      //jtype = moltypeH[j];
      jmoltype = jtype - 1;

//      if (rsq < cutsq[itype][jtype] && lj1[itype][jtype] != 0) {
      if (rsq < cutsq[itype][jtype]) {


          imol = molecule[i];
          jmol = molecule[j];
          if (molmap_H) {
              imol = molmap_H[imol-idlo];
              jmol = molmap_H[jmol-idlo];
          }
          else {
              imol--;
              jmol--;
          }


         if(((iLambda==1 && jLambda==1) || AllCoarseGrained)){

                r2inv = 1.0/rsq;
                r6inv = r2inv*r2inv*r2inv;
                forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);

                fpair = factor_lj*forcelj*r2inv;
                mol_f_H[imol][0] += delx*fpair;
                mol_f_H[imol][1] += dely*fpair;
                mol_f_H[imol][2] += delz*fpair;

                if (newton_pair || j < nlocal) {
                  mol_f_H[jmol][0] -= delx*fpair;
                  mol_f_H[jmol][1] -= dely*fpair;
                  mol_f_H[jmol][2] -= delz*fpair;
                }

                if (eflag) {
                    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                            offset[itype][jtype];
                  evdwl *= factor_lj;
                }

                if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                     evdwl,0.0,fpair,delx,dely,delz);

            }

            else if(iLambda != 0 || jLambda != 0){

                ijLambda = 0.5 * (iLambda + jLambda);
                r2inv = 1.0/rsq;
                r6inv = r2inv*r2inv*r2inv;
                forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
                fpair = factor_lj*(forcelj*ijLambda)*r2inv;

                Vij = 0.5*(r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                        offset[itype][jtype]);


                ibin = floor(iLambda/CG_lambda_Increment);
                if(ibin==CG_Bin_Num)ibin = CG_Bin_Num - 1;

                if(CG_Pressure_Compensation_Run != 0 && iLambda != 0 && iLambda != 1)Comp_Energy_H[ibin][imoltype] += Vij;

                mol_f_H[imol][0] += delx*fpair + Vij*gradlambdaH[i][0];
                mol_f_H[imol][1] += dely*fpair + Vij*gradlambdaH[i][1];
                mol_f_H[imol][2] += delz*fpair + Vij*gradlambdaH[i][2];

                if (newton_pair || j < nlocal) {

                  jbin = floor(jLambda/CG_lambda_Increment);
                  if(jbin==CG_Bin_Num)jbin = CG_Bin_Num - 1;
                  if(CG_Pressure_Compensation_Run != 0 && jLambda != 0 && jLambda != 1)Comp_Energy_H[jbin][jmoltype] += Vij;

                  mol_f_H[jmol][0] -= delx*fpair - Vij*gradlambdaH[j][0];
                  mol_f_H[jmol][1] -= dely*fpair - Vij*gradlambdaH[j][1];
                  mol_f_H[jmol][2] -= delz*fpair - Vij*gradlambdaH[j][2];

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


  MPI_Allreduce(&mol_f_H[0][0],&mol_f_all_H[0][0],3*nmolecules,MPI_DOUBLE,MPI_SUM,world);

  if(CG_Pressure_Compensation_Run != 0 && CG_Pressure_Comp_Flag != 0){
      for (int i = 0; i < nlocal; i++){
          iLambda = 1 - lambdaH[i];
          if(replambdaH[i] != 0 && lambdaH[i] != 0 && lambdaH[i] != 1){
              ibin = floor(iLambda/CG_lambda_Increment);
              if(ibin==CG_Bin_Num)ibin = CG_Bin_Num - 1;
              itype = moltypeH[i] - 1;
              Comp_Energy_Num_H[ibin][itype]++;

          }
      }

      if(This_Step % CG_Update_Frequency == 0 && This_Step > CG_Update_Time_Begin)CG_Update_Compensation_Energy();
  }


  double mol_mass,mass_frac;
  double Grad_Density,r;
  if(AllCoarseGrained != 1 && (CG_Density_Comp_Flag != 0 || CG_Pressure_Comp_Flag != 0)){

      for (int i = 0; i < nlocal; i++){
            imol = molecule[i];
            if (molmap_H)imol = molmap_H[imol-idlo];
            else imol--;
            mol_mass = masstotal_H[imol];
            mass_frac = mass[type[i]] / mol_mass;
            iLambda = 1 - lambdaH[i];

            if(iLambda != 0 && iLambda != 1){

                ibin = floor(iLambda/CG_lambda_Increment);
                itype = moltypeH[i] - 1;

                f[i][0] += mass_frac*(mol_f_all_H[imol][0]-gradlambdaH[i][0]*Mean_Comp_Energy_H[ibin][itype]);
                f[i][1] += mass_frac*(mol_f_all_H[imol][1]-gradlambdaH[i][1]*Mean_Comp_Energy_H[ibin][itype]);
                f[i][2] += mass_frac*(mol_f_all_H[imol][2]-gradlambdaH[i][2]*Mean_Comp_Energy_H[ibin][itype]);

                if (evflag) ev_tally(i,i,nlocal,newton_pair,
                                     -0.5*Int_Mean_Energy_H[ibin][itype],0.0,0.0,0.0,0.0,0.0);


                if(CG_Density_Comp_Flag != 0){

                    if(CG_Hybrid_Style == 0){
                        ibin = floor((comH[i][0]-CG_x0lo)/CG_Density_Bin_Size);
                        f[i][0] += mass_frac*(-1.0*CG_Mean_grad_Comp_Density_Conv_H[ibin][itype]);
                    }
                    else if(CG_Hybrid_Style == 1){
                        delx = comH[i][0] - CG_center_box[0];
                        dely = comH[i][1] - CG_center_box[1];
                        delz = comH[i][2] - CG_center_box[2];
                        r = sqrt(delx*delx + dely*dely + delz*delz);
                        ibin = floor(r/CG_Density_Bin_Size);

                        Grad_Density = CG_Mean_grad_Comp_Density_Conv_H[ibin][itype] / r;
                        f[i][0] += -mass_frac * Grad_Density * delx;
                        f[i][1] += -mass_frac * Grad_Density * dely;
                        f[i][2] += -mass_frac * Grad_Density * delz;
                    }
                    else if(CG_Hybrid_Style == 2){
                        delx = comH[i][0] - CG_center_box[0];
                        dely = comH[i][1] - CG_center_box[1];
                        r = sqrt(delx*delx + dely*dely);
                        ibin = floor(r/CG_Density_Bin_Size);

                        Grad_Density = CG_Mean_grad_Comp_Density_Conv_H[ibin][itype] / r;
                        f[i][0] += -mass_frac * Grad_Density * delx;
                        f[i][1] += -mass_frac * Grad_Density * dely;
                    }

                }
            }
            else{
                f[i][0] += mass_frac*mol_f_all_H[imol][0];
                f[i][1] += mass_frac*mol_f_all_H[imol][1];
                f[i][2] += mass_frac*mol_f_all_H[imol][2];
            }

     }

  }
  else{

      for (int i = 0; i < nlocal; i++){
            imol = molecule[i];
            if (molmap_H)imol = molmap_H[imol-idlo];
            else imol--;
            mol_mass = masstotal_H[imol];
            mass_frac = mass[type[i]] / mol_mass;

            f[i][0] += mass_frac*mol_f_all_H[imol][0];
            f[i][1] += mass_frac*mol_f_all_H[imol][1];
            f[i][2] += mass_frac*mol_f_all_H[imol][2];
      }

  }



  if(This_Step == CG_Update_Time_End){
      CG_Pressure_Compensation_Run = 0;
      if(me == 0){
          if(screen)fprintf(screen,"\nEnd of constant-pressure route\n");
          if(logfile)fprintf(logfile,"\nEnd of constant-pressure route\n");

      }
  }


  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLJCutHARSCG::compute_inner()
{

    error->all(FLERR,"Rrespa has not been included!");

}

/* ---------------------------------------------------------------------- */

void PairLJCutHARSCG::compute_middle()
{
    error->all(FLERR,"Rrespa has not been included!");

}

/* ---------------------------------------------------------------------- */

void PairLJCutHARSCG::compute_outer(int eflag, int vflag)
{
    error->all(FLERR,"Rrespa has not been included!");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutHARSCG::allocate()
{

  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pairLJHCG:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pairLJHCG:cutsq");
  memory->create(cut,n+1,n+1,"pairLJHCG:cut");
  memory->create(epsilon,n+1,n+1,"pairLJHCG:epsilon");
  memory->create(sigma,n+1,n+1,"pairLJHCG:sigma");
  memory->create(lj1,n+1,n+1,"pairLJHCG:lj1");
  memory->create(lj2,n+1,n+1,"pairLJHCG:lj2");
  memory->create(lj3,n+1,n+1,"pairLJHCG:lj3");
  memory->create(lj4,n+1,n+1,"pairLJHCG:lj4");
  memory->create(offset,n+1,n+1,"pairLJHCG:offset");



}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutHARSCG::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  AllCoarseGrained  = force->numeric(FLERR,arg[1]);

  Load_File_Flag = force->numeric(FLERR,arg[2]);
  // reset cutoffs that have been explicitly set

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

void PairLJCutHARSCG::coeff(int narg, char **arg)
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

void PairLJCutHARSCG::init_style()
{
  // request regular or rRESPA neighbor lists
    if(me == 0){
        if (screen)fprintf(screen,"CG_H_AdResS_allocated flag = %d\n",H_AdResS_allocated);
        if (logfile)fprintf(logfile,"CG_H_AdResS_allocated flag = %d\n",H_AdResS_allocated);
    }

    if(!H_AdResS_allocated)H_AdResS_Allocation();


    int This_Step = update->ntimestep;
    CG_Restart_Time_Step = This_Step;

    if((This_Step > CG_Update_Time_Begin || Load_File_Flag)  && CG_Pressure_Comp_Flag != 0)Load_Compensation_Pressure();

    if(This_Step < CG_Update_Time_End && This_Step >= CG_Update_Time_Begin)Comp_Counter_H = floor((This_Step-CG_Update_Time_Begin)/CG_Update_Frequency);

    if(me==0 && This_Step < CG_Update_Time_End && This_Step > CG_Update_Time_Begin){
        if(screen)fprintf(screen,"CG_Pressure componsation forces are again being updated after previous %d times\n",Comp_Counter_H);
        if(logfile)fprintf(logfile,"CG_Pressure componsation forces are again being updated after previous %d times\n",Comp_Counter_H);
    }

  int irequest;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0)irequest = neighbor->request(this,instance_me);
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

  } else irequest = neighbor->request(this,instance_me);

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

void PairLJCutHARSCG::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutHARSCG::init_one(int i, int j)
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

void PairLJCutHARSCG::write_restart(FILE *fp)
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

void PairLJCutHARSCG::read_restart(FILE *fp)
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

void PairLJCutHARSCG::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutHARSCG::read_restart_settings(FILE *fp)
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

void PairLJCutHARSCG::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCutHARSCG::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJCutHARSCG::single(int i, int j, int itype, int jtype, double rsq,
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

void *PairLJCutHARSCG::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return NULL;
}



int PairLJCutHARSCG::molecules_in_group(tagint &idlo, tagint &idhi)
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

  memory->create(molmap_H,nlen,"pair:molmap_H");
  for (i = 0; i < nlen; i++) molmap_H[i] = 0;

  for (i = 0; i < nlocal; i++)
   // if (mask[i] & groupbit)
       if (mask[i])
      molmap_H[molecule[i]-idlo] = 1;

  int *molmapall;
  memory->create(molmapall,nlen,"pair:molmapall");
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
//    if (mask[i] & groupbit) continue;
      if (mask[i]) continue;
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



void PairLJCutHARSCG::CG_Print_Compensation_Energy(){

    FILE *fp1;

    fp1 = fopen("Mean_Comp_Energy_CG.txt","w");
    if (fp1 == NULL) {
        char str[128];
        sprintf(str,"Cannot open Mean_Comp_Energy_CG.txt file %s","Mean_Comp_Energy_CG.txt");
        error->one(FLERR,str);

    }

    for(int i = 0;i < CG_Bin_Num; i++){
        fprintf(fp1,"%d",i+1);
        for(int k = 0; k < atom->nmoltypesH;k++)
            fprintf(fp1,"\t%.10f",Mean_Comp_Energy_H[i][k]);

        fprintf(fp1,"\n");
    }

    fclose(fp1);

}



void PairLJCutHARSCG::CG_Update_Compensation_Energy(){

    MPI_Allreduce(&Comp_Energy_H[0][0],&Comp_Energy_all_H[0][0],CG_Bin_Num*(atom->nmoltypesH+1),MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&Comp_Energy_Num_H[0][0],&Comp_Energy_Num_all_H[0][0],CG_Bin_Num*(atom->nmoltypesH+1),MPI_INT,MPI_SUM,world);

    for(int k = 0; k < atom->nmoltypesH;k++)for(int i = 0;i < CG_Bin_Num; i++)Mean_Energy_H[i][k] = Comp_Energy_all_H[i][k] / Comp_Energy_Num_all_H[i][k];
    for(int k = 0; k < atom->nmoltypesH;k++)for(int i = 0;i < CG_Bin_Num; i++)Mean_Comp_Energy_H[i][k] = (Comp_Counter_H * Mean_Comp_Energy_H[i][k] + Mean_Energy_H[i][k]) / (Comp_Counter_H + 1);
    for(int k = 0; k < atom->nmoltypesH;k++)for(int i = 0;i < CG_Bin_Num; i++)Int_Mean_Energy_H[i][k]=0;
    for(int k = 0; k < atom->nmoltypesH;k++)for(int i = 0;i < CG_Bin_Num; i++)for(int j = 0;j <= i; j++)Int_Mean_Energy_H[i][k] += Mean_Comp_Energy_H[j][k] * CG_lambda_Increment;

    Comp_Counter_H++;

    for(int k = 0; k < atom->nmoltypesH;k++){
        for(int i = 0;i < CG_Bin_Num; i++){
            Comp_Energy_Num_H[i][k] = 0;
            Comp_Energy_Num_all_H[i][k] = 0;
            Comp_Energy_H[i][k] = 0;
            Comp_Energy_all_H[i][k] =0;
        }
    }

    if (me == 0)CG_Print_Compensation_Energy();


}


void PairLJCutHARSCG::Load_Compensation_Pressure(){

    if(me == 0){
        FILE *fp1;
        char str[128];

        fp1 = fopen("Mean_Comp_Energy_CG.txt","r");
        if (fp1 == NULL) {
            sprintf(str,"Cannot open fix Mean_Comp_Energy_CG.txt file %s","Mean_Comp_Energy_CG.txt");
            error->one(FLERR,str);
        }

        int i1;
        float i2;

        while (!feof(fp1)){
            fscanf (fp1,"%d",&i1);
            for(int k = 0; k < atom->nmoltypesH;k++){
                fscanf (fp1,"\t%f",&i2);
                Mean_Comp_Energy_H[i1-1][k] = i2;
                if(i1 > CG_Bin_Num){
                    sprintf(str,"CG drift force compensation bin number mismatches %d != %d",CG_Bin_Num,i1);
                    error->one(FLERR,str);
                }
            }

        }
    }


    if(me==0){
        if(screen)fprintf(screen,"CG_Pressure componsation forces distributed successfully!\n");
        if(logfile)fprintf(logfile,"CG_Pressure componsation forces distributed successfully!\n");
    }

    MPI_Bcast(Mean_Comp_Energy_H,CG_Bin_Num*(atom->nmoltypesH+1),MPI_DOUBLE,0,world);
}


void PairLJCutHARSCG::H_AdResS_Allocation(){


    for (int i = 0; i < modify->nfix; i++){

        if (strcmp(modify->fix[i]->style,"lambdah/calc") == 0){

            lambda_H_fix = (FixLambdaHCalc *) modify->fix[i];
            CG_lambda_Increment = lambda_H_fix->Pressure_lambda_Increment;
            CG_Bin_Num = lambda_H_fix->Pressure_Bin_Num;
            CG_Update_Frequency = lambda_H_fix->Pressure_Update_Frequency;
            CG_Update_Time_Begin = lambda_H_fix->Pressure_Update_Time_Begin;
            CG_Update_Time_End = lambda_H_fix->Pressure_Update_Time_End;

            CG_Density_Bin_Num = lambda_H_fix->Density_Bin_Num;
            CG_Density_Bin_Size = lambda_H_fix->Density_Bin_Size;
            CG_Density_Update_Frequency = lambda_H_fix->Density_Update_Frequency;
            CG_Density_Update_Time_Begin = lambda_H_fix->Density_Update_Time_Begin;
            CG_Density_Update_Time_End = lambda_H_fix->Density_Update_Time_End;

            CG_Pressure_Comp_Flag = lambda_H_fix->Pressure_Comp_Flag;
            CG_Density_Comp_Flag = lambda_H_fix->Density_Comp_Flag;
            CG_center_box = lambda_H_fix->center_box;
            CG_Hybrid_Style = lambda_H_fix->Hybrid_Style;
            CG_x0lo = lambda_H_fix->x0lo;
            CG_x0BoxSize = lambda_H_fix->x0BoxSize;
        }
    }

    if(me == 0){
        if (screen){
            fprintf(screen,"CG_lambda_Increment= %f\n",CG_lambda_Increment);
            fprintf(screen,"CG_Bin_Num= %d\n",CG_Bin_Num);
            fprintf(screen,"CG_Update_Frequency= %d\n",CG_Update_Frequency);
            fprintf(screen,"CG_Update_Time_Begin= %d\n",CG_Update_Time_Begin);
            fprintf(screen,"CG_Update_Time_End= %d\n",CG_Update_Time_End);
            fprintf(screen,"CG_Pressure_Comp_Flag= %d\n",CG_Pressure_Comp_Flag);
            fprintf(screen,"CG_Density_Comp_Flag= %d\n",CG_Density_Comp_Flag);
            fprintf(screen,"CG_Hybrid_Style= %d\n",CG_Hybrid_Style);

        }

        if (logfile){
            fprintf(logfile,"CG_lambda_Increment= %f\n",CG_lambda_Increment);
            fprintf(logfile,"CG_Bin_Num= %d\n",CG_Bin_Num);
            fprintf(logfile,"CG_Update_Frequency= %d\n",CG_Update_Frequency);
            fprintf(logfile,"CG_Update_Time_Begin= %d\n",CG_Update_Time_Begin);
            fprintf(logfile,"CG_Update_Time_End= %d\n",CG_Update_Time_End);
            fprintf(logfile,"CG_Pressure_Comp_Flag= %d\n",CG_Pressure_Comp_Flag);
            fprintf(logfile,"CG_Density_Comp_Flag= %d\n",CG_Density_Comp_Flag);
            fprintf(logfile,"CG_Hybrid_Style= %d\n",CG_Hybrid_Style);

        }


    }


    memory->create(Comp_Energy_Num_H,CG_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:Comp_Energy_Num_H");
    memory->create(Comp_Energy_Num_all_H,CG_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:Comp_Energy_Num_all_H");
    memory->create(Int_Mean_Energy_H,CG_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:Int_Mean_Energy_H");
    memory->create(Comp_Energy_H,CG_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:Comp_Energy_H");
    memory->create(Comp_Energy_all_H,CG_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:Comp_Energy_all_H");
    memory->create(Mean_Comp_Energy_H,CG_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:Mean_Comp_Energy_H");
    memory->create(Mean_Energy_H,CG_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:Mean_Energy_H");

    memory->create(CG_Mean_grad_Comp_Density_Conv_H,CG_Density_Bin_Num,atom->nmoltypesH+1,"pairLJHCG:CG_Mean_grad_Comp_Density_Conv_H");
    CG_Mean_grad_Comp_Density_Conv_H = lambda_H_fix->Mean_grad_Comp_Density_Conv_H;


    for(int i = 0;i < CG_Bin_Num; i++){
        for(int j = 0; j < atom->nmoltypesH; j++){
            Int_Mean_Energy_H[i][j] = 0;
            Comp_Energy_H[i][j] = 0;
            Comp_Energy_all_H[i][j] = 0;
            Mean_Comp_Energy_H[i][j] = 0;
            Comp_Energy_Num_H[i][j] = 0;
            Comp_Energy_Num_all_H[i][j] = 0;
        }
    }


    H_AdResS_allocated = 1;

}
