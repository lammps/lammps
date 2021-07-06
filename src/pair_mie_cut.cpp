// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Cassiano Aimoli (aimoli@gmail.com)
------------------------------------------------------------------------- */

#include "pair_mie_cut.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairMIECut::PairMIECut(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

PairMIECut::~PairMIECut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(gamR);
    memory->destroy(gamA);
    memory->destroy(Cmie);
    memory->destroy(mie1);
    memory->destroy(mie2);
    memory->destroy(mie3);
    memory->destroy(mie4);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairMIECut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,rgamR,rgamA,forcemie,factor_mie;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_mie = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_mie = special_mie[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rgamA = pow(r2inv,(gamA[itype][jtype]/2.0));
        rgamR = pow(r2inv,(gamR[itype][jtype]/2.0));
        forcemie =  (mie1[itype][jtype]*rgamR - mie2[itype][jtype]*rgamA);
        fpair = factor_mie*forcemie*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = (mie3[itype][jtype]*rgamR - mie4[itype][jtype]*rgamA) -
            offset[itype][jtype];
          evdwl *= factor_mie;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairMIECut::compute_inner()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,r2inv,rgamR,rgamA,forcemie,factor_mie,rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_mie = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum_inner;
  ilist = list->ilist_inner;
  numneigh = list->numneigh_inner;
  firstneigh = list->firstneigh_inner;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];

  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_mie = special_mie[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_out_off_sq) {
        jtype = type[j];
        r2inv = 1.0/rsq;
        rgamA = pow(r2inv,(gamA[itype][jtype]/2.0));
        rgamR = pow(r2inv,(gamR[itype][jtype]/2.0));
        forcemie =  (mie1[itype][jtype]*rgamR - mie2[itype][jtype]*rgamA);
        fpair = factor_mie*forcemie*r2inv;
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair *= 1.0 - rsw*rsw*(3.0 - 2.0*rsw);
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMIECut::compute_middle()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,r2inv,rgamR,rgamA,forcemie,factor_mie,rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_mie = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum_middle;
  ilist = list->ilist_middle;
  numneigh = list->numneigh_middle;
  firstneigh = list->firstneigh_middle;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_mie = special_mie[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
        jtype = type[j];
        r2inv = 1.0/rsq;
        rgamA = pow(r2inv,(gamA[itype][jtype]/2.0));
        rgamR = pow(r2inv,(gamR[itype][jtype]/2.0));
        forcemie =  (mie1[itype][jtype]*rgamR - mie2[itype][jtype]*rgamA);
        fpair = factor_mie*forcemie*r2inv;
        if (rsq < cut_in_on_sq) {
          rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
          fpair *= rsw*rsw*(3.0 - 2.0*rsw);
        }
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair *= 1.0 + rsw*rsw*(2.0*rsw - 3.0);
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMIECut::compute_outer(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,rgamR,rgamA,forcemie,factor_mie,rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_mie = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_mie = special_mie[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        if (rsq > cut_in_off_sq) {
        r2inv = 1.0/rsq;
        rgamA = pow(r2inv,(gamA[itype][jtype]/2.0));
        rgamR = pow(r2inv,(gamR[itype][jtype]/2.0));
        forcemie =  (mie1[itype][jtype]*rgamR - mie2[itype][jtype]*rgamA);
          fpair = factor_mie*forcemie*r2inv;
          if (rsq < cut_in_on_sq) {
            rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
            fpair *= rsw*rsw*(3.0 - 2.0*rsw);
          }

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
        }

        if (eflag) {
          r2inv = 1.0/rsq;
        rgamA = pow(r2inv,(gamA[itype][jtype]/2.0));
        rgamR = pow(r2inv,(gamR[itype][jtype]/2.0));
          evdwl = (mie3[itype][jtype]*rgamR - mie4[itype][jtype]*rgamA) -
            offset[itype][jtype];
          evdwl *= factor_mie;
        }

        if (vflag) {
          if (rsq <= cut_in_off_sq) {
        r2inv = 1.0/rsq;
        rgamA = pow(r2inv,(gamA[itype][jtype]/2.0));
        rgamR = pow(r2inv,(gamR[itype][jtype]/2.0));
        forcemie =  (mie1[itype][jtype]*rgamR - mie2[itype][jtype]*rgamA);
            fpair = factor_mie*forcemie*r2inv;
          } else if (rsq < cut_in_on_sq)
            fpair = factor_mie*forcemie*r2inv;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMIECut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(gamR,n+1,n+1,"pair:gamR");
  memory->create(gamA,n+1,n+1,"pair:gamA");
  memory->create(Cmie,n+1,n+1,"pair:Cmie");
  memory->create(mie1,n+1,n+1,"pair:mie1");
  memory->create(mie2,n+1,n+1,"pair:mie2");
  memory->create(mie3,n+1,n+1,"pair:mie3");
  memory->create(mie4,n+1,n+1,"pair:mie4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMIECut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMIECut::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);
  double gamR_one = utils::numeric(FLERR,arg[4],false,lmp);
  double gamA_one = utils::numeric(FLERR,arg[5],false,lmp);

  double cut_one = cut_global;
  if (narg == 7) cut_one = utils::numeric(FLERR,arg[6],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      gamR[i][j] = gamR_one;
      gamA[i][j] = gamA_one;
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

void PairMIECut::init_style()
{
  // request regular or rRESPA neighbor list

  int irequest;
  int respa = 0;

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
  }

  irequest = neighbor->request(this,instance_me);

  if (respa >= 1) {
    neighbor->requests[irequest]->respaouter = 1;
    neighbor->requests[irequest]->respainner = 1;
  }
  if (respa == 2) neighbor->requests[irequest]->respamiddle = 1;

  // set rRESPA cutoffs

  if (utils::strmatch(update->integrate_style,"^respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = nullptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMIECut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    gamR[i][j] = mix_distance(gamR[i][i],gamR[j][j]);
    gamA[i][j] = mix_distance(gamA[i][i],gamA[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  gamA[j][i] = gamA[i][j];
  gamR[j][i] = gamR[i][j];
  Cmie[i][j] = (gamR[i][j]/(gamR[i][j]-gamA[i][j]) *
                pow((gamR[i][j]/gamA[i][j]),
                    (gamA[i][j]/(gamR[i][j]-gamA[i][j]))));
  mie1[i][j] = Cmie[i][j] * gamR[i][j]* epsilon[i][j] *
    pow(sigma[i][j],gamR[i][j]);
  mie2[i][j] = Cmie[i][j] * gamA[i][j] * epsilon[i][j] *
    pow(sigma[i][j],gamA[i][j]);
  mie3[i][j] = Cmie[i][j] * epsilon[i][j] * pow(sigma[i][j],gamR[i][j]);
  mie4[i][j] = Cmie[i][j] * epsilon[i][j] * pow(sigma[i][j],gamA[i][j]);

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = Cmie[i][j] * epsilon[i][j] *
      (pow(ratio,gamR[i][j]) - pow(ratio,gamA[i][j]));
  } else offset[i][j] = 0.0;

  mie1[j][i] = mie1[i][j];
  mie2[j][i] = mie2[i][j];
  mie3[j][i] = mie3[i][j];
  mie4[j][i] = mie4[i][j];
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

    double siggamA = pow(sigma[i][j],gamA[i][j]);
    double siggamR = pow(sigma[i][j],gamR[i][j]);
    double rcgamA = pow(cut[i][j],(gamA[i][j]-3.0));
    double rcgamR = pow(cut[i][j],(gamR[i][j]-3.0));
    etail_ij = Cmie[i][j]*2.0*MY_PI*all[0]*all[1]*epsilon[i][j]*
      (siggamR/((gamR[i][j]-3.0)*rcgamR)-siggamA/((gamA[i][j]-3.0)*rcgamA));
    ptail_ij = Cmie[i][j]*2.0*MY_PI*all[0]*all[1]*epsilon[i][j]/3.0*
      ((gamR[i][j]/(gamR[i][j]-3.0))*siggamR/rcgamR-
       (gamA[i][j]/(gamA[i][j]-3.0))*siggamA/rcgamA);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMIECut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&gamR[i][j],sizeof(double),1,fp);
        fwrite(&gamA[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMIECut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamR[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamA[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamR[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamA[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMIECut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMIECut::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairMIECut::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                           double /*factor_coul*/, double factor_mie,
                           double &fforce)
{
  double r2inv,rgamR,rgamA,forcemie,phimie;

  r2inv = 1.0/rsq;
  rgamA = pow(r2inv,(gamA[itype][jtype]/2.0));
  rgamR = pow(r2inv,(gamR[itype][jtype]/2.0));
  forcemie =  (mie1[itype][jtype]*rgamR - mie2[itype][jtype]*rgamA);
  fforce = factor_mie*forcemie*r2inv;

  phimie = (mie3[itype][jtype]*rgamR - mie4[itype][jtype]*rgamA) -
    offset[itype][jtype];
  return factor_mie*phimie;
}

/* ---------------------------------------------------------------------- */

void *PairMIECut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"gamR") == 0) return (void *) gamR;
  if (strcmp(str,"gamA") == 0) return (void *) gamA;
  return nullptr;
}
