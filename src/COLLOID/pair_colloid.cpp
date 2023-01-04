// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include "pair_colloid.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

PairColloid::PairColloid(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairColloid::~PairColloid()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(form);
    memory->destroy(a12);
    memory->destroy(sigma);
    memory->destroy(d1);
    memory->destroy(d2);
    memory->destroy(a1);
    memory->destroy(a2);
    memory->destroy(diameter);
    memory->destroy(cut);
    memory->destroy(offset);
    memory->destroy(sigma3);
    memory->destroy(sigma6);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
  }
}

/* ---------------------------------------------------------------------- */

void PairColloid::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r,forcelj,factor_lj;
  double r2inv,r6inv,c1,c2,fR,dUR,dUA;
  double K[9],h[4],g[4];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
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
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      switch (form[itype][jtype]) {
      case SMALL_SMALL:
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;
        if (eflag) evdwl = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
                     offset[itype][jtype];
        break;

      case SMALL_LARGE:
        c2 = a2[itype][jtype];
        K[1] = c2*c2;
        K[2] = rsq;
        K[0] = K[1] - rsq;
        K[4] = rsq*rsq;
        K[3] = K[1] - K[2];
        K[3] *= K[3]*K[3];
        K[6] = K[3]*K[3];
        fR = sigma3[itype][jtype]*a12[itype][jtype]*c2*K[1]/K[3];
        fpair = 4.0/15.0*fR*factor_lj *
          (2.0*(K[1]+K[2]) * (K[1]*(5.0*K[1]+22.0*K[2])+5.0*K[4]) *
           sigma6[itype][jtype]/K[6]-5.0) / K[0];
        if (eflag)
          evdwl = 2.0/9.0*fR *
            (1.0-(K[1]*(K[1]*(K[1]/3.0+3.0*K[2])+4.2*K[4])+K[2]*K[4]) *
             sigma6[itype][jtype]/K[6]) - offset[itype][jtype];
        if (rsq <= K[1])
          error->one(FLERR,"Overlapping small/large in pair colloid");
        break;

      case LARGE_LARGE:
        r = sqrt(rsq);
        c1 = a1[itype][jtype];
        c2 = a2[itype][jtype];
        K[0] = c1*c2;
        K[1] = c1+c2;
        K[2] = c1-c2;
        K[3] = K[1]+r;
        K[4] = K[1]-r;
        K[5] = K[2]+r;
        K[6] = K[2]-r;
        K[7] = 1.0/(K[3]*K[4]);
        K[8] = 1.0/(K[5]*K[6]);
        g[0] = powint(K[3],-7);
        g[1] = powint(K[4],-7);
        g[2] = powint(K[5],-7);
        g[3] = powint(K[6],-7);
        h[0] = ((K[3]+5.0*K[1])*K[3]+30.0*K[0])*g[0];
        h[1] = ((K[4]+5.0*K[1])*K[4]+30.0*K[0])*g[1];
        h[2] = ((K[5]+5.0*K[2])*K[5]-30.0*K[0])*g[2];
        h[3] = ((K[6]+5.0*K[2])*K[6]-30.0*K[0])*g[3];
        g[0] *= 42.0*K[0]/K[3]+6.0*K[1]+K[3];
        g[1] *= 42.0*K[0]/K[4]+6.0*K[1]+K[4];
        g[2] *= -42.0*K[0]/K[5]+6.0*K[2]+K[5];
        g[3] *= -42.0*K[0]/K[6]+6.0*K[2]+K[6];

        fR = a12[itype][jtype]*sigma6[itype][jtype]/r/37800.0;
        evdwl = fR * (h[0]-h[1]-h[2]+h[3]);
        dUR = evdwl/r + 5.0*fR*(g[0]+g[1]-g[2]-g[3]);
        dUA = -a12[itype][jtype]/3.0*r*((2.0*K[0]*K[7]+1.0)*K[7] +
                                        (2.0*K[0]*K[8]-1.0)*K[8]);
        fpair = factor_lj * (dUR+dUA)/r;
        if (eflag)
          evdwl += a12[itype][jtype]/6.0 *
            (2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7])) - offset[itype][jtype];
        if (r <= K[1])
          error->one(FLERR,"Overlapping large/large in pair colloid");
        break;
      }

      if (eflag) evdwl *= factor_lj;

      f[i][0] += delx*fpair;
      f[i][1] += dely*fpair;
      f[i][2] += delz*fpair;
      if (newton_pair || j < nlocal) {
        f[j][0] -= delx*fpair;
        f[j][1] -= dely*fpair;
        f[j][2] -= delz*fpair;
      }

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairColloid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(form,n+1,n+1,"pair:form");
  memory->create(a12,n+1,n+1,"pair:a12");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(d1,n+1,n+1,"pair:d1");
  memory->create(d2,n+1,n+1,"pair:d2");
  memory->create(a1,n+1,n+1,"pair:a1");
  memory->create(a2,n+1,n+1,"pair:a2");
  memory->create(diameter,n+1,n+1,"pair:diameter");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(offset,n+1,n+1,"pair:offset");
  memory->create(sigma3,n+1,n+1,"pair:sigma3");
  memory->create(sigma6,n+1,n+1,"pair:sigma6");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairColloid::settings(int narg, char **arg)
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

void PairColloid::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a12_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);
  double d1_one = utils::numeric(FLERR,arg[4],false,lmp);
  double d2_one = utils::numeric(FLERR,arg[5],false,lmp);

  double cut_one = cut_global;
  if (narg == 7) cut_one = utils::numeric(FLERR,arg[6],false,lmp);

  if (d1_one < 0.0 || d2_one < 0.0)
    error->all(FLERR,"Invalid d1 or d2 value for pair colloid coeff");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a12[i][j] = a12_one;
      sigma[i][j] = sigma_one;
      if (i == j && d1_one != d2_one)
        error->all(FLERR,"Invalid d1 or d2 value for pair colloid coeff");
      d1[i][j] = d1_one;
      d2[i][j] = d2_one;
      diameter[i][j] = 0.5*(d1_one+d2_one);
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairColloid::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    a12[i][j] = mix_energy(a12[i][i],a12[j][j],sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    d1[i][j] = mix_distance(d1[i][i],d1[j][j]);
    d2[i][j] = mix_distance(d2[i][i],d2[j][j]);
    diameter[i][j] = 0.5 * (d1[i][j] + d2[i][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  sigma3[i][j] = sigma[i][j]*sigma[i][j]*sigma[i][j];
  sigma6[i][j] = sigma3[i][j]*sigma3[i][j];

  if (d1[i][j] == 0.0 && d2[i][j] == 0.0) form[i][j] = SMALL_SMALL;
  else if (d1[i][j] == 0.0 || d2[i][j] == 0.0) form[i][j] = SMALL_LARGE;
  else form[i][j] = LARGE_LARGE;

  // for SMALL_SMALL, a1/a2 do not need to be set
  // for SMALL_LARGE, a1 does not need to be set, a2 = radius for i,j and j,i
  // for LARGE_LARGE, a1/a2 are radii, swap them for j,i

  if (form[i][j] == SMALL_LARGE) {
    if (d1[i][j] > 0.0) a2[i][j] = 0.5*d1[i][j];
    else a2[i][j] = 0.5*d2[i][j];
    a2[j][i] = a2[i][j];
  } else if (form[i][j] == LARGE_LARGE) {
    a2[j][i] = a1[i][j] = 0.5*d1[i][j];
    a1[j][i] = a2[i][j] = 0.5*d2[i][j];
  }

  form[j][i] = form[i][j];
  a12[j][i] = a12[i][j];
  sigma[j][i] = sigma[i][j];
  sigma3[j][i] = sigma3[i][j];
  sigma6[j][i] = sigma6[i][j];
  diameter[j][i] = diameter[i][j];

  double epsilon = a12[i][j]/144.0;
  lj1[j][i] = lj1[i][j] = 48.0 * epsilon * sigma6[i][j] * sigma6[i][j];
  lj2[j][i] = lj2[i][j] = 24.0 * epsilon * sigma6[i][j];
  lj3[j][i] = lj3[i][j] = 4.0 * epsilon * sigma6[i][j] * sigma6[i][j];
  lj4[j][i] = lj4[i][j] = 4.0 * epsilon * sigma6[i][j];

  offset[j][i] = offset[i][j] = 0.0;
  if (offset_flag && (cut[i][j] > 0.0)) {
    double tmp;
    offset[j][i] = offset[i][j] =
      single(0,0,i,j,cut[i][j]*cut[i][j],0.0,1.0,tmp);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairColloid::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a12[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&d1[i][j],sizeof(double),1,fp);
        fwrite(&d2[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairColloid::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (comm->me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (comm->me == 0) {
          utils::sfread(FLERR,&a12[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&d1[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&d2[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&a12[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&d1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&d2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairColloid::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairColloid::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairColloid::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,a12[i][i],sigma[i][i],d1[i][i],d2[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairColloid::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %g %g %g %g %g\n",i,
              a12[i][j],sigma[i][j],d1[i][j],d2[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairColloid::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                           double /*factor_coul*/, double factor_lj,
                           double &fforce)
{
  double K[9],h[4],g[4];
  double r,r2inv,r6inv,forcelj,c1,c2,phi,fR,dUR,dUA;

  phi = 0.0;
  switch (form[itype][jtype]) {
  case SMALL_SMALL:
    r2inv = 1.0/rsq;
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    fforce = factor_lj*forcelj*r2inv;
    phi = r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
      offset[itype][jtype];
    break;

  case SMALL_LARGE:
    c2 = a2[itype][jtype];
    K[1] = c2*c2;
    K[2] = rsq;
    K[0] = K[1] - rsq;
    K[4] = rsq*rsq;
    K[3] = K[1] - K[2];
    K[3] *= K[3]*K[3];
    K[6] = K[3]*K[3];
    fR = sigma3[itype][jtype]*a12[itype][jtype]*c2*K[1]/K[3];
    fforce = 4.0/15.0*fR*factor_lj *
      (2.0*(K[1]+K[2])*(K[1]*(5.0*K[1]+22.0*K[2])+5.0*K[4]) *
       sigma6[itype][jtype]/K[6] - 5.0)/K[0];
    phi = 2.0/9.0*fR *
      (1.0-(K[1]*(K[1]*(K[1]/3.0+3.0*K[2])+4.2*K[4])+K[2]*K[4]) *
       sigma6[itype][jtype]/K[6]) - offset[itype][jtype];
    break;

  case LARGE_LARGE:
    r = sqrt(rsq);
    c1 = a1[itype][jtype];
    c2 = a2[itype][jtype];
    K[0] = c1*c2;
    K[1] = c1+c2;
    K[2] = c1-c2;
    K[3] = K[1]+r;
    K[4] = K[1]-r;
    K[5] = K[2]+r;
    K[6] = K[2]-r;
    K[7] = 1.0/(K[3]*K[4]);
    K[8] = 1.0/(K[5]*K[6]);
    g[0] = powint(K[3],-7);
    g[1] = powint(K[4],-7);
    g[2] = powint(K[5],-7);
    g[3] = powint(K[6],-7);
    h[0] = ((K[3]+5.0*K[1])*K[3]+30.0*K[0])*g[0];
    h[1] = ((K[4]+5.0*K[1])*K[4]+30.0*K[0])*g[1];
    h[2] = ((K[5]+5.0*K[2])*K[5]-30.0*K[0])*g[2];
    h[3] = ((K[6]+5.0*K[2])*K[6]-30.0*K[0])*g[3];
    g[0] *= 42.0*K[0]/K[3]+6.0*K[1]+K[3];
    g[1] *= 42.0*K[0]/K[4]+6.0*K[1]+K[4];
    g[2] *= -42.0*K[0]/K[5]+6.0*K[2]+K[5];
    g[3] *= -42.0*K[0]/K[6]+6.0*K[2]+K[6];

    fR = a12[itype][jtype]*sigma6[itype][jtype]/r/37800.0;
    phi = fR * (h[0]-h[1]-h[2]+h[3]);
    dUR = phi/r + 5.0*fR*(g[0]+g[1]-g[2]-g[3]);
    dUA = -a12[itype][jtype]/3.0*r*((2.0*K[0]*K[7]+1.0)*K[7] +
                                    (2.0*K[0]*K[8]-1.0)*K[8]);
    fforce = factor_lj*(dUR+dUA)/r;
    phi += a12[itype][jtype]/6.0*(2.0*K[0]*(K[7]+K[8])-log(K[8]/K[7])) -
      offset[itype][jtype];
    break;
  }

  return factor_lj*phi;
}
