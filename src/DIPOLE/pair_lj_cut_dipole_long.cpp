/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_lj_cut_dipole_long.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "update.h"


using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutDipoleLong::PairLJCutDipoleLong(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  ewaldflag = dipoleflag = 1;
  respa_enable = 0;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJCutDipoleLong::~PairLJCutDipoleLong()
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

void PairLJCutDipoleLong::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,r,rinv,r2inv,r6inv;
  double forcecoulx,forcecouly,forcecoulz,fforce;
  double tixcoul,tiycoul,tizcoul,tjxcoul,tjycoul,tjzcoul;
  double fx,fy,fz,fdx,fdy,fdz,fax,fay,faz;
  double pdotp,pidotr,pjdotr,pre1,pre2,pre3;
  double grij,expm2,t,erfc;
  double g0,g1,g2,b0,b1,b2,b3,d0,d1,d2,d3;
  double zdix,zdiy,zdiz,zdjx,zdjy,zdjz,zaix,zaiy,zaiz,zajx,zajy,zajz;
  double g0b1_g1b2_g2b3,g0d1_g1d2_g2d3;
  double forcelj,factor_coul,factor_lj,facm1;
  double evdwl,ecoul;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  pre1 = 2.0 * g_ewald / MY_PIS;
  pre2 = 4.0 * pow(g_ewald,3.0) / MY_PIS;
  pre3 = 8.0 * pow(g_ewald,5.0) / MY_PIS;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = atom->q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          pdotp = mu[i][0]*mu[j][0] + mu[i][1]*mu[j][1] + mu[i][2]*mu[j][2];
          pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
          pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;

          g0 = qtmp*q[j];
          g1 = qtmp*pjdotr - q[j]*pidotr + pdotp;
          g2 = -pidotr*pjdotr;

          if (factor_coul > 0.0) {
            b0 = erfc * rinv;
            b1 = (b0 + pre1*expm2) * r2inv;
            b2 = (3.0*b1 + pre2*expm2) * r2inv;
            b3 = (5.0*b2 + pre3*expm2) * r2inv;

            g0b1_g1b2_g2b3 = g0*b1 + g1*b2 + g2*b3;
            fdx = delx * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][0] - q[j]*mu[i][0]) +
              b2 * (pjdotr*mu[i][0] + pidotr*mu[j][0]);
            fdy = dely * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][1] - q[j]*mu[i][1]) +
              b2 * (pjdotr*mu[i][1] + pidotr*mu[j][1]);
            fdz = delz * g0b1_g1b2_g2b3 -
              b1 * (qtmp*mu[j][2] - q[j]*mu[i][2]) +
              b2 * (pjdotr*mu[i][2] + pidotr*mu[j][2]);

            zdix = delx * (q[j]*b1 + b2*pjdotr) - b1*mu[j][0];
            zdiy = dely * (q[j]*b1 + b2*pjdotr) - b1*mu[j][1];
            zdiz = delz * (q[j]*b1 + b2*pjdotr) - b1*mu[j][2];
            zdjx = delx * (-qtmp*b1 + b2*pidotr) - b1*mu[i][0];
            zdjy = dely * (-qtmp*b1 + b2*pidotr) - b1*mu[i][1];
            zdjz = delz * (-qtmp*b1 + b2*pidotr) - b1*mu[i][2];

            if (factor_coul < 1.0) {
              fdx *= factor_coul;
              fdy *= factor_coul;
              fdz *= factor_coul;
              zdix *= factor_coul;
              zdiy *= factor_coul;
              zdiz *= factor_coul;
              zdjx *= factor_coul;
              zdjy *= factor_coul;
              zdjz *= factor_coul;
            }
          } else {
            fdx = fdy = fdz = 0.0;
            zdix = zdiy = zdiz = 0.0;
            zdjx = zdjy = zdjz = 0.0;
          }

          if (factor_coul < 1.0) {
            d0 = (erfc - 1.0) * rinv;
            d1 = (d0 + pre1*expm2) * r2inv;
            d2 = (3.0*d1 + pre2*expm2) * r2inv;
            d3 = (5.0*d2 + pre3*expm2) * r2inv;

            g0d1_g1d2_g2d3 = g0*d1 + g1*d2 + g2*d3;
            fax = delx * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][0] - q[j]*mu[i][0]) +
              d2 * (pjdotr*mu[i][0] + pidotr*mu[j][0]);
            fay = dely * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][1] - q[j]*mu[i][1]) +
              d2 * (pjdotr*mu[i][1] + pidotr*mu[j][1]);
            faz = delz * g0d1_g1d2_g2d3 -
              d1 * (qtmp*mu[j][2] - q[j]*mu[i][2]) +
              d2 * (pjdotr*mu[i][2] + pidotr*mu[j][2]);

            zaix = delx * (q[j]*d1 + d2*pjdotr) - d1*mu[j][0];
            zaiy = dely * (q[j]*d1 + d2*pjdotr) - d1*mu[j][1];
            zaiz = delz * (q[j]*d1 + d2*pjdotr) - d1*mu[j][2];
            zajx = delx * (-qtmp*d1 + d2*pidotr) - d1*mu[i][0];
            zajy = dely * (-qtmp*d1 + d2*pidotr) - d1*mu[i][1];
            zajz = delz * (-qtmp*d1 + d2*pidotr) - d1*mu[i][2];

            if (factor_coul > 0.0) {
              facm1 = 1.0 - factor_coul;
              fax *= facm1;
              fay *= facm1;
              faz *= facm1;
              zaix *= facm1;
              zaiy *= facm1;
              zaiz *= facm1;
              zajx *= facm1;
              zajy *= facm1;
              zajz *= facm1;
            }
          } else {
            fax = fay = faz = 0.0;
            zaix = zaiy = zaiz = 0.0;
            zajx = zajy = zajz = 0.0;
          }

          forcecoulx = fdx + fax;
          forcecouly = fdy + fay;
          forcecoulz = fdz + faz;

          tixcoul = mu[i][1]*(zdiz + zaiz) - mu[i][2]*(zdiy + zaiy);
          tiycoul = mu[i][2]*(zdix + zaix) - mu[i][0]*(zdiz + zaiz);
          tizcoul = mu[i][0]*(zdiy + zaiy) - mu[i][1]*(zdix + zaix);
          tjxcoul = mu[j][1]*(zdjz + zajz) - mu[j][2]*(zdjy + zajy);
          tjycoul = mu[j][2]*(zdjx + zajx) - mu[j][0]*(zdjz + zajz);
          tjzcoul = mu[j][0]*(zdjy + zajy) - mu[j][1]*(zdjx + zajx);

        } else {
          forcecoulx = forcecouly = forcecoulz = 0.0;
          tixcoul = tiycoul = tizcoul = 0.0;
          tjxcoul = tjycoul = tjzcoul = 0.0;
        }

        // LJ interaction

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          fforce = factor_lj * forcelj*r2inv;
        } else fforce = 0.0;

        // total force

        fx = qqrd2e*forcecoulx + delx*fforce;
        fy = qqrd2e*forcecouly + dely*fforce;
        fz = qqrd2e*forcecoulz + delz*fforce;

        // force & torque accumulation

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        torque[i][0] += qqrd2e*tixcoul;
        torque[i][1] += qqrd2e*tiycoul;
        torque[i][2] += qqrd2e*tizcoul;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] += qqrd2e*tjxcoul;
          torque[j][1] += qqrd2e*tjycoul;
          torque[j][2] += qqrd2e*tjzcoul;
        }

        if (eflag) {
          if (rsq < cut_coulsq && factor_coul > 0.0) {
            ecoul = qqrd2e*(b0*g0 + b1*g1 + b2*g2);
            if (factor_coul < 1.0) {
              ecoul *= factor_coul;
              ecoul += (1-factor_coul) * qqrd2e * (d0*g0 + d1*g1 + d2*g2);
            }
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,ecoul,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutDipoleLong::allocate()
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

void PairLJCutDipoleLong::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect args in pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);
  if (narg == 1) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutDipoleLong::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

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
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutDipoleLong::init_one(int i, int j)
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

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutDipoleLong::init_style()
{
  if (!atom->q_flag || !atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair dipole/long requires atom attributes q, mu, torque");

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");

  g_ewald = force->kspace->g_ewald;

  cut_coulsq = cut_coul * cut_coul;

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutDipoleLong::write_restart(FILE *fp)
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

void PairLJCutDipoleLong::read_restart(FILE *fp)
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

void PairLJCutDipoleLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutDipoleLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void *PairLJCutDipoleLong::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  } else if (strcmp(str,"ewald_order") == 0) {
    ewald_order = 0;
    ewald_order |= 1<<1;
    ewald_order |= 1<<3;
    dim = 0;
    return (void *) &ewald_order;
  } else if (strcmp(str,"ewald_mix") == 0) {
    dim = 0;
    return (void *) &mix_flag;
  }
  return NULL;
}
