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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_thole.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "fix.h"
#include "fix_store.h"
#include "domain.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairThole::PairThole(LAMMPS *lmp) : Pair(lmp) {
    fix_drude = NULL;
}

/* ---------------------------------------------------------------------- */

PairThole::~PairThole()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(polar);
    memory->destroy(thole);
    memory->destroy(ascreen);
    memory->destroy(cut);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairThole::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qi,qj,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,rsq,r2inv,rinv,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_f,factor_e;
  int di,dj;
  double dcoul,asr,exp_asr;

  ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  int *drudetype = fix_drude->drudetype;
  tagint *drudeid = fix_drude->drudeid;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    // only on core-drude pair
    if (drudetype[type[i]] == NOPOL_TYPE)
      continue;

    di = domain->closest_image(i, atom->map(drudeid[i]));
    // get dq of the core via the drude charge
    if (drudetype[type[i]] == DRUDE_TYPE)
      qi = q[i];
    else
      qi = -q[di];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // only on core-drude pair, but not into the same pair
      if (drudetype[type[j]] == NOPOL_TYPE || j == di)
        continue;

      // get dq of the core via the drude charge
      if (drudetype[type[j]] == DRUDE_TYPE)
        qj = q[j];
      else {
        dj = domain->closest_image(j, atom->map(drudeid[j]));
        qj = -q[dj];
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        r = sqrt(rsq);
        asr = ascreen[itype][jtype] * r;
        exp_asr = exp(-asr);
        dcoul = qqrd2e * qi * qj *scale[itype][jtype] * rinv;
        factor_f = 0.5*(2. + (exp_asr * (-2. - asr * (2. + asr)))) - factor_coul;
        if(eflag) factor_e = 0.5*(2. - (exp_asr * (2. + asr))) - factor_coul;
        fpair = factor_f * dcoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag)
          ecoul = factor_e * dcoul;

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairThole::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(scale,n+1,n+1,"pair:scale");
  memory->create(ascreen,n+1,n+1,"pair:ascreen");
  memory->create(thole,n+1,n+1,"pair:thole");
  memory->create(polar,n+1,n+1,"pair:polar");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairThole::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  thole_global = force->numeric(FLERR,arg[0]);
  cut_global = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          thole[i][j] = thole_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairThole::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double polar_one = force->numeric(FLERR,arg[2]);
  double thole_one = thole_global;
  double cut_one = cut_global;
  if (narg >=4) thole_one = force->numeric(FLERR,arg[3]);
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      polar[i][j] = polar_one;
      thole[i][j] = thole_one;
      ascreen[i][j] = thole[i][j] / pow(polar[i][j], 1./3.);
      cut[i][j] = cut_one;
      scale[i][j] = 1.0;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairThole::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style thole requires atom attribute q");
  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(modify->fix[ifix]->style,"drude") == 0) break;
  if (ifix == modify->nfix) error->all(FLERR, "Pair thole requires fix drude");
  fix_drude = (FixDrude *) modify->fix[ifix];

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairThole::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);

  polar[j][i] = polar[i][j];
  thole[j][i] = thole[i][j];
  ascreen[j][i] = ascreen[i][j];
  scale[j][i] = scale[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairThole::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&polar[i][j],sizeof(double),1,fp);
        fwrite(&thole[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairThole::read_restart(FILE *fp)
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
          fread(&polar[i][j],sizeof(double),1,fp);
          fread(&thole[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
          ascreen[i][j] = thole[i][j] / pow(polar[i][j], 1./3.);
          }
        MPI_Bcast(&polar[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&thole[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&ascreen[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairThole::write_restart_settings(FILE *fp)
{
  fwrite(&thole_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairThole::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&thole_global,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&thole_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairThole::single(int i, int j, int itype, int jtype,
                         double rsq, double factor_coul, double /*factor_lj*/,
                         double &fforce)
{
  double r2inv,rinv,r,phicoul;
  double qi,qj,factor_f,factor_e,dcoul,asr,exp_asr;
  int di, dj;

  int *drudetype = fix_drude->drudetype;
  tagint *drudeid = fix_drude->drudeid;
  int *type = atom->type;

  // only on core-drude pair, but not on the same pair
  if (drudetype[type[i]] == NOPOL_TYPE || drudetype[type[j]] == NOPOL_TYPE ||
      j == i)
    return 0.0;

  // get dq of the core via the drude charge
  if (drudetype[type[i]] == DRUDE_TYPE)
    qi = atom->q[i];
  else {
    di = domain->closest_image(i, atom->map(drudeid[i]));
    qi = -atom->q[di];
  }
  if (drudetype[type[j]] == DRUDE_TYPE)
    qj = atom->q[j];
  else {
    dj = domain->closest_image(j, atom->map(drudeid[j]));
    qj = -atom->q[dj];
  }

  r2inv = 1.0/rsq;
  fforce = phicoul = 0.0;
  if (rsq < cutsq[itype][jtype]) {
    rinv = sqrt(r2inv);
    r = sqrt(rsq);
    asr = ascreen[itype][jtype] * r;
    exp_asr = exp(-asr);
    dcoul = force->qqrd2e * qi * qj * scale[itype][jtype] * rinv;
    factor_f = 0.5*(2. + (exp_asr * (-2. - asr * (2. + asr)))) - factor_coul;
    fforce = factor_f * dcoul * r2inv;
    factor_e = 0.5*(2. - (exp_asr * (2. + asr))) - factor_coul;
    phicoul = factor_e * dcoul;
  }

  return phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairThole::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  if (strcmp(str,"polar") == 0) return (void *) polar;
  if (strcmp(str,"thole") == 0) return (void *) thole;
  if (strcmp(str,"ascreen") == 0) return (void *) ascreen;
  return NULL;
}
