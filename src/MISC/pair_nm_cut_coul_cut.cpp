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
   Contributing Author: Julien Devemy (ICCF)
------------------------------------------------------------------------- */

#include "pair_nm_cut_coul_cut.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairNMCutCoulCut::PairNMCutCoulCut(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairNMCutCoulCut::~PairNMCutCoulCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    memory->destroy(e0);
    memory->destroy(r0);
    memory->destroy(nn);
    memory->destroy(mm);
    memory->destroy(nm);
    memory->destroy(e0nm);
    memory->destroy(r0n);
    memory->destroy(r0m);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairNMCutCoulCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,factor_coul,factor_lj;
  double r,forcecoul,forcenm,rminv,rninv;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
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

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
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
        if (rsq < cut_coulsq[itype][jtype])
          forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r = sqrt(rsq);
          rminv = pow(r2inv,mm[itype][jtype]/2.0);
          rninv = pow(r2inv,nn[itype][jtype]/2.0);
          forcenm = e0nm[itype][jtype]*nm[itype][jtype] *
            (r0n[itype][jtype]/pow(r,nn[itype][jtype]) -
             r0m[itype][jtype]/pow(r,mm[itype][jtype]));
        } else forcenm = 0.0;

        fpair = (factor_coul*forcecoul + factor_lj*forcenm) * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype])
            ecoul = factor_coul * qqrd2e * qtmp*q[j]*sqrt(r2inv);
          else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = e0nm[itype][jtype]*(mm[itype][jtype] *
                                        r0n[itype][jtype]*rninv -
                                        nn[itype][jtype] *
                                        r0m[itype][jtype]*rminv) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairNMCutCoulCut::allocate()
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
  memory->create(cut_coul,n+1,n+1,"pair:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  memory->create(e0,n+1,n+1,"pair:e0");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(nn,n+1,n+1,"pair:nn");
  memory->create(mm,n+1,n+1,"pair:mm");
  memory->create(nm,n+1,n+1,"pair:nm");
  memory->create(e0nm,n+1,n+1,"pair:e0nm");
  memory->create(r0n,n+1,n+1,"pair:r0n");
  memory->create(r0m,n+1,n+1,"pair:r0m");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairNMCutCoulCut::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = force->numeric(FLERR,arg[0]);
  if (narg == 1) cut_coul_global = cut_lj_global;
  else cut_coul_global = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_lj[i][j] = cut_lj_global;
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairNMCutCoulCut::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double e0_one = force->numeric(FLERR,arg[2]);
  double r0_one = force->numeric(FLERR,arg[3]);
  double nn_one = force->numeric(FLERR,arg[4]);
  double mm_one = force->numeric(FLERR,arg[5]);

  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;
  if (narg >= 7) cut_coul_one = cut_lj_one = force->numeric(FLERR,arg[4]);
  if (narg == 8) cut_coul_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      e0[i][j] = e0_one;
      r0[i][j] = r0_one;
      nn[i][j] = nn_one;
      mm[i][j] = mm_one;
      cut_lj[i][j] = cut_lj_one;
      cut_coul[i][j] = cut_coul_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairNMCutCoulCut::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style nm/cut/coul/cut requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNMCutCoulCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j],cut_coul[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

  nm[i][j] = nn[i][j]*mm[i][j];
  e0nm[i][j] = e0[i][j]/(nn[i][j]-mm[i][j]);
  r0n[i][j] = pow(r0[i][j],nn[i][j]);
  r0m[i][j] = pow(r0[i][j],mm[i][j]);

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    offset[i][j] = e0nm[i][j] *
      ((mm[i][j]*r0n[i][j] / pow(cut_lj[i][j],nn[i][j])) -
       (nn[i][j]*r0m[i][j] / pow(cut_lj[i][j],mm[i][j])));
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  e0[j][i] = e0[i][j];
  nn[j][i] = nn[i][j];
  mm[j][i] = mm[i][j];
  nm[j][i] = nm[i][j];
  r0[j][i] = r0[i][j];
  e0nm[j][i] = e0nm[i][j];
  r0n[j][i] = r0n[i][j];
  r0m[j][i] = r0m[i][j];
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

    double cut_lj3 = cut_lj[i][j]*cut_lj[i][j]*cut_lj[i][j];
    ptail_ij = 2.*MY_PI/3.*all[0]*all[1]*e0nm[i][j]*nm[i][j]*cut_lj3 *
      (pow(r0[i][j]/cut_lj[i][j],nn[i][j])/(nn[i][j]-3) - pow(r0[i][j]/cut_lj[i][j],mm[i][j])/(mm[i][j]-3));
    etail_ij = 2.*MY_PI*all[0]*all[1]*e0nm[i][j]*cut_lj3 *
      (mm[i][j]*pow(r0[i][j]/cut_lj[i][j],nn[i][j])/(nn[i][j]-3) - nn[i][j]*pow(r0[i][j]/cut_lj[i][j],mm[i][j])/(mm[i][j]-3));

  }

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNMCutCoulCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&e0[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&nn[i][j],sizeof(double),1,fp);
        fwrite(&mm[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
        fwrite(&cut_coul[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNMCutCoulCut::read_restart(FILE *fp)
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
          fread(&e0[i][j],sizeof(double),1,fp);
          fread(&r0[i][j],sizeof(double),1,fp);
          fread(&nn[i][j],sizeof(double),1,fp);
          fread(&mm[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
          fread(&cut_coul[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&e0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&nn[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&mm[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNMCutCoulCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNMCutCoulCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairNMCutCoulCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,e0[i][i],r0[i][i],nn[i][i],mm[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairNMCutCoulCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,
              e0[i][j],r0[i][j],nn[i][j],mm[i][j],cut_lj[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairNMCutCoulCut::single(int i, int j, int itype, int jtype,
                                double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r,forcecoul,forcenm,phicoul,phinm;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq[itype][jtype])
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j]*sqrt(r2inv);
  else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r = sqrt(rsq);
    forcenm = e0nm[itype][jtype]*nm[itype][jtype] *
      (r0n[itype][jtype]/pow(r,nn[itype][jtype]) -
       r0m[itype][jtype]/pow(r,mm[itype][jtype]));
  } else forcenm = 0.0;
  fforce = (factor_coul*forcecoul + factor_lj*forcenm) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j]*sqrt(r2inv);
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    phinm = e0nm[itype][jtype] *
      (mm[itype][jtype]*r0n[itype][jtype]/pow(r,nn[itype][jtype]) -
       nn[itype][jtype]*r0m[itype][jtype]/pow(r,mm[itype][jtype])) -
      offset[itype][jtype];
    eng += factor_lj*phinm;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairNMCutCoulCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  if (strcmp(str,"e0") == 0) return (void *) e0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"nn") == 0) return (void *) nn;
  if (strcmp(str,"mm") == 0) return (void *) mm;
  return NULL;
}
