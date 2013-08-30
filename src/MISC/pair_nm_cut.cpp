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

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_nm_cut.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairNMCut::PairNMCut(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairNMCut::~PairNMCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
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

void PairNMCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,factor_lj;
  double r,forcenm,rminv,rninv;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);

        rminv = pow(r2inv,mm[itype][jtype]/2.0);
        rninv = pow(r2inv,nn[itype][jtype]/2.0);

        forcenm = e0nm[itype][jtype]*nm[itype][jtype] * 
          (r0n[itype][jtype]/pow(r,nn[itype][jtype]) - 
           r0m[itype][jtype]/pow(r,mm[itype][jtype]));
        fpair = factor_lj*forcenm*r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = e0nm[itype][jtype] * 
            (mm[itype][jtype]*r0n[itype][jtype]*rninv - 
             nn[itype][jtype]*r0m[itype][jtype]*rminv) - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairNMCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
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

void PairNMCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

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

void PairNMCut::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double e0_one = force->numeric(FLERR,arg[2]);
  double r0_one = force->numeric(FLERR,arg[3]);
  double nn_one = force->numeric(FLERR,arg[4]);
  double mm_one = force->numeric(FLERR,arg[5]);

  double cut_one = cut_global;
  if (narg == 7) cut_one = force->numeric(FLERR,arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      e0[i][j] = e0_one;
      r0[i][j] = r0_one;
      nn[i][j] = nn_one;
      mm[i][j] = mm_one;
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

double PairNMCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  nm[i][j] = nn[i][j]*mm[i][j];
  e0nm[i][j] = e0[i][j]/(nn[i][j]-mm[i][j]);
  r0n[i][j] = pow(r0[i][j],nn[i][j]);
  r0m[i][j] = pow(r0[i][j],mm[i][j]); 

  if (offset_flag) {
    offset[i][j] = e0nm[i][j] * 
      ((mm[i][j]*r0n[i][j] / pow(cut[i][j],nn[i][j])) - 
       (nn[i][j]*r0m[i][j] / pow(cut[i][j],mm[i][j])));
  } else offset[i][j] = 0.0;

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

    double rr1 = mm[i][j]*(nn[i][j]-1)*pow(r0[i][j],nn[i][j]);
    double rr2 = nn[i][j]*(mm[i][j]-1)*pow(r0[i][j],mm[i][j]);
    double p1 = 1-nn[i][j];
    double p2 = 1-mm[i][j];

    double rrr1 = pow(r0[i][j],nn[i][j])*(1-nn[i][j]);
    double rrr2 = pow(r0[i][j],mm[i][j])*(1-mm[i][j]);

    etail_ij = 2.0*MY_PI*all[0]*all[1]*e0nm[i][j] * 
      (rr1*pow(cut[i][j],p1)-rr2*pow(cut[i][j],p2));
    ptail_ij = 2.0*MY_PI*all[0]*all[1]*e0nm[i][j] * 
      nn[i][j]*mm[i][j]*(rrr1*pow(cut[i][j],p1)-rrr2*pow(cut[i][j],p2));
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNMCut::write_restart(FILE *fp)
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
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNMCut::read_restart(FILE *fp)
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
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&e0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&nn[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&mm[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNMCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNMCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
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

void PairNMCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,e0[i][i],r0[i][i],nn[i][i],mm[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairNMCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,
              e0[i][j],r0[i][j],nn[i][j],mm[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairNMCut::single(int i, int j, int itype, int jtype,
                      double rsq, double factor_coul, double factor_lj,
                      double &fforce)
{
  double r2inv,r,forcenm,phinm;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);

  forcenm = e0nm[itype][jtype]*nm[itype][jtype] * 
    (r0n[itype][jtype]/pow(r,nn[itype][jtype]) - 
     r0m[itype][jtype]/pow(r,mm[itype][jtype]));
  fforce = factor_lj*forcenm*r2inv;

  phinm = e0nm[itype][jtype] * 
    (mm[itype][jtype] * r0n[itype][jtype]/pow(r,nn[itype][jtype]) - 
     nn[itype][jtype]*r0m[itype][jtype] /pow(r,mm[itype][jtype])) - 
    offset[itype][jtype];
  return factor_lj*phinm;
}

/* ---------------------------------------------------------------------- */

void *PairNMCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"e0") == 0) return (void *) e0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"nn") == 0) return (void *) nn;
  if (strcmp(str,"mm") == 0) return (void *) mm;
  return NULL;
}
