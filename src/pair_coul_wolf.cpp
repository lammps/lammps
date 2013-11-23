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
   Contributing author: Yongfeng Zhang (INL), yongfeng.zhang@inl.gov
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_coul_wolf.h"
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

PairCoulWolf::PairCoulWolf(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairCoulWolf::~PairCoulWolf()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulWolf::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,forcecoul,factor_coul;
  double prefactor;
  double r,rexp;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double erfcc,erfcd,v_sh,dvdrr,e_self,e_shift,f_shift,qisq;

  ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  // self and shifted coulombic energy

  e_self = v_sh = 0.0;
  e_shift = erfc(alf*cut_coul)/cut_coul;
  f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /
    cut_coul;

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

    qisq = qtmp*qtmp;
    e_self = -(e_shift/2.0 + alf/MY_PIS) * qisq*qqrd2e;
    if (evflag) ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        r = sqrt(rsq);
        prefactor = qqrd2e*qtmp*q[j]/r;
        erfcc = erfc(alf*r);
        erfcd = exp(-alf*alf*r*r);
        v_sh = (erfcc - e_shift*r) * prefactor;
        dvdrr = (erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r) + f_shift;
        forcecoul = dvdrr*rsq*prefactor;
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        fpair = forcecoul / rsq;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = v_sh;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
        }

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

void PairCoulWolf::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
   unlike other pair styles,
     there are no individual pair settings that these override
------------------------------------------------------------------------- */

void PairCoulWolf::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  alf = force->numeric(FLERR,arg[0]);
  cut_coul = force->numeric(FLERR,arg[1]);
}

/* ----------------------------------------------------------------------
   set cutoffs for one or more type pairs, optional
------------------------------------------------------------------------- */

void PairCoulWolf::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulWolf::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair coul/wolf requires atom attribute q");

  int irequest = neighbor->request(this);

  cut_coulsq = cut_coul*cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulWolf::init_one(int i, int j)
{
  return cut_coul;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulWolf::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fwrite(&setflag[i][j],sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulWolf::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulWolf::write_restart_settings(FILE *fp)
{
  fwrite(&alf,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulWolf::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&alf,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&alf,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   only the pair part is calculated here
------------------------------------------------------------------------- */

double PairCoulWolf::single(int i, int j, int itype, int jtype, double rsq,
                            double factor_coul, double factor_lj,
                            double &fforce)
{
  double r6inv,r,prefactor,rexp;
  double forcecoul,forceborn,phicoul;
  double e_shift,f_shift,dvdrr,erfcc,erfcd;

  e_shift = erfc(alf*cut_coul) / cut_coul;
  f_shift = -(e_shift+ 2.0*alf/MY_PIS * exp(-alf*alf*cut_coul*cut_coul)) /
    cut_coul;

  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
    erfcc = erfc(alf*r);
    erfcd = exp(-alf*alf*r*r);
    dvdrr = (erfcc/rsq + 2.0*alf/MY_PIS * erfcd/r) + f_shift;
    forcecoul = dvdrr*rsq*prefactor;
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  } else forcecoul = 0.0;
  fforce = forcecoul / rsq;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    phicoul = prefactor * (erfcc-e_shift*r);
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    eng += phicoul;
  }
  return eng;
}
