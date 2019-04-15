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
#include <cstdlib>
#include <cstring>
#include "pair_lj_cut_dipole_cut.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutDipoleCut::PairLJCutDipoleCut(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairLJCutDipoleCut::~PairLJCutDipoleCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
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

void PairLJCutDipoleCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fx,fy,fz;
  double rsq,rinv,r2inv,r6inv,r3inv,r5inv,r7inv;
  double forcecoulx,forcecouly,forcecoulz,crossx,crossy,crossz;
  double tixcoul,tiycoul,tizcoul,tjxcoul,tjycoul,tjzcoul;
  double fq,pdotp,pidotr,pjdotr,pre1,pre2,pre3,pre4;
  double forcelj,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

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
        rinv = sqrt(r2inv);

        // atom can have both a charge and dipole
        // i,j = charge-charge, dipole-dipole, dipole-charge, or charge-dipole

        forcecoulx = forcecouly = forcecoulz = 0.0;
        tixcoul = tiycoul = tizcoul = 0.0;
        tjxcoul = tjycoul = tjzcoul = 0.0;

        if (rsq < cut_coulsq[itype][jtype]) {

          if (qtmp != 0.0 && q[j] != 0.0) {
            r3inv = r2inv*rinv;
            pre1 = qtmp*q[j]*r3inv;

            forcecoulx += pre1*delx;
            forcecouly += pre1*dely;
            forcecoulz += pre1*delz;
          }

          if (mu[i][3] > 0.0 && mu[j][3] > 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            r7inv = r5inv*r2inv;

            pdotp = mu[i][0]*mu[j][0] + mu[i][1]*mu[j][1] + mu[i][2]*mu[j][2];
            pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
            pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;

            pre1 = 3.0*r5inv*pdotp - 15.0*r7inv*pidotr*pjdotr;
            pre2 = 3.0*r5inv*pjdotr;
            pre3 = 3.0*r5inv*pidotr;
            pre4 = -1.0*r3inv;

            forcecoulx += pre1*delx + pre2*mu[i][0] + pre3*mu[j][0];
            forcecouly += pre1*dely + pre2*mu[i][1] + pre3*mu[j][1];
            forcecoulz += pre1*delz + pre2*mu[i][2] + pre3*mu[j][2];

            crossx = pre4 * (mu[i][1]*mu[j][2] - mu[i][2]*mu[j][1]);
            crossy = pre4 * (mu[i][2]*mu[j][0] - mu[i][0]*mu[j][2]);
            crossz = pre4 * (mu[i][0]*mu[j][1] - mu[i][1]*mu[j][0]);

            tixcoul += crossx + pre2 * (mu[i][1]*delz - mu[i][2]*dely);
            tiycoul += crossy + pre2 * (mu[i][2]*delx - mu[i][0]*delz);
            tizcoul += crossz + pre2 * (mu[i][0]*dely - mu[i][1]*delx);
            tjxcoul += -crossx + pre3 * (mu[j][1]*delz - mu[j][2]*dely);
            tjycoul += -crossy + pre3 * (mu[j][2]*delx - mu[j][0]*delz);
            tjzcoul += -crossz + pre3 * (mu[j][0]*dely - mu[j][1]*delx);
          }

          if (mu[i][3] > 0.0 && q[j] != 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
            pre1 = 3.0*q[j]*r5inv * pidotr;
            pre2 = q[j]*r3inv;

            forcecoulx += pre2*mu[i][0] - pre1*delx;
            forcecouly += pre2*mu[i][1] - pre1*dely;
            forcecoulz += pre2*mu[i][2] - pre1*delz;
            tixcoul += pre2 * (mu[i][1]*delz - mu[i][2]*dely);
            tiycoul += pre2 * (mu[i][2]*delx - mu[i][0]*delz);
            tizcoul += pre2 * (mu[i][0]*dely - mu[i][1]*delx);
          }

          if (mu[j][3] > 0.0 && qtmp != 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;
            pre1 = 3.0*qtmp*r5inv * pjdotr;
            pre2 = qtmp*r3inv;

            forcecoulx += pre1*delx - pre2*mu[j][0];
            forcecouly += pre1*dely - pre2*mu[j][1];
            forcecoulz += pre1*delz - pre2*mu[j][2];
            tjxcoul += -pre2 * (mu[j][1]*delz - mu[j][2]*dely);
            tjycoul += -pre2 * (mu[j][2]*delx - mu[j][0]*delz);
            tjzcoul += -pre2 * (mu[j][0]*dely - mu[j][1]*delx);
          }
        }

        // LJ interaction

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          forcelj *= factor_lj * r2inv;
        } else forcelj = 0.0;

        // total force

        fq = factor_coul*qqrd2e;
        fx = fq*forcecoulx + delx*forcelj;
        fy = fq*forcecouly + dely*forcelj;
        fz = fq*forcecoulz + delz*forcelj;

        // force & torque accumulation

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        torque[i][0] += fq*tixcoul;
        torque[i][1] += fq*tiycoul;
        torque[i][2] += fq*tizcoul;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] += fq*tjxcoul;
          torque[j][1] += fq*tjycoul;
          torque[j][2] += fq*tjzcoul;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype]) {
            ecoul = qtmp*q[j]*rinv;
            if (mu[i][3] > 0.0 && mu[j][3] > 0.0)
              ecoul += r3inv*pdotp - 3.0*r5inv*pidotr*pjdotr;
            if (mu[i][3] > 0.0 && q[j] != 0.0)
              ecoul += -q[j]*r3inv*pidotr;
            if (mu[j][3] > 0.0 && qtmp != 0.0)
              ecoul += qtmp*r3inv*pjdotr;
            ecoul *= factor_coul*qqrd2e;
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

void PairLJCutDipoleCut::allocate()
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

void PairLJCutDipoleCut::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect args in pair_style command");

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

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

void PairLJCutDipoleCut::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;
  if (narg >= 5) cut_coul_one = cut_lj_one = force->numeric(FLERR,arg[4]);
  if (narg == 6) cut_coul_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
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

void PairLJCutDipoleCut::init_style()
{
  if (!atom->q_flag || !atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair dipole/cut requires atom attributes q, mu, torque");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutDipoleCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i],cut_lj[j][j]);
    cut_coul[i][j] = mix_distance(cut_coul[i][i],cut_coul[j][j]);
  }

  double cut = MAX(cut_lj[i][j],cut_coul[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutDipoleCut::write_restart(FILE *fp)
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
        fwrite(&cut_coul[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutDipoleCut::read_restart(FILE *fp)
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
          fread(&cut_coul[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutDipoleCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutDipoleCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
