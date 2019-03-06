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
   Contributing author: Hendrik Heenen (Technical University of Munich)
                        and Rochus Schmid (Ruhr-Universitaet Bochum)
   References: Bureekaew and Schmid, Phys. Status Solidi B 250, 1128 (2013)
               Fennell and Gezelter, JCP 124, 234104 (2006)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_buck6d_coul_gauss_dsf.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairBuck6dCoulGaussDSF::PairBuck6dCoulGaussDSF(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairBuck6dCoulGaussDSF::~PairBuck6dCoulGaussDSF()
{
  if (!copymode) {
    if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);

      memory->destroy(cut_lj);
      memory->destroy(cut_ljsq);
      memory->destroy(alpha_ij);
      memory->destroy(f_shift_ij);
      memory->destroy(e_shift_ij);
      memory->destroy(buck6d1);
      memory->destroy(buck6d2);
      memory->destroy(buck6d3);
      memory->destroy(buck6d4);
      memory->destroy(c0);
      memory->destroy(c1);
      memory->destroy(c2);
      memory->destroy(c3);
      memory->destroy(c4);
      memory->destroy(c5);
      memory->destroy(rsmooth_sq);
      memory->destroy(offset);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,r14inv,rexp,forcecoul,forcebuck6d,factor_coul,factor_lj;
  double term1,term2,term3,term4,term5;
  double rcu,rqu,sme,smf,ebuck6d;
  double prefactor,erfcc,erfcd,arg;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  double *special_coul = force->special_coul;
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
        r = sqrt(rsq);

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          r14inv = r6inv*r6inv*r2inv;
          rexp = exp(-r*buck6d2[itype][jtype]);
          term1 = buck6d3[itype][jtype]*r6inv;
          term2 = buck6d4[itype][jtype]*r14inv;
          term3 = term2*term2;
          term4 = 1.0/(1.0 + term2);
          term5 = 1.0/(1.0 + 2.0*term2 + term3);
          forcebuck6d = buck6d1[itype][jtype]*buck6d2[itype][jtype]*r*rexp;
          forcebuck6d -= term1*(6.0*term4 - term5*14.0*term2);
          ebuck6d = buck6d1[itype][jtype]*rexp - term1*term4;

          // smoothing term
          if (rsq > rsmooth_sq[itype][jtype]) {
            rcu = r*rsq;
            rqu = rsq*rsq;
            sme = c5[itype][jtype]*rqu*r + c4[itype][jtype]*rqu + c3[itype][jtype]*rcu +
                  c2[itype][jtype]*rsq + c1[itype][jtype]*r + c0[itype][jtype];
            smf = 5.0*c5[itype][jtype]*rqu + 4.0*c4[itype][jtype]*rcu +
                  3.0*c3[itype][jtype]*rsq + 2.0*c2[itype][jtype]*r + c1[itype][jtype];
            // forcebuck6d is -dE/dr*r
            forcebuck6d = forcebuck6d*sme - ebuck6d*smf*r; //RS was here: changed this from +E*smf to -E*smf*r
            ebuck6d *= sme;
          }
        } else forcebuck6d = 0.0;

        if (rsq < cut_coulsq) {
          prefactor = qqrd2e*qtmp*q[j]/r;

        //below: parametrization for erfcc = erf(alpha_ij[itype][jtype]*r);
          arg = alpha_ij[itype][jtype]*r;
          erfcd = MathSpecial::expmsq(arg);
          erfcc = 1 - (MathSpecial::my_erfcx(arg) * erfcd);

          forcecoul = prefactor * ((erfcc/r) - (2.0/MY_PIS*alpha_ij[itype][jtype]*erfcd) +
                                                r*f_shift_ij[itype][jtype]) * r;

          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        fpair = (forcecoul + factor_lj*forcebuck6d) * r2inv;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = ebuck6d - offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;

          if (rsq < cut_coulsq) {
            ecoul = prefactor * (erfcc - r*e_shift_ij[itype][jtype] -
                                 rsq*f_shift_ij[itype][jtype]);
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
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

void PairBuck6dCoulGaussDSF::allocate()
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
  memory->create(alpha_ij,n+1,n+1,"pair:alpha_ij");
  memory->create(f_shift_ij,n+1,n+1,"pair:f_shift_ij");
  memory->create(e_shift_ij,n+1,n+1,"pair:e_shift_ij");
  memory->create(buck6d1,n+1,n+1,"pair:buck6d1");
  memory->create(buck6d2,n+1,n+1,"pair:buck6d2");
  memory->create(buck6d3,n+1,n+1,"pair:buck6d3");
  memory->create(buck6d4,n+1,n+1,"pair:buck6d4");
  memory->create(c0,n+1,n+1,"pair:c0");
  memory->create(c1,n+1,n+1,"pair:c1");
  memory->create(c2,n+1,n+1,"pair:c2");
  memory->create(c3,n+1,n+1,"pair:c3");
  memory->create(c4,n+1,n+1,"pair:c4");
  memory->create(c5,n+1,n+1,"pair:c5");
  memory->create(rsmooth_sq,n+1,n+1,"pair:rsmooth_sq");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR,"Illegal pair_style command");

  vdwl_smooth = force->numeric(FLERR,arg[0]);

  cut_lj_global = force->numeric(FLERR,arg[1]);
  if (narg == 2) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[2]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double buck6d1_one = force->numeric(FLERR,arg[2]);
  double buck6d2_one = force->numeric(FLERR,arg[3]);
  double buck6d3_one = force->numeric(FLERR,arg[4]);
  double buck6d4_one = force->numeric(FLERR,arg[5]);
  double alpha_one = force->numeric(FLERR,arg[6]);

  double cut_lj_one = cut_lj_global;
  if (narg == 8) cut_lj_one = force->numeric(FLERR,arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      buck6d1[i][j] = buck6d1_one;
      buck6d2[i][j] = buck6d2_one;
      buck6d3[i][j] = buck6d3_one;
      buck6d4[i][j] = buck6d4_one;
      alpha_ij[i][j] = alpha_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style buck6d/coul/dsf requires atom attribute q");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBuck6dCoulGaussDSF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j],cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  //calculation of smoothing coefficients c0-c5
  c0[i][j] = c1[i][j] = c2[i][j] = c3[i][j] = c4[i][j] = c5[i][j] = 0.0;
  rsmooth_sq[i][j] = cut_ljsq[i][j];
  if (vdwl_smooth < 1.0) {
    double rsm = vdwl_smooth * cut_lj[i][j];
    double rsm_sq = rsm * rsm;
    double denom = pow((cut_lj[i][j]-rsm),5.0);
    c0[i][j] = cut_lj[i][j]*cut_ljsq[i][j]*(cut_ljsq[i][j]-
               5.0*cut_lj[i][j]*rsm+10.0*rsm_sq)/denom;
    c1[i][j] = -30.0*(cut_ljsq[i][j]*rsm_sq)/denom;
    c2[i][j] = 30.0*(cut_ljsq[i][j]*rsm + cut_lj[i][j]*rsm_sq)/denom;
    c3[i][j] = -10.0*(cut_ljsq[i][j] + 4.0*cut_lj[i][j]*rsm + rsm_sq)/denom;
    c4[i][j] = 15.0*(cut_lj[i][j]+rsm)/denom;
    c5[i][j] = -6.0/denom;
    rsmooth_sq[i][j] = rsm_sq;
  }

  // if offset_flag, shift is only invoked if there is not already smoothing
  if (offset_flag && vdwl_smooth >= 1.0) {
    double term1 = buck6d3[i][j]/pow(cut_lj[i][j],6.0);
    double term4 = 1.0/(1.0 + (buck6d4[i][j]/pow(cut_lj[i][j],14.0)));
    double rexp = exp(-cut_lj[i][j]*buck6d2[i][j]);
    offset[i][j] = buck6d1[i][j]*rexp - term1*term4;
  } else offset[i][j] = 0.0;

  double erfcd_c = exp(-alpha_ij[i][j]*alpha_ij[i][j]*cut_coul*cut_coul);
  double erfcc_c = erf(alpha_ij[i][j]*cut_coul);
  f_shift_ij[i][j] = -erfcc_c/cut_coulsq + 2.0/MY_PIS*alpha_ij[i][j]*erfcd_c/cut_coul;
  e_shift_ij[i][j] = erfcc_c/cut_coul - f_shift_ij[i][j]*cut_coul;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  alpha_ij[j][i] = alpha_ij[i][j];
  f_shift_ij[j][i] = f_shift_ij[i][j];
  e_shift_ij[j][i] = e_shift_ij[i][j];
  buck6d1[j][i] = buck6d1[i][j];
  buck6d2[j][i] = buck6d2[i][j];
  buck6d3[j][i] = buck6d3[i][j];
  buck6d4[j][i] = buck6d4[i][j];
  c0[j][i] = c0[i][j];
  c1[j][i] = c1[i][j];
  c2[j][i] = c2[i][j];
  c3[j][i] = c3[i][j];
  c4[j][i] = c4[i][j];
  c5[j][i] = c5[i][j];
  rsmooth_sq[j][i] = rsmooth_sq[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&buck6d1[i][j],sizeof(double),1,fp);
        fwrite(&buck6d2[i][j],sizeof(double),1,fp);
        fwrite(&buck6d3[i][j],sizeof(double),1,fp);
        fwrite(&buck6d4[i][j],sizeof(double),1,fp);
        fwrite(&alpha_ij[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
            }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::read_restart(FILE *fp)
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
          fread(&buck6d1[i][j],sizeof(double),1,fp);
          fread(&buck6d2[i][j],sizeof(double),1,fp);
          fread(&buck6d3[i][j],sizeof(double),1,fp);
          fread(&buck6d4[i][j],sizeof(double),1,fp);
          fread(&alpha_ij[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&buck6d1[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&buck6d2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&buck6d3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&buck6d4[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha_ij[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::write_restart_settings(FILE *fp)
{
  fwrite(&vdwl_smooth,sizeof(double),1,fp);
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&vdwl_smooth,sizeof(double),1,fp);
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&vdwl_smooth,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,
            buck6d1[i][i],buck6d2[i][i],buck6d3[i][i],
            buck6d4[i][i],alpha_ij[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussDSF::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g\n",i,j,
              buck6d1[i][j],buck6d2[i][j],buck6d3[i][j],
              buck6d4[i][j],alpha_ij[i][j],cut_lj[i][j]);
}
/* ---------------------------------------------------------------------- */

double PairBuck6dCoulGaussDSF::single(int i, int j, int itype, int jtype, double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r,r2inv,r6inv,r14inv,rexp,term1,term2,term3,term4,term5;
  double rcu,rqu,sme,smf;
  double erfcc,erfcd,prefactor,arg;
  double forcecoul,forcebuck6d,ebuck6d,phicoul,phibuck6d;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    r14inv = r6inv*r6inv*r2inv;
    rexp = exp(-r*buck6d2[itype][jtype]);
    term1 = buck6d3[itype][jtype]*r6inv;
    term2 = buck6d4[itype][jtype]*r14inv;
    term3 = term2*term2;
    term4 = 1.0/(1.0 + term2);
    term5 = 1.0/(1.0 + 2.0*term2 + term3);
    forcebuck6d = buck6d1[itype][jtype]*buck6d2[itype][jtype]*r*rexp;
    forcebuck6d -= term1*(6.0*term4 - term5*14.0*term2);
    ebuck6d = buck6d1[itype][jtype]*rexp - term1*term4;
    // smoothing term
    if (rsq > rsmooth_sq[itype][jtype]) {
      rcu = r*rsq;
      rqu = rsq*rsq;
      sme = c5[itype][jtype]*rqu*r + c4[itype][jtype]*rqu + c3[itype][jtype]*rcu +
            c2[itype][jtype]*rsq + c1[itype][jtype]*r + c0[itype][jtype];
      smf = 5.0*c5[itype][jtype]*rqu + 4.0*c4[itype][jtype]*rcu +
            3.0*c3[itype][jtype]*rsq + 2.0*c2[itype][jtype]*r + c1[itype][jtype];
      // forcebuck6d is -dE/dr*r
      forcebuck6d = forcebuck6d*sme - ebuck6d*smf*r; //RS was here: changed this from +E*smf to -E*smf*r
      ebuck6d *= sme;
    }
  } else forcebuck6d = 0.0;

  if (rsq < cut_coulsq) {
    prefactor = factor_coul * force->qqrd2e * atom->q[i] * atom->q[j] / r;
    arg = alpha_ij[itype][jtype]*r;
    erfcd = MathSpecial::expmsq(arg);
    erfcc = 1 - (MathSpecial::my_erfcx(arg) * erfcd);
    forcecoul = prefactor * ((erfcc/r) - (2.0/MY_PIS*alpha_ij[itype][jtype]*erfcd) +
                                          r*f_shift_ij[itype][jtype]) * r;
  } else forcecoul = 0.0;

  fforce = (forcecoul + factor_lj*forcebuck6d) * r2inv;

  double eng = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    phibuck6d = ebuck6d - offset[itype][jtype];
    eng += factor_lj*phibuck6d;
  }

  if (rsq < cut_coulsq) {
    phicoul = prefactor * (erfcc - r*e_shift_ij[itype][jtype] -
                                 rsq*f_shift_ij[itype][jtype]);
    eng += phicoul;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairBuck6dCoulGaussDSF::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_ljsq") == 0) return (void *) cut_ljsq;
  if (strcmp(str,"buck6d1") == 0) return (void *) buck6d1;
  if (strcmp(str,"buck6d2") == 0) return (void *) buck6d2;
  if (strcmp(str,"buck6d3") == 0) return (void *) buck6d3;
  if (strcmp(str,"buck6d4") == 0) return (void *) buck6d4;
  if (strcmp(str,"rsmooth_sq") == 0) return (void *) rsmooth_sq;
  if (strcmp(str,"c0") == 0) return (void *) c0;
  if (strcmp(str,"c1") == 0) return (void *) c1;
  if (strcmp(str,"c2") == 0) return (void *) c2;
  if (strcmp(str,"c3") == 0) return (void *) c3;
  if (strcmp(str,"c4") == 0) return (void *) c4;
  if (strcmp(str,"c5") == 0) return (void *) c5;
  if (strcmp(str,"offset") == 0) return (void *) offset;

  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  return NULL;
}
