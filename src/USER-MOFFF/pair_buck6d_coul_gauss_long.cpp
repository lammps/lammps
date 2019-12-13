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
------------------------------------------------------------------------- */

#include "pair_buck6d_coul_gauss_long.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "math_special.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917

/* ---------------------------------------------------------------------- */

PairBuck6dCoulGaussLong::PairBuck6dCoulGaussLong(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  single_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairBuck6dCoulGaussLong::~PairBuck6dCoulGaussLong()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(alpha_ij);
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

/* ---------------------------------------------------------------------- */

void PairBuck6dCoulGaussLong::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,r14inv,rexp,forcecoul,forcebuck6d,factor_coul,factor_lj;
  double grij,expm2,erf;
  double term1,term2,term3,term4,term5;
  double rcu,rqu,sme,smf,ebuck6d,ealpha;
  double prefactor,erfa,expa,arg,falpha;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

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
            forcebuck6d = forcebuck6d*sme - ebuck6d*smf*r;
            ebuck6d *= sme;
          }
        } else forcebuck6d = 0.0;

        if (rsq < cut_coulsq) {
        // long range - real space
          grij = g_ewald * r;
          expm2 = MathSpecial::expmsq(grij);
          erf = 1 - (MathSpecial::my_erfcx(grij) * expm2);

        // gaussian for 1/r alpha_ij contribution
          arg = alpha_ij[itype][jtype]*r;
          expa = MathSpecial::expmsq(arg);
          erfa = 1 - (MathSpecial::my_erfcx(arg) * expa);

          prefactor = qqrd2e*qtmp*q[j]/r;
          falpha = erfa - EWALD_F*arg*expa;
          forcecoul = prefactor * (falpha - erf + EWALD_F*grij*expm2);
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor*falpha;

          // (q*q/r) * (gauss(alpha_ij) - gauss(alpha_long)
          ealpha = prefactor * (erfa-erf);
          // smoothing term - NOTE: ingnored in special_bonds correction
          // since likely rsmooth_sq_c >> d(special)
          if (rsq > rsmooth_sq_c) {
            rcu = r*rsq;
            rqu = rsq*rsq;
            sme = c5_c*rqu*r + c4_c*rqu + c3_c*rcu + c2_c*rsq + c1_c*r + c0_c;
            smf = 5.0*c5_c*rqu + 4.0*c4_c*rcu + 3.0*c3_c*rsq + 2.0*c2_c*r + c1_c;
            forcecoul = forcecoul*sme - ealpha*smf*r;
            ealpha *= sme;
          }
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
            ecoul = ealpha;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor*erfa;
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

void PairBuck6dCoulGaussLong::allocate()
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

void PairBuck6dCoulGaussLong::settings(int narg, char **arg)
{
  if (narg < 3 || narg > 4) error->all(FLERR,"Illegal pair_style command");

  vdwl_smooth = force->numeric(FLERR,arg[0]);
  coul_smooth = force->numeric(FLERR,arg[1]);

  cut_lj_global = force->numeric(FLERR,arg[2]);
  if (narg == 3) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[3]);

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

void PairBuck6dCoulGaussLong::coeff(int narg, char **arg)
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

void PairBuck6dCoulGaussLong::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style buck6d/coul/dsf requires atom attribute q");

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;

  //calculation of smoothing coefficients c0_c-c5_c for coulomb smoothing
  c0_c = c1_c = c2_c = c3_c = c4_c = c5_c = 0.0;
  rsmooth_sq_c = cut_coulsq;
  if (coul_smooth < 1.0) {
    double rsm = coul_smooth * cut_coul;
    double rsm_sq = rsm * rsm;
    double denom = pow((cut_coul-rsm),5.0);
    c0_c = cut_coul*cut_coulsq*(cut_coulsq-
           5.0*cut_coul*rsm+10.0*rsm_sq)/denom;
    c1_c = -30.0*(cut_coulsq*rsm_sq)/denom;
    c2_c = 30.0*(cut_coulsq*rsm + cut_coul*rsm_sq)/denom;
    c3_c = -10.0*(cut_coulsq + 4.0*cut_coul*rsm + rsm_sq)/denom;
    c4_c = 15.0*(cut_coul+rsm)/denom;
    c5_c = -6.0/denom;
    rsmooth_sq_c = rsm_sq;
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBuck6dCoulGaussLong::init_one(int i, int j)
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

  cut_ljsq[j][i] = cut_ljsq[i][j];
  alpha_ij[j][i] = alpha_ij[i][j];
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

void PairBuck6dCoulGaussLong::write_restart(FILE *fp)
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

void PairBuck6dCoulGaussLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&buck6d1[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&buck6d2[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&buck6d3[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&buck6d4[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&alpha_ij[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut_lj[i][j],sizeof(double),1,fp,NULL,error);
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

void PairBuck6dCoulGaussLong::write_restart_settings(FILE *fp)
{
  fwrite(&vdwl_smooth,sizeof(double),1,fp);
  fwrite(&coul_smooth,sizeof(double),1,fp);
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&vdwl_smooth,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&coul_smooth,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_lj_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&vdwl_smooth,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&coul_smooth,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussLong::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,
            buck6d1[i][i],buck6d2[i][i],buck6d3[i][i],
            buck6d4[i][i],alpha_ij[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBuck6dCoulGaussLong::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g\n",i,j,
              buck6d1[i][j],buck6d2[i][j],buck6d3[i][j],
              buck6d4[i][j],alpha_ij[i][j],cut_lj[i][j]);
}
/* ---------------------------------------------------------------------- */

double PairBuck6dCoulGaussLong::single(int i, int j, int itype, int jtype, double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r,r2inv,r6inv,r14inv,rexp,term1,term2,term3,term4,term5;
  double rcu,rqu,sme,smf;
  double erfa,expa,prefactor,arg,falpha,ealpha;
  double grij,expm2,erf;
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
  // long range - real space
    grij = g_ewald * r;
    expm2 = MathSpecial::expmsq(grij);
    erf = 1 - (MathSpecial::my_erfcx(grij) * expm2);

    // gaussian for 1/r alpha_ij contribution
    arg = alpha_ij[itype][jtype]*r;
    expa = MathSpecial::expmsq(arg);
    erfa = 1 - (MathSpecial::my_erfcx(arg) * expa);

    prefactor = force->qqrd2e * atom->q[i] * atom->q[j] / r;
    falpha = erfa - EWALD_F*arg*expa;
    forcecoul = prefactor * (falpha - erf + EWALD_F*grij*expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor*falpha;

    ealpha = prefactor * (erfa-erf);
    // smoothing term
    if (rsq > rsmooth_sq_c) {
      rcu = r*rsq;
      rqu = rsq*rsq;
      sme = c5_c*rqu*r + c4_c*rqu + c3_c*rcu + c2_c*rsq + c1_c*r + c0_c;
      smf = 5.0*c5_c*rqu + 4.0*c4_c*rcu + 3.0*c3_c*rsq + 2.0*c2_c*r + c1_c;
      forcecoul = forcecoul*sme - ealpha*smf*r;
      ealpha *= sme;
    }
  } else forcecoul = 0.0;

  fforce = (forcecoul + factor_lj*forcebuck6d) * r2inv;

  double eng = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    phibuck6d = ebuck6d - offset[itype][jtype];
    eng += factor_lj*phibuck6d;
  }

  if (rsq < cut_coulsq) {
    phicoul = ealpha;
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor*erfa;
    eng += phicoul;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairBuck6dCoulGaussLong::extract(const char *str, int &dim)
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
