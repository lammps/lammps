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
   Contributing Author: Julien Devemy (ICCF)
------------------------------------------------------------------------- */

#include "pair_nm_cut_coul_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairNMCutCoulLong::PairNMCutCoulLong(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  ftable = nullptr;
  qdist = 0.0;
}

/* ---------------------------------------------------------------------- */

PairNMCutCoulLong::~PairNMCutCoulLong()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
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
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

void PairNMCutCoulLong::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,r2inv,factor_coul,factor_lj;
  double forcecoul,forcenm,rminv,rninv;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

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

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r = sqrt(rsq);
          rminv = pow(r2inv,mm[itype][jtype]/2.0);
          rninv = pow(r2inv,nn[itype][jtype]/2.0);
          forcenm = e0nm[itype][jtype]*nm[itype][jtype] *
            (r0n[itype][jtype]/pow(r,nn[itype][jtype]) -
             r0m[itype][jtype]/pow(r,mm[itype][jtype]));
        } else forcenm = 0.0;

        fpair = (forcecoul + factor_lj*forcenm) * r2inv;

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
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = e0nm[itype][jtype] *
              (mm[itype][jtype]*r0n[itype][jtype]*rninv -
               nn[itype][jtype]*r0m[itype][jtype]*rminv) - offset[itype][jtype];
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

void PairNMCutCoulLong::allocate()
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

void PairNMCutCoulLong::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_lj_global = utils::numeric(FLERR,arg[0],false,lmp);
  if (narg == 1) cut_coul = cut_lj_global;
  else cut_coul = utils::numeric(FLERR,arg[1],false,lmp);

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

void PairNMCutCoulLong::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double e0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double r0_one = utils::numeric(FLERR,arg[3],false,lmp);
  double nn_one = utils::numeric(FLERR,arg[4],false,lmp);
  double mm_one = utils::numeric(FLERR,arg[5],false,lmp);

  double cut_lj_one = cut_lj_global;
  if (narg == 7) cut_lj_one = utils::numeric(FLERR,arg[6],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      e0[i][j] = e0_one;
      r0[i][j] = r0_one;
      nn[i][j] = nn_one;
      mm[i][j] = mm_one;
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

void PairNMCutCoulLong::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style nm/cut/coul/long requires atom attribute q");

  neighbor->add_request(this);

  cut_coulsq = cut_coul * cut_coul;
  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul,nullptr);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairNMCutCoulLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  // include TIP4P qdist in full cutoff, qdist = 0.0 if not TIP4P

  double cut = MAX(cut_lj[i][j],cut_coul+2.0*qdist);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

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

void PairNMCutCoulLong::write_restart(FILE *fp)
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
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNMCutCoulLong::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&e0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&r0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&nn[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&mm[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_lj[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&e0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&nn[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&mm[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairNMCutCoulLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
  fwrite(&ncoultablebits,sizeof(int),1,fp);
  fwrite(&tabinner,sizeof(double),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairNMCutCoulLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_lj_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&ncoultablebits,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&tabinner,sizeof(double),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
  MPI_Bcast(&ncoultablebits,1,MPI_INT,0,world);
  MPI_Bcast(&tabinner,1,MPI_DOUBLE,0,world);
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairNMCutCoulLong::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,e0[i][i],r0[i][i],nn[i][i],mm[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairNMCutCoulLong::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,
              e0[i][j],r0[i][j],nn[i][j],mm[i][j],cut_lj[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairNMCutCoulLong::single(int i, int j, int itype, int jtype,
                                 double rsq,
                                 double factor_coul, double factor_lj,
                                 double &fforce)
{
  double r2inv,r,grij,expm2,t,erfc,prefactor;
  double fraction,table,forcecoul,forcenm,phicoul,phinm;
  int itable;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq) {
      r = sqrt(rsq);
      grij = g_ewald * r;
      expm2 = exp(-grij*grij);
      t = 1.0 / (1.0 + EWALD_P*grij);
      erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
      prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
      forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
      if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
    } else {
      union_int_float_t rsq_lookup_single;
      rsq_lookup_single.f = rsq;
      itable = rsq_lookup_single.i & ncoulmask;
      itable >>= ncoulshiftbits;
      fraction = (rsq_lookup_single.f - rtable[itable]) * drtable[itable];
      table = ftable[itable] + fraction*dftable[itable];
      forcecoul = atom->q[i]*atom->q[j] * table;
      if (factor_coul < 1.0) {
        table = ctable[itable] + fraction*dctable[itable];
        prefactor = atom->q[i]*atom->q[j] * table;
        forcecoul -= (1.0-factor_coul)*prefactor;
      }
    }
  } else forcecoul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {
    r = sqrt(rsq);
    forcenm = e0nm[itype][jtype]*nm[itype][jtype] *
      (r0n[itype][jtype]/pow(r,nn[itype][jtype]) -
       r0m[itype][jtype]/pow(r,mm[itype][jtype]));
  } else forcenm = 0.0;

  fforce = (forcecoul + factor_lj*forcenm) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq)
      phicoul = prefactor*erfc;
    else {
      table = etable[itable] + fraction*detable[itable];
      phicoul = atom->q[i]*atom->q[j] * table;
    }
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    eng += phicoul;
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

void *PairNMCutCoulLong::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  dim = 2;
  if (strcmp(str,"e0") == 0) return (void *) e0;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"nn") == 0) return (void *) nn;
  if (strcmp(str,"mm") == 0) return (void *) mm;
  return nullptr;
}
