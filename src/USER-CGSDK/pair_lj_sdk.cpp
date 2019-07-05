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
   Contributing author: Axel Kohlmeyer (Temple U)
   This style is a simplified re-implementation of the CG/CMM pair style
------------------------------------------------------------------------- */

#include "pair_lj_sdk.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#define LMP_NEED_SDK_FIND_LJ_TYPE 1
#include "lj_sdk_common.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace LJSDKParms;

/* ---------------------------------------------------------------------- */

PairLJSDK::PairLJSDK(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 0;
  single_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJSDK::~PairLJSDK()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(lj_type);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);

    memory->destroy(rminsq);
    memory->destroy(emin);

    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJSDK::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) eval<1,1,1>();
      else eval<1,1,0>();
    } else {
      if (force->newton_pair) eval<1,0,1>();
      else eval<1,0,0>();
    }
  } else {
    if (force->newton_pair) eval<0,0,1>();
    else eval<0,0,0>();
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJSDK::eval()
{
  int i,j,ii,jj,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,forcelj,factor_lj;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = atom->f;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_lj = force->special_lj;
  double fxtmp,fytmp,fztmp;

  const int inum = list->inum;
  const int * const ilist = list->ilist;
  const int * const numneigh = list->numneigh;
  const int * const * const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp=fytmp=fztmp=0.0;

    const int itype = type[i];
    const int * const jlist = firstneigh[i];
    const int jnum = numneigh[i];

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
        const int ljt = lj_type[itype][jtype];

        if (ljt == LJ12_4) {
          const double r4inv=r2inv*r2inv;
          forcelj = r4inv*(lj1[itype][jtype]*r4inv*r4inv
                           - lj2[itype][jtype]);

          if (EFLAG)
            evdwl = r4inv*(lj3[itype][jtype]*r4inv*r4inv
                           - lj4[itype][jtype]) - offset[itype][jtype];

        } else if (ljt == LJ9_6) {
          const double r3inv = r2inv*sqrt(r2inv);
          const double r6inv = r3inv*r3inv;
          forcelj = r6inv*(lj1[itype][jtype]*r3inv
                           - lj2[itype][jtype]);
          if (EFLAG)
            evdwl = r6inv*(lj3[itype][jtype]*r3inv
                           - lj4[itype][jtype]) - offset[itype][jtype];

        } else if (ljt == LJ12_6) {
          const double r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv*(lj1[itype][jtype]*r6inv
                          - lj2[itype][jtype]);
          if (EFLAG)
            evdwl = r6inv*(lj3[itype][jtype]*r6inv
                           - lj4[itype][jtype]) - offset[itype][jtype];
        } else continue;

        fpair = factor_lj*forcelj*r2inv;

        fxtmp += delx*fpair;
        fytmp += dely*fpair;
        fztmp += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (EFLAG) evdwl *= factor_lj;
        if (EVFLAG) ev_tally(i,j,nlocal,NEWTON_PAIR,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJSDK::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(lj_type,n+1,n+1,"pair:lj_type");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
      lj_type[i][j] = LJ_NOT_SET;
    }
  }

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");

  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");

  memory->create(offset,n+1,n+1,"pair:offset");

  memory->create(rminsq,n+1,n+1,"pair:rminsq");
  memory->create(emin,n+1,n+1,"pair:emin");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJSDK::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

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

void PairLJSDK::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  int lj_type_one = find_lj_type(arg[2],lj_type_list);
  if (lj_type_one == LJ_NOT_SET)
    error->all(FLERR,"Cannot parse LJ type flag.");

  double epsilon_one = force->numeric(FLERR,arg[3]);
  double sigma_one = force->numeric(FLERR,arg[4]);

  double cut_one = cut_global;
  if (narg == 6) cut_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      lj_type[i][j] = lj_type_one;
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
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

double PairLJSDK::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR,"No mixing support for lj/sdk. "
               "Coefficients for all pairs need to be set explicitly.");

  const int ljt = lj_type[i][j];

  if (ljt == LJ_NOT_SET)
    error->all(FLERR,"unrecognized LJ parameter flag");

  lj1[i][j] = lj_prefact[ljt] * lj_pow1[ljt] * epsilon[i][j] * pow(sigma[i][j],lj_pow1[ljt]);
  lj2[i][j] = lj_prefact[ljt] * lj_pow2[ljt] * epsilon[i][j] * pow(sigma[i][j],lj_pow2[ljt]);
  lj3[i][j] = lj_prefact[ljt] * epsilon[i][j] * pow(sigma[i][j],lj_pow1[ljt]);
  lj4[i][j] = lj_prefact[ljt] * epsilon[i][j] * pow(sigma[i][j],lj_pow2[ljt]);

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = lj_prefact[ljt] * epsilon[i][j] * (pow(ratio,lj_pow1[ljt]) - pow(ratio,lj_pow2[ljt]));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  cut[j][i] = cut[i][j];
  cutsq[j][i] = cutsq[i][j];
  offset[j][i] = offset[i][j];
  lj_type[j][i] = lj_type[i][j];

  // compute derived parameters for SDK angle potential

  const double eps = epsilon[i][j];
  const double sig = sigma[i][j];
  const double rmin = sig*exp(1.0/(lj_pow1[ljt]-lj_pow2[ljt])
                              *log(lj_pow1[ljt]/lj_pow2[ljt]) );
  rminsq[j][i] = rminsq[i][j] = rmin*rmin;

  const double ratio = sig/rmin;
  const double emin_one = lj_prefact[ljt] * eps * (pow(ratio,lj_pow1[ljt])
                                                   - pow(ratio,lj_pow2[ljt]));
  emin[j][i] = emin[i][j] = emin_one;

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag)
    error->all(FLERR,"Tail flag not supported by lj/sdk pair style");

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJSDK::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&lj_type[i][j],sizeof(int),1,fp);
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJSDK::read_restart(FILE *fp)
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
          fread(&lj_type[i][j],sizeof(int),1,fp);
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&lj_type[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJSDK::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJSDK::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
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
   lj/sdk does not support per atom type output with mixing
------------------------------------------------------------------------- */

void PairLJSDK::write_data(FILE *)
{
  error->one(FLERR, "Pair style lj/sdk requires using "
             "write_data with the 'pair ij' option");
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJSDK::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %s %g %g %g\n",i,j,lj_type_list[lj_type[i][j]],
              epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJSDK::single(int, int, int itype, int jtype, double rsq,
                         double, double factor_lj, double &fforce)
{

  if (rsq < cutsq[itype][jtype]) {

    const int ljt = lj_type[itype][jtype];
    const double ljpow1 = lj_pow1[ljt];
    const double ljpow2 = lj_pow2[ljt];
    const double ljpref = lj_prefact[ljt];

    const double ratio = sigma[itype][jtype]/sqrt(rsq);
    const double eps = epsilon[itype][jtype];

    fforce = factor_lj * ljpref*eps * (ljpow1*pow(ratio,ljpow1)
                          - ljpow2*pow(ratio,ljpow2))/rsq;
    return factor_lj * (ljpref*eps * (pow(ratio,ljpow1) - pow(ratio,ljpow2))
                        - offset[itype][jtype]);

  } else fforce=0.0;

  return 0.0;
}

/* ---------------------------------------------------------------------- */

void *PairLJSDK::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"lj_type") == 0) return (void *) lj_type;
  if (strcmp(str,"lj1") == 0) return (void *) lj1;
  if (strcmp(str,"lj2") == 0) return (void *) lj2;
  if (strcmp(str,"lj3") == 0) return (void *) lj3;
  if (strcmp(str,"lj4") == 0) return (void *) lj4;
  if (strcmp(str,"rminsq") == 0) return (void *) rminsq;
  if (strcmp(str,"emin") == 0) return (void *) emin;
  return NULL;
}

/* ---------------------------------------------------------------------- */

double PairLJSDK::memory_usage()
{
  double bytes = Pair::memory_usage();
  int n = atom->ntypes;

  // setflag/lj_type
  bytes += 2 * (n+1)*(n+1)*sizeof(int);
  // cut/cutsq/epsilon/sigma/offset/lj1/lj2/lj3/lj4/rminsq/emin
  bytes += 11 * (n+1)*(n+1)*sizeof(double);

  return bytes;
}
