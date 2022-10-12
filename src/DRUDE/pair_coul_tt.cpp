// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_coul_tt.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "fix.h"
#include "fix_drude.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "modify.h"
#include "error.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulTT::PairCoulTT(LAMMPS *lmp) : Pair(lmp) {
    fix_drude = nullptr;
}

/* ---------------------------------------------------------------------- */

PairCoulTT::~PairCoulTT()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(b);
    memory->destroy(c);
    memory->destroy(ntt);
    memory->destroy(cut);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulTT::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qi,qj,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,rsq,r2inv,rinv,factor_coul;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_f,factor_e;
  int di,dj;
  double dcoul;
  double beta, gamma, betaprime, gammaprime, gammatmp;

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

    qi = q[i];

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

      if (drudetype[type[i]] == drudetype[type[j]] && drudetype[type[j]] != CORE_TYPE)
        continue;

      qj = q[j];

      if (drudetype[type[i]] == CORE_TYPE) {
        di = domain->closest_image(i, atom->map(drudeid[i]));
        if (di == j)
          continue;
        switch (drudetype[type[j]]) {
          case DRUDE_TYPE:
            qi = q[i]+q[di];
            break;
          case NOPOL_TYPE:
            qi = -q[di];
            break;
        }
      }

      if (drudetype[type[j]] == CORE_TYPE) {
        dj = domain->closest_image(j, atom->map(drudeid[j]));
        if (dj == i)
          continue;
        switch (drudetype[type[i]]) {
          case DRUDE_TYPE:
            qj = q[j]+q[dj];
            break;
          case NOPOL_TYPE:
            qj = -q[dj];
            break;
        }
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
        beta = c[itype][jtype] * exp(-b[itype][jtype] * r);
        betaprime = -b[itype][jtype] * beta;
        gamma = 1.0 + b[itype][jtype] * r;
        gammaprime = b[itype][jtype];
        gammatmp = 1.0;
        for (int k = 2; k <= ntt[itype][jtype]; k++) {
          gammatmp *= b[itype][jtype] * r / k;
          gamma += gammatmp * b[itype][jtype] * r;
          gammaprime += gammatmp * b[itype][jtype] * k;
        }

        if (drudetype[type[i]] == CORE_TYPE && drudetype[type[j]] == CORE_TYPE)
          dcoul = qqrd2e * ( -(q[i]+q[di])*q[dj] - q[di]*(q[j]+q[dj]) ) * scale[itype][jtype] * rinv;
        else
          dcoul = qqrd2e * qi * qj *scale[itype][jtype] * rinv;

        factor_f = (-beta*gamma + r*betaprime*gamma + r*beta*gammaprime)*factor_coul;
        if (eflag) factor_e = - beta*gamma*factor_coul;
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

void PairCoulTT::allocate()
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
  memory->create(b, n + 1, n + 1, "pair:b");
  memory->create(c, n + 1, n + 1, "pair:c");
  memory->create(ntt, n + 1, n + 1, "pair:ntt");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulTT::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal pair_style command");

  n_global = utils::inumeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
            ntt[i][j] = n_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulTT::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;

  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double b_one = utils::numeric(FLERR,arg[2],false,lmp);
  double c_one = utils::numeric(FLERR,arg[3],false,lmp);
  int n_one = n_global;
  double cut_one = cut_global;
  if (narg >= 5) n_one = utils::inumeric(FLERR,arg[4],false,lmp);
  if (narg == 6) cut_one = utils::numeric(FLERR,arg[5],false,lmp);

  if (n_one > n_global)
    error->all(FLERR, "Incorrect coefficients for pair style coul/tt: n should not be larger than global setting");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      b[i][j] = b_one;
      c[i][j] = c_one;
      ntt[i][j] = n_one;
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

void PairCoulTT::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/tt requires atom attribute q");
  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (utils::strmatch(modify->fix[ifix]->style,"^drude")) break;
  if (ifix == modify->nfix) error->all(FLERR, "Pair coul/tt requires fix drude");
  fix_drude = dynamic_cast<FixDrude *>(modify->fix[ifix]);

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulTT::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);

  b[j][i] = b[i][j];
  c[j][i] = c[i][j];
  ntt[j][i] = ntt[i][j];
  scale[j][i] = scale[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulTT::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&b[i][j], sizeof(double), 1, fp);
        fwrite(&c[i][j], sizeof(double), 1, fp);
        fwrite(&ntt[i][j], sizeof(int), 1, fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulTT::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&b[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&c[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&ntt[i][j],sizeof(int),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          }
        MPI_Bcast(&b[i][j],  1, MPI_DOUBLE,0,world);
        MPI_Bcast(&c[i][j],  1, MPI_DOUBLE,0,world);
        MPI_Bcast(&ntt[i][j],  1, MPI_INT,0,world);
        MPI_Bcast(&cut[i][j],1, MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulTT::write_restart_settings(FILE *fp)
{
  fwrite(&n_global, sizeof(int), 1, fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulTT::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&n_global,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&n_global, 1, MPI_INT, 0, world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulTT::single(int i, int j, int itype, int jtype,
                         double rsq, double factor_coul, double /*factor_lj*/,
                         double &fforce)
{
  double r2inv,rinv,r,phicoul;
  double qi,qj,factor_f,factor_e,dcoul;
  double beta, betaprime, gamma, gammaprime, gammatmp;

  // single() is used in Pair::write_file where there is no information about topology
  // it corresponds to pair_write command in input file
  // charges qi and qj are defined by the user (or 1.0 by defaut)

  qi = atom->q[i];
  qj = atom->q[j];

  r2inv = 1.0/rsq;
  fforce = phicoul = 0.0;
  if (rsq < cutsq[itype][jtype]) {
    rinv = sqrt(r2inv);
    r = sqrt(rsq);
    beta = c[itype][jtype] * exp(-b[itype][jtype] * r);
    betaprime = -b[itype][jtype] * beta;
    gamma = 1 + b[itype][jtype] * r;
    gammaprime = b[itype][jtype];
    gammatmp = 1;
    for (int k = 2; k <= ntt[itype][jtype]; k++) {
      gammatmp *= b[itype][jtype] * r / k;
      gamma += gammatmp * b[itype][jtype] * r;
      gammaprime += gammatmp * b[itype][jtype] * k;
    }
    dcoul = force->qqrd2e * qi * qj * scale[itype][jtype] * rinv;
    factor_f = (-beta*gamma+r*betaprime*gamma+r*beta*gammaprime)*factor_coul;
    fforce = factor_f * dcoul * r2inv;
    factor_e = - beta*gamma*factor_coul;
    phicoul = factor_e * dcoul;
  }

  return phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulTT::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  if (strcmp(str,"b") == 0) return (void *) b;
  if (strcmp(str,"c") == 0) return (void *) c;
  if (strcmp(str,"ntt") == 0) return (void *) ntt;
  return nullptr;
}
