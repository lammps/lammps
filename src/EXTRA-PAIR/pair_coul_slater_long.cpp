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
 *     Contributing author:  Evangelos Voyiatzis (Royal DSM)
 * ------------------------------------------------------------------------- */

#include "pair_coul_slater_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairCoulSlaterLong::PairCoulSlaterLong(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  qdist = 0.0;
}

/* ---------------------------------------------------------------------- */

PairCoulSlaterLong::~PairCoulSlaterLong()
{
  if (!copymode) {
    if (allocated) {
      memory->destroy(setflag);
      memory->destroy(cutsq);

      memory->destroy(scale);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulSlaterLong::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,r2inv,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;
  double slater_term;

  ev_init(eflag,vflag);
  ecoul = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
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
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        r2inv = 1.0/rsq;
        r = sqrt(rsq);
        grij = g_ewald * r;
        expm2 = exp(-grij*grij);
        t = 1.0 / (1.0 + EWALD_P*grij);
        erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
        slater_term = exp(-2*r/lamda)*(1 + (2*r/lamda*(1+r/lamda)));
        prefactor = qqrd2e * scale[itype][jtype] * qtmp*q[j]/r;
        forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - slater_term);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor*(1-slater_term);

        fpair = forcecoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          ecoul = prefactor*(erfc - (1 + r/lamda)*exp(-2*r/lamda));
          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor*(1.0-(1 + r/lamda)*exp(-2*r/lamda));
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

void PairCoulSlaterLong::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(scale,n+1,n+1,"pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulSlaterLong::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  lamda = utils::numeric(FLERR,arg[0],false,lmp);
  cut_coul = utils::numeric(FLERR,arg[1],false,lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulSlaterLong::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
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

void PairCoulSlaterLong::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/slater/long requires atom attribute q");

  neighbor->add_request(this);

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

 if (force->kspace == nullptr)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  // if (ncoultablebits) init_tables(cut_coul,nullptr);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulSlaterLong::init_one(int i, int j)
{
  scale[j][i] = scale[i][j];
  return cut_coul+2.0*qdist;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulSlaterLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j])
        fwrite(&scale[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulSlaterLong::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j],sizeof(int),1,fp, nullptr, error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) utils::sfread(FLERR, &scale[i][j],sizeof(double),1,fp, nullptr, error);
        MPI_Bcast(&scale[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulSlaterLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&lamda,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulSlaterLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_coul,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR, &lamda,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR, &offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR, &mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&lamda,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulSlaterLong::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                                  double factor_coul, double /*factor_lj*/, double &fforce)
{
  double r2inv,r,grij,expm2,t,erfc,prefactor;
  double slater_term;
  double forcecoul,phicoul;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  grij = g_ewald * r;
  expm2 = exp(-grij*grij);
  t = 1.0 / (1.0 + EWALD_P*grij);
  erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
  slater_term = exp(-2*r/lamda)*(1 + (2*r/lamda*(1+r/lamda)));
  prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
  forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - slater_term);
  if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  fforce = forcecoul * r2inv;

  phicoul = prefactor*(erfc - (1 + r/lamda)*exp(-2*r/lamda));
  if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;

  return phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulSlaterLong::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  if (strcmp(str,"lamda") == 0) {
    dim = 0;
    return (void *) &lamda;
  }
  if (strcmp(str,"scale") == 0) {
    dim = 2;
    return (void *) scale;
  }
  return nullptr;
}
