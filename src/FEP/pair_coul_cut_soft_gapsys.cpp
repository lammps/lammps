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

#include "pair_coul_cut_soft_gapsys.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulCutSoftGapsys::PairCoulCutSoftGapsys(LAMMPS *lmp) : Pair(lmp), cut(nullptr), lambda(nullptr)
{
  centroidstressflag = CENTROID_SAME;
  cut_global = sigmaq = alphaq = 0.0;
}

/* ---------------------------------------------------------------------- */

PairCoulCutSoftGapsys::~PairCoulCutSoftGapsys()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(lambda);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,forcecoul,factor_coul;
  double cut_inner;
  int *ilist,*jlist,*numneigh,**firstneigh;

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

      cut_inner = (1.0 + sigmaq * fabs(qtmp * q[j])) * alphaq * pow(lambda[itype][jtype], 1.0 / 6.0);

      if (rsq < cut_inner * cut_inner) {

        forcecoul = factor_coul * qqrd2e * qtmp * q[j];
        fpair = (- 2.0 / pow(cut_inner, 3) + 3.0 / (pow(cut_inner, 2) * sqrt(rsq))) * forcecoul;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag)
          ecoul = factor_coul * qqrd2e * qtmp * q[j] *
          (rsq / pow(cut_inner, 3) - 3 * sqrt(rsq) / pow(cut_inner, 2) + 3 / cut_inner);

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);

      } else if (rsq < cutsq[itype][jtype]) {

        double r2inv = 1.0 / rsq;
        double rinv = sqrt(r2inv);
        forcecoul = qqrd2e * qtmp * q[j] * rinv;
        fpair = factor_coul * forcecoul * r2inv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag)
          ecoul = factor_coul * qqrd2e * qtmp * q[j] * rinv;

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

void PairCoulCutSoftGapsys::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(lambda,n+1,n+1,"pair:lambda");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  sigmaq = utils::numeric(FLERR,arg[0],false,lmp);
  alphaq = utils::numeric(FLERR,arg[1],false,lmp);

  cut_global = utils::numeric(FLERR,arg[2],false,lmp);

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

void PairCoulCutSoftGapsys::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double lambda_one = utils::numeric(FLERR,arg[2],false,lmp);

  double cut_one = cut_global;
  if (narg == 4) cut_one = utils::numeric(FLERR,arg[3],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      lambda[i][j] = lambda_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/cut/soft/gapsys requires atom attribute q");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulCutSoftGapsys::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    if (lambda[i][i] != lambda[j][j])
      error->all(FLERR,"Pair coul/cut/soft/gapsys different lambda values in mix");
    lambda[i][j] = lambda[i][i];
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  cut[j][i] = cut[i][j];
  lambda[j][i] = lambda[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&lambda[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&lambda[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&lambda[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::write_restart_settings(FILE *fp)
{
  fwrite(&sigmaq,sizeof(double),1,fp);
  fwrite(&alphaq,sizeof(double),1,fp);

  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&sigmaq,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&alphaq,sizeof(double),1,fp,nullptr,error);

    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&sigmaq,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alphaq,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g\n",i,lambda[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairCoulCutSoftGapsys::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g\n",i,j,lambda[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairCoulCutSoftGapsys::single(int i, int j, int itype, int jtype,
                           double rsq, double factor_coul, double /*factor_lj*/,
                           double &fforce)
{
  double cut_inner;

  if (rsq > cutsq[itype][jtype]) {
    fforce = 0.0;
    return 0.0;
  }

  cut_inner = (1.0 + sigmaq * fabs(atom->q[i] * atom->q[j])) * alphaq * pow(lambda[itype][jtype], 1.0 / 6.0);

  fforce = factor_coul * force->qqrd2e * atom->q[i] * atom->q[j];
  double phicoul = factor_coul * force->qqrd2e * atom->q[i] * atom->q[j];

  if (rsq > cut_inner * cut_inner) {

    double r2inv = 1.0 / rsq;
    double rinv = sqrt(r2inv);
    fforce *= rinv * r2inv;
    phicoul *= rinv;

  } else {

    fforce *= - 2.0 / pow(cut_inner, 3) + 3.0 / (pow(cut_inner, 2) * sqrt(rsq));
    phicoul *= rsq / pow(cut_inner, 3) - 3 * sqrt(rsq) / pow(cut_inner, 2) + 3 / cut_inner;

  }

  return phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulCutSoftGapsys::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  return nullptr;
}
