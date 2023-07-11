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
   Contributing Author: Xipeng Wang, Simon Ramirez-Hinestrosa
------------------------------------------------------------------------- */

#include "pair_wf_cut.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

PairWFCut::PairWFCut(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairWFCut::~PairWFCut()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(nu);
    memory->destroy(mu);
    memory->destroy(nm);
    memory->destroy(e0nm);  //Alpha * epsilon
    memory->destroy(rcmu);
    memory->destroy(sigma_mu);
  }
}

/* ---------------------------------------------------------------------- */

void PairWFCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,factor_lj;
  double forcenm,rminv, rm, rn;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);
  evdwl = 0.0;

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
        rminv = powint(r2inv,mu[itype][jtype]);
        rm = sigma_mu[itype][jtype]*rminv - 1.0;
        rn = rcmu[itype][jtype]*rminv - 1.0;

        forcenm = 2.0*mu[itype][jtype] *sigma_mu[itype][jtype]*powint(rn,2*nu[itype][jtype])
                + 4.0*nm[itype][jtype] *rcmu[itype][jtype]*rm*powint(rn,2*nu[itype][jtype]-1);
        fpair = factor_lj*e0nm[itype][jtype]*forcenm*powint(r2inv,mu[itype][jtype]+1);

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          evdwl = e0nm[itype][jtype] * (rm*powint(rn,2*nu[itype][jtype]));
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

void PairWFCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(nu,n+1,n+1,"pair:nu");
  memory->create(mu,n+1,n+1,"pair:mu");
  memory->create(nm,n+1,n+1,"pair:nm");
  memory->create(e0nm,n+1,n+1,"pair:e0nm");
  memory->create(rcmu,n+1,n+1,"pair:rcmu");
  memory->create(sigma_mu,n+1,n+1,"pair:sigma_mu");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairWFCut::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

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

void PairWFCut::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);
  int nu_one = utils::inumeric(FLERR,arg[4],false,lmp);
  int mu_one = utils::inumeric(FLERR,arg[5],false,lmp);

  double cut_one = cut_global;
  if (narg == 7) cut_one = utils::numeric(FLERR,arg[6],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      nu[i][j] = nu_one;
      mu[i][j] = mu_one;
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

double PairWFCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  nm[i][j] = nu[i][j]*mu[i][j];
  e0nm[i][j] = epsilon[i][j]*2.0*nu[i][j]*powint(cut[i][j]/sigma[i][j],2*mu[i][j])
                       *powint((1+2.0*nu[i][j])/(2.0*nu[i][j])/(MathSpecial::powint(cut[i][j]/sigma[i][j],2*mu[i][j])-1.0),
                              2*nu[i][j]+1);
  rcmu[i][j] = powint(cut[i][j],2*mu[i][j]);
  sigma_mu[i][j] = powint(sigma[i][j], 2*mu[i][j]);

  epsilon[j][i] = epsilon[i][j];
  nu[j][i] = nu[i][j];
  mu[j][i] = mu[i][j];
  nm[j][i] = nm[i][j];
  sigma[j][i] = sigma[i][j];
  e0nm[j][i] = e0nm[i][j];
  rcmu[j][i] = rcmu[i][j];
  sigma_mu[j][i] = sigma_mu[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWFCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&nu[i][j],sizeof(int),1,fp);
        fwrite(&mu[i][j],sizeof(int),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWFCut::read_restart(FILE *fp)
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
        if (me == 0) {
          utils::sfread(FLERR, &epsilon[i][j],sizeof(double),1,fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j],sizeof(double),1,fp, nullptr, error);
          utils::sfread(FLERR, &nu[i][j],sizeof(int),1,fp, nullptr, error);
          utils::sfread(FLERR, &mu[i][j],sizeof(int),1,fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j],sizeof(double),1,fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&nu[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&mu[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWFCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWFCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR, &mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairWFCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %d %d\n",i,epsilon[i][i],sigma[i][i],nu[i][i],mu[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairWFCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %d %d %g\n",i,j,
              epsilon[i][j],sigma[i][j],nu[i][j],mu[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairWFCut::single(int /*i*/, int /*j*/, int itype, int jtype,
                      double rsq, double /*factor_coul*/, double factor_lj,
                      double &fforce)
{
  double r2inv,rminv,rm,rn,forcenm,phinm;

  r2inv = 1.0/rsq;
  rminv =powint(r2inv,mu[itype][jtype]);
  rm = sigma_mu[itype][jtype]*rminv - 1.0;
  rn = rcmu[itype][jtype]*rminv - 1.0;
  forcenm = 2.0*mu[itype][jtype] *sigma_mu[itype][jtype]*powint(rn,2*nu[itype][jtype])
                + 4.0*nm[itype][jtype] *rcmu[itype][jtype]*rm*powint(rn,2*nu[itype][jtype]-1);
  fforce = factor_lj*e0nm[itype][jtype]*forcenm*powint(r2inv,mu[itype][jtype]+1);

  phinm = e0nm[itype][jtype] * rm*powint(rn,2*nu[itype][jtype]);
  return factor_lj*phinm;
}

/* ---------------------------------------------------------------------- */

void *PairWFCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"nu") == 0) return (void *) nu;
  if (strcmp(str,"mu") == 0) return (void *) mu;
  return nullptr;
}
