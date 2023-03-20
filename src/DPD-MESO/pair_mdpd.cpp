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
   Contributing author: Zhen Li (Clemson University)
   Email: zli7@clemson.edu
------------------------------------------------------------------------- */

#include "pair_mdpd.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

#define EPSILON 1.0e-10

static const char cite_pair_mdpd[] =
  "pair mdpd command: doi:10.1063/1.4812366\n\n"
  "@Article{ZLi2013_POF,\n"
  " author = {Li, Z. and Hu, G. H. and Wang, Z. L. and Ma Y. B. and Zhou, Z. W.},\n"
  " title = {Three Dimensional Flow Structures in a Moving Droplet on Substrate: a Dissipative Particle Dynamics Study},\n"
  " journal = {Physics of Fluids},\n"
  " year = {2013},\n"
  " volume = {25},\n"
  " number = {7},\n"
  " pages = {072103}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairMDPD::PairMDPD(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_mdpd);

  writedata = 1;
  random = nullptr;
}

/* ---------------------------------------------------------------------- */

PairMDPD::~PairMDPD()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_r);
    memory->destroy(A_att);
    memory->destroy(B_rep);
    memory->destroy(gamma);
    memory->destroy(sigma);
  }
  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairMDPD::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wc,wc_r, wr,randnum,factor_dpd;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rhoi, rhoj;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rho= atom->rho;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);

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
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    rhoi = rho[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in MDPD systems
        rinv = 1.0/r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx*delvx + dely*delvy + delz*delvz;

        wc = 1.0 - r/cut[itype][jtype];
        wc_r = 1.0 - r/cut_r[itype][jtype];
        wc_r = MAX(wc_r,0.0);
        wr = wc;

        rhoj = rho[j];
        randnum = random->gaussian();

        // conservative force = A_att * wc + B_rep*(rhoi+rhoj)*wc_r
        // drag force = -gamma * wr^2 * (delx dot delv) / r
        // random force = sigma * wr * rnd * dtinvsqrt;

        fpair = A_att[itype][jtype]*wc + B_rep[itype][jtype]*(rhoi+rhoj)*wc_r;
        fpair -= gamma[itype][jtype]*wr*wr*dot*rinv;
        fpair += sigma[itype][jtype]*wr*randnum*dtinvsqrt;
        fpair *= factor_dpd*rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          // unshifted eng of conservative term:
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*A_att[itype][jtype]*cut[itype][jtype] * wr*wr + 0.5*B_rep[itype][jtype]*cut_r[itype][jtype]*(rhoi+rhoj)*wc_r*wc_r;
          evdwl *= factor_dpd;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMDPD::allocate()
{
  int i,j;
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_r,n+1,n+1,"pair:cut_r");
  memory->create(A_att,n+1,n+1,"pair:A_att");
  memory->create(B_rep,n+1,n+1,"pair:B_rep");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(sigma,n+1,n+1,"pair:sigma");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMDPD::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  temperature = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);
  seed = utils::inumeric(FLERR,arg[2],false,lmp);

  // initialize Marsaglia RNG with processor-unique seed
  // create a positive seed based on the system clock, if requested.

  if (seed <= 0) {
    constexpr double LARGE_NUM = 2<<30;
    seed = int(fmod(platform::walltime() * LARGE_NUM, LARGE_NUM)) + 1;
  }

  delete random;
  random = new RanMars(lmp,(seed + comm->me) % 900000000);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMDPD::coeff(int narg, char **arg)
{
  if (narg != 7 ) error->all(FLERR,"Incorrect args for pair coefficients\n itype jtype A B gamma cutA cutB");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double A_one = utils::numeric(FLERR,arg[2],false,lmp);
  double B_one = utils::numeric(FLERR,arg[3],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[4],false,lmp);
  double cut_one = utils::numeric(FLERR,arg[5],false,lmp);
  double cut_two = utils::numeric(FLERR,arg[6],false,lmp);

  if (cut_one < cut_two) error->all(FLERR,"Incorrect args for pair coefficients\n cutA should be larger than cutB.");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      A_att[i][j] = A_one;
      B_rep[i][j] = B_one;
      gamma[i][j] = gamma_one;
      cut[i][j] = cut_one;
      cut_r[i][j] = cut_two;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairMDPD::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair mdpd requires ghost atoms store velocity");

  if (!atom->rho_flag)
    error->all(FLERR,"Pair style mdpd requires atom attribute rho");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR, "Pair mdpd needs newton pair on for momentum conservation");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMDPD::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);

  cut[j][i] = cut[i][j];
  cut_r[j][i] = cut_r[i][j];
  A_att[j][i] = A_att[i][j];
  B_rep[j][i] = B_rep[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMDPD::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&A_att[i][j],sizeof(double),1,fp);
        fwrite(&B_rep[i][j],sizeof(double),1,fp);
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
        fwrite(&cut_r[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMDPD::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&A_att[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&B_rep[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_r[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&A_att[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&B_rep[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_r[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMDPD::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMDPD::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&temperature,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&temperature,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMDPD::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,A_att[i][i],B_rep[i][i],gamma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMDPD::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,A_att[i][j],B_rep[i][j],gamma[i][j],cut[i][j],cut_r[i][j]);
}

