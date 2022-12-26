// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
  Contributing authors: Martin Svoboda (ICPF, UJEP), Martin Lísal (ICPF, UJEP)
  based on pair style dpd by: Kurt Smith (U Pittsburgh)
------------------------------------------------------------------------- */

#include "pair_dpd_ext.h"

#include "atom.h"
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

/* ---------------------------------------------------------------------- */

PairDPDExt::PairDPDExt(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  random = nullptr;
}

/* ---------------------------------------------------------------------- */

PairDPDExt::~PairDPDExt()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(a0);
    memory->destroy(gamma);
    memory->destroy(gammaT);
    memory->destroy(sigma);
    memory->destroy(sigmaT);
    memory->destroy(ws);
    memory->destroy(wsT);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDPDExt::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpairx,fpairy,fpairz,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,wdPar,wdPerp,randnum,randnumx,randnumy,randnumz;
  double prefactor_g,prefactor_s,factor_dpd,factor_sqrt;
  double P[3][3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
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

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_dpd = special_lj[sbmask(j)];
      factor_sqrt = special_sqrt[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        if (r < EPSILON) continue;     // r can be 0.0 in DPD systems
        rinv = 1.0/r;
        delvx = vxtmp - v[j][0];
        delvy = vytmp - v[j][1];
        delvz = vztmp - v[j][2];
        dot = delx*delvx + dely*delvy + delz*delvz;

        P[0][0] = 1.0 - delx*delx*rinv*rinv;
        P[0][1] =     - delx*dely*rinv*rinv;
        P[0][2] =     - delx*delz*rinv*rinv;

        P[1][0] = P[0][1];
        P[1][1] = 1.0 - dely*dely*rinv*rinv;
        P[1][2] =     - dely*delz*rinv*rinv;

        P[2][0] = P[0][2];
        P[2][1] = P[1][2];
        P[2][2] = 1.0 - delz*delz*rinv*rinv;

        wd = 1.0 - r/cut[itype][jtype];
        wdPar = pow(wd,ws[itype][jtype]);
        wdPerp = pow(wd,wsT[itype][jtype]);

        randnum = random->gaussian();
        randnumx = random->gaussian();
        randnumy = random->gaussian();
        randnumz = random->gaussian();

        // conservative force
        fpair = a0[itype][jtype]*wd;

        // drag force - parallel
        fpair -= gamma[itype][jtype]*wdPar*wdPar*dot*rinv;
        fpair *= factor_dpd;

        // random force - parallel
        fpair += factor_sqrt*sigma[itype][jtype]*wdPar*randnum*dtinvsqrt;

        fpairx = fpair*rinv*delx;
        fpairy = fpair*rinv*dely;
        fpairz = fpair*rinv*delz;

        // drag force - perpendicular
        prefactor_g = factor_dpd*gammaT[itype][jtype]*wdPerp*wdPerp;
        fpairx -= prefactor_g * (P[0][0]*delvx + P[0][1]*delvy + P[0][2]*delvz);
        fpairy -= prefactor_g * (P[1][0]*delvx + P[1][1]*delvy + P[1][2]*delvz);
        fpairz -= prefactor_g * (P[2][0]*delvx + P[2][1]*delvy + P[2][2]*delvz);

        // random force - perpendicular
        prefactor_s = factor_sqrt * sigmaT[itype][jtype]*wdPerp * dtinvsqrt;
        fpairx += prefactor_s * (P[0][0]*randnumx + P[0][1]*randnumy + P[0][2]*randnumz);
        fpairy += prefactor_s * (P[1][0]*randnumx + P[1][1]*randnumy + P[1][2]*randnumz);
        fpairz += prefactor_s * (P[2][0]*randnumx + P[2][1]*randnumy + P[2][2]*randnumz);

        f[i][0] += fpairx;
        f[i][1] += fpairy;
        f[i][2] += fpairz;
        if (newton_pair || j < nlocal) {
          f[j][0] -= fpairx;
          f[j][1] -= fpairy;
          f[j][2] -= fpairz;
        }

        if (eflag) {
          // unshifted eng of conservative term:
          // evdwl = -a0[itype][jtype]*r * (1.0-0.5*r/cut[itype][jtype]);
          // eng shifted to 0.0 at cutoff
          evdwl = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
          evdwl *= factor_dpd;
        }
        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,evdwl,0.0,fpairx,fpairy,fpairz,
                                 delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairDPDExt::allocate()
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
  memory->create(a0,n+1,n+1,"pair:a0");
  memory->create(gamma,n+1,n+1,"pair:gamma");
  memory->create(gammaT,n+1,n+1,"pair:gammaT");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(sigmaT,n+1,n+1,"pair:sigmaT");
  memory->create(ws,n+1,n+1,"pair:ws");
  memory->create(wsT,n+1,n+1,"pair:wsT");
  for (i = 0; i <= atom->ntypes; i++)
  {
    for (j = 0; j <= atom->ntypes; j++)
    {
      sigma[i][j] = gamma[i][j] =sigmaT[i][j] = gammaT[i][j] = 0.0;
      ws[i][j] = wsT[i][j] = 1.0;
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDPDExt::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  temperature = utils::numeric(FLERR,arg[0],false,lmp);
  cut_global = utils::numeric(FLERR,arg[1],false,lmp);
  seed = utils::inumeric(FLERR,arg[2],false,lmp);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanMars(lmp,seed + comm->me);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut[i][j] = cut_global;
          cutsq[i][j] = cut_global*cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDExt::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double gamma_one = utils::numeric(FLERR,arg[3],false,lmp);
  double gammaT_one = utils::numeric(FLERR,arg[4],false,lmp);
  double ws_one = utils::numeric(FLERR,arg[5],false,lmp);
  double wsT_one = utils::numeric(FLERR,arg[6],false,lmp);

  double cut_one = cut_global;
  if (narg == 8) cut_one = utils::numeric(FLERR,arg[7],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      gamma[i][j] = gamma_one;
      gammaT[i][j] = gammaT_one;
      ws[i][j] = ws_one;
      wsT[i][j] = wsT_one;
      cut[i][j] = cut_one;
      cutsq[i][j] = cut_one*cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDExt::init_style()
{
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair dpd requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0)
    error->warning(FLERR, "Pair dpd needs newton pair on for momentum conservation");

  neighbor->add_request(this);

  // precompute random force scaling factors

  for (int i = 0; i < 4; ++i) special_sqrt[i] = sqrt(force->special_lj[i]);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDPDExt::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  sigma[i][j] = sqrt(2.0*force->boltz*temperature*gamma[i][j]);
  sigmaT[i][j] = sqrt(2.0*force->boltz*temperature*gammaT[i][j]);

  cut[j][i] = cut[i][j];
  cutsq[j][i] = cutsq[i][j];
  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  gammaT[j][i] = gammaT[i][j];
  sigma[j][i] = sigma[i][j];
  sigmaT[j][i] = sigmaT[i][j];
  ws[j][i] = ws[i][j];
  wsT[j][i] = wsT[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDExt::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a0[i][j],sizeof(double),1,fp);
        fwrite(&gamma[i][j],sizeof(double),1,fp);
        fwrite(&gammaT[i][j],sizeof(double),1,fp);
        fwrite(&ws[i][j],sizeof(double),1,fp);
        fwrite(&wsT[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDExt::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&a0[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gammaT[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&ws[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&wsT[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&gammaT[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&ws[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&wsT[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDExt::write_restart_settings(FILE *fp)
{
  fwrite(&temperature,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDExt::read_restart_settings(FILE *fp)
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

void PairDPDExt::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,
      a0[i][i],gamma[i][i],gammaT[i][i],ws[i][i],wsT[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDPDExt::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g\n",i,j,
        a0[i][j],gamma[i][j],gammaT[i][j],ws[i][j],wsT[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairDPDExt::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                       double /*factor_coul*/, double factor_dpd, double &fforce)
{
  double r,rinv,wd,phi;

  r = sqrt(rsq);
  if (r < EPSILON) {
    fforce = 0.0;
    return 0.0;
  }

  rinv = 1.0/r;
  wd = 1.0 - r/cut[itype][jtype];
  fforce = a0[itype][jtype]*wd * factor_dpd*rinv;

  phi = 0.5*a0[itype][jtype]*cut[itype][jtype] * wd*wd;
  return factor_dpd*phi;
}
