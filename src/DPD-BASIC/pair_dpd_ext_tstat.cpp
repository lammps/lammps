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
------------------------------------------------------------------------- */

#include "pair_dpd_ext_tstat.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr double EPSILON = 1.0e-10;

/* ---------------------------------------------------------------------- */

PairDPDExtTstat::PairDPDExtTstat(LAMMPS *lmp) : PairDPDExt(lmp)
{
  single_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

void PairDPDExtTstat::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpairx,fpairy,fpairz,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,wdPar,wdPerp,randnum,randnumx,randnumy,randnumz;
  double prefactor_g,prefactor_s,factor_dpd,factor_sqrt;
  double P[3][3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  ev_init(eflag,vflag);

  // adjust sigma if target T is changing

  if (t_start != t_stop) {
    double delta = update->ntimestep - update->beginstep;
    if (delta != 0.0) delta /= update->endstep - update->beginstep;
    temperature = t_start + delta * (t_stop-t_start);
    double boltz = force->boltz;
    for (i = 1; i <= atom->ntypes; i++) {
      for (j = i; j <= atom->ntypes; j++) {
        sigma[i][j] = sigma[j][i] = sqrt(2.0*boltz*temperature*gamma[i][j]);
        sigmaT[i][j] = sigmaT[j][i] = sqrt(2.0*boltz*temperature*gammaT[i][j]);
      }
    }
  }

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

        // drag force - parallel
        fpair = -factor_dpd * gamma[itype][jtype]*wdPar*wdPar*dot*rinv;

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

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                      0.0,0.0,
                      fpairx, fpairy, fpairz,
                      delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairDPDExtTstat::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  t_start = utils::numeric(FLERR,arg[0],false,lmp);
  t_stop = utils::numeric(FLERR,arg[1],false,lmp);
  cut_global = utils::numeric(FLERR,arg[2],false,lmp);
  seed = utils::inumeric(FLERR,arg[3],false,lmp);

  temperature = t_start;

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all(FLERR,"Illegal pair_style command");
  delete random;
  random = new RanMars(lmp,seed + comm->me);

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

void PairDPDExtTstat::coeff(int narg, char **arg)
{
  if (narg < 6 || narg > 7)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double a0_one = 0.0;
  double gamma_one = utils::numeric(FLERR,arg[2],false,lmp);
  double gammaT_one = utils::numeric(FLERR,arg[3],false,lmp);
  double ws_one = utils::numeric(FLERR,arg[4],false,lmp);
  double wsT_one = utils::numeric(FLERR,arg[5],false,lmp);

  double cut_one = cut_global;
  if (narg == 7) cut_one = utils::numeric(FLERR,arg[6],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a0[i][j] = a0_one;
      gamma[i][j] = gamma_one;
      gammaT[i][j] = gammaT_one;
      ws[i][j] = ws_one;
      wsT[i][j] = wsT_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDExtTstat::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
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

void PairDPDExtTstat::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&gamma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&gammaT[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&ws[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&wsT[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
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

void PairDPDExtTstat::write_restart_settings(FILE *fp)
{
  fwrite(&t_start,sizeof(double),1,fp);
  fwrite(&t_stop,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDExtTstat::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&t_start,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&t_stop,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&seed,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&t_start,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&t_stop,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  temperature = t_start;

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairDPDExtTstat::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,gamma[i][i],gammaT[i][i],ws[i][i],wsT[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairDPDExtTstat::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,
        gamma[i][j],gammaT[i][j],ws[i][j],wsT[i][j],cut[i][j]);
}
