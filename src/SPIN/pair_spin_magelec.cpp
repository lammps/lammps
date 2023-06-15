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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "pair_spin_magelec.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSpinMagelec::~PairSpinMagelec()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_magelec);
    memory->destroy(ME);
    memory->destroy(ME_mech);
    memory->destroy(v_mex);
    memory->destroy(v_mey);
    memory->destroy(v_mez);
    memory->destroy(cutsq); // to be deteled
    memory->destroy(emag);
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinMagelec::settings(int narg, char **arg)
{

  PairSpin::settings(narg,arg);

  cut_spin_magelec_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_spin_magelec[i][j] = cut_spin_magelec_global;
        }
  }

}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs (only one for now)
------------------------------------------------------------------------- */

void PairSpinMagelec::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // check if args correct

  if (strcmp(arg[2],"magelec") != 0)
    error->all(FLERR,"Incorrect args in pair_style command");
  if (narg != 8)
    error->all(FLERR,"Incorrect args in pair_style command");

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  const double rij = utils::numeric(FLERR,arg[3],false,lmp);
  const double magelec = utils::numeric(FLERR,arg[4],false,lmp);
  double mex = utils::numeric(FLERR,arg[5],false,lmp);
  double mey = utils::numeric(FLERR,arg[6],false,lmp);
  double mez = utils::numeric(FLERR,arg[7],false,lmp);

  double inorm = 1.0/(mex*mex+mey*mey+mez*mez);
  mex *= inorm;
  mey *= inorm;
  mez *= inorm;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_spin_magelec[i][j] = rij;
      ME[i][j] = magelec/hbar;
      ME_mech[i][j] = magelec;
      v_mex[i][j] = mex;
      v_mey[i][j] = mey;
      v_mez[i][j] = mez;
      setflag[i][j] = 1;
      count++;
    }
  }
  if (count == 0)
    error->all(FLERR,"Incorrect args in pair_style command");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinMagelec::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  ME[j][i] = ME[i][j];
  ME_mech[j][i] = ME_mech[i][j];
  v_mex[j][i] = v_mex[i][j];
  v_mey[j][i] = v_mey[i][j];
  v_mez[j][i] = v_mez[i][j];
  cut_spin_magelec[j][i] = cut_spin_magelec[i][j];

  return cut_spin_magelec_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinMagelec::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut") == 0) return (void *) &cut_spin_magelec_global;
  return nullptr;
}


/* ---------------------------------------------------------------------- */

void PairSpinMagelec::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl, ecoul;
  double xi[3], eij[3];
  double delx,dely,delz;
  double spi[3], spj[3];
  double fi[3], fmi[3];
  double local_cut2;
  double rsq, inorm;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **fm = atom->fm;
  double **sp = atom->sp;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // checking size of emag

  if (nlocal_max < nlocal) {                    // grow emag lists if necessary
    nlocal_max = nlocal;
    memory->grow(emag,nlocal_max,"pair/spin:emag");
  }

  // magneto-electric computation
  // loop over atoms and their neighbors

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    spi[0] = sp[i][0];
    spi[1] = sp[i][1];
    spi[2] = sp[i][2];
    emag[i] = 0.0;

    // loop on neighbors

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      evdwl = 0.0;
      fi[0] = fi[1] = fi[2] = 0.0;
      fmi[0] = fmi[1] = fmi[2] = 0.0;

      delx = xi[0] - x[j][0];
      dely = xi[1] - x[j][1];
      delz = xi[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      inorm = 1.0/sqrt(rsq);
      eij[0] = -inorm*delx;
      eij[1] = -inorm*dely;
      eij[2] = -inorm*delz;

      local_cut2 = cut_spin_magelec[itype][jtype]*cut_spin_magelec[itype][jtype];

      // compute me interaction

      if (rsq <= local_cut2) {
        compute_magelec(i,j,eij,fmi,spj);

        if (lattice_flag)
          compute_magelec_mech(i,j,fi,spi,spj);

        if (eflag) {
          evdwl -= (spi[0]*fmi[0] + spi[1]*fmi[1] + spi[2]*fmi[2]);
          evdwl *= 0.5*hbar;
          emag[i] += evdwl;
        } else evdwl = 0.0;

        f[i][0] += fi[0];
        f[i][1] += fi[1];
        f[i][2] += fi[2];
        if (newton_pair || j < nlocal) {
          f[j][0] -= fi[0];
          f[j][1] -= fi[1];
          f[j][2] -= fi[2];
        }
        fm[i][0] += fmi[0];
        fm[i][1] += fmi[1];
        fm[i][2] += fmi[2];

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
            evdwl,ecoul,fi[0],fi[1],fi[2],delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   update the pair interactions fmi acting on the spin ii
------------------------------------------------------------------------- */

void PairSpinMagelec::compute_single_pair(int ii, double fmi[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2;
  double xi[3], eij[3];
  double delx,dely,delz;
  double spj[3];

  int j,jnum,itype,jtype,ntypes;
  int k,locflag;
  int *jlist,*numneigh,**firstneigh;

  double rsq, inorm;

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // check if interaction applies to type of ii

  itype = type[ii];
  ntypes = atom->ntypes;
  locflag = 0;
  k = 1;
  while (k <= ntypes) {
    if (k <= itype) {
      if (setflag[k][itype] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else if (k > itype) {
      if (setflag[itype][k] == 1) {
        locflag =1;
        break;
      }
      k++;
    } else error->all(FLERR,"Wrong type number");
  }

  // if interaction applies to type ii,
  // locflag = 1 and compute pair interaction

  if (locflag == 1) {

    xi[0] = x[ii][0];
    xi[1] = x[ii][1];
    xi[2] = x[ii][2];

    jlist = firstneigh[ii];
    jnum = numneigh[ii];

    for (int jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      local_cut2 = cut_spin_magelec[itype][jtype]*cut_spin_magelec[itype][jtype];

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      delx = xi[0] - x[j][0];
      dely = xi[1] - x[j][1];
      delz = xi[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      inorm = 1.0/sqrt(rsq);
      eij[0] = -inorm*delx;
      eij[1] = -inorm*dely;
      eij[2] = -inorm*delz;

      if (rsq <= local_cut2) {
        compute_magelec(ii,j,eij,fmi,spj);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSpinMagelec::compute_magelec(int i, int j, double eij[3], double fmi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  double meix,meiy,meiz;
  double vx,vy,vz;
  itype = type[i];
  jtype = type[j];

  vx = v_mex[itype][jtype];
  vy = v_mey[itype][jtype];
  vz = v_mez[itype][jtype];

  meix = (vy*eij[2] - vz*eij[1]);
  meiy = (vz*eij[0] - vx*eij[2]);
  meiz = (vx*eij[1] - vy*eij[0]);

  meix *= ME[itype][jtype];
  meiy *= ME[itype][jtype];
  meiz *= ME[itype][jtype];

  fmi[0] += (spj[1]*meiz - spj[2]*meiy);
  fmi[1] += (spj[2]*meix - spj[0]*meiz);
  fmi[2] += (spj[0]*meiy - spj[1]*meix);
}

/* ---------------------------------------------------------------------- */

void PairSpinMagelec::compute_magelec_mech(int i, int j, double fi[3], double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype, jtype;
  itype = type[i];
  jtype = type[j];

  double meix,meiy,meiz;
  double vx, vy, vz;

  vx = v_mex[itype][jtype];
  vy = v_mey[itype][jtype];
  vz = v_mez[itype][jtype];

  meix = (spi[1]*spi[2] - spi[2]*spj[1]);
  meiy = (spi[2]*spi[0] - spi[0]*spj[2]);
  meiz = (spi[0]*spi[1] - spi[1]*spj[0]);

  meix *= ME_mech[itype][jtype];
  meiy *= ME_mech[itype][jtype];
  meiz *= ME_mech[itype][jtype];

  fi[0] += 0.5*(meiy*vz - meiz*vy);
  fi[1] += 0.5*(meiz*vx - meix*vz);
  fi[2] += 0.5*(meix*vy - meiy*vx);

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinMagelec::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cut_spin_magelec,n+1,n+1,"pair/spin/me:cut_spin_magelec");
  memory->create(ME,n+1,n+1,"pair/spin/me:ME");
  memory->create(ME_mech,n+1,n+1,"pair/spin/me:ME_mech");
  memory->create(v_mex,n+1,n+1,"pair/spin/me:ME_vector_x");
  memory->create(v_mey,n+1,n+1,"pair/spin/me:ME_vector_y");
  memory->create(v_mez,n+1,n+1,"pair/spin/me:ME_vector_z");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinMagelec::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&ME[i][j],sizeof(double),1,fp);
        fwrite(&v_mex[i][j],sizeof(double),1,fp);
        fwrite(&v_mey[i][j],sizeof(double),1,fp);
        fwrite(&v_mez[i][j],sizeof(double),1,fp);
        fwrite(&cut_spin_magelec[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinMagelec::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&ME[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&v_mex[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&v_mey[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&v_mez[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_spin_magelec[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&ME[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_mex[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_mey[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&v_mez[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spin_magelec[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinMagelec::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_magelec_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinMagelec::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_spin_magelec_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_spin_magelec_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
