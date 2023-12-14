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

#include "pair_spin_exchange_biquadratic.h"

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

PairSpinExchangeBiquadratic::PairSpinExchangeBiquadratic(LAMMPS *lmp) :
  PairSpin(lmp)
{
  e_offset = 0;
}

/* ---------------------------------------------------------------------- */

PairSpinExchangeBiquadratic::~PairSpinExchangeBiquadratic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cut_spin_exchange);
    memory->destroy(J1_mag);
    memory->destroy(J1_mech);
    memory->destroy(J2);
    memory->destroy(J3);
    memory->destroy(K1_mag);
    memory->destroy(K1_mech);
    memory->destroy(K2);
    memory->destroy(K3);
    memory->destroy(cutsq); // to be implemented
    memory->destroy(emag);
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::settings(int narg, char **arg)
{
  PairSpin::settings(narg,arg);

  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_spin_exchange_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_spin_exchange[i][j] = cut_spin_exchange_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type spin pairs
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // check if args correct

  if (strcmp(arg[2],"biquadratic") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if ((narg != 10) && (narg != 12))
    error->all(FLERR,"Incorrect args for pair coefficients");

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // get exchange arguments from input command

  int iarg = 10;
  const double rc = utils::numeric(FLERR,arg[3],false,lmp);
  const double j1 = utils::numeric(FLERR,arg[4],false,lmp);
  const double j2 = utils::numeric(FLERR,arg[5],false,lmp);
  const double j3 = utils::numeric(FLERR,arg[6],false,lmp);
  const double k1 = utils::numeric(FLERR,arg[7],false,lmp);
  const double k2 = utils::numeric(FLERR,arg[8],false,lmp);
  const double k3 = utils::numeric(FLERR,arg[9],false,lmp);

  // read energy offset flag if specified

  while (iarg < narg) {
    if (strcmp(arg[iarg],"offset") == 0) {
      e_offset = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_spin_exchange[i][j] = rc;
      J1_mag[i][j] = j1/hbar;
      J1_mech[i][j] = j1;
      J2[i][j] = j2;
      J3[i][j] = j3;
      K1_mag[i][j] = k1/hbar;
      K1_mech[i][j] = k1;
      K2[i][j] = k2;
      K3[i][j] = k3;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args in pair_style command");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSpinExchangeBiquadratic::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  J1_mag[j][i] = J1_mag[i][j];
  J1_mech[j][i] = J1_mech[i][j];
  J2[j][i] = J2[i][j];
  J3[j][i] = J3[i][j];
  K1_mag[j][i] = K1_mag[i][j];
  K1_mech[j][i] = K1_mech[i][j];
  K2[j][i] = K2[i][j];
  K3[j][i] = K3[i][j];
  cut_spin_exchange[j][i] = cut_spin_exchange[i][j];

  return cut_spin_exchange_global;
}

/* ----------------------------------------------------------------------
   extract the larger cutoff
------------------------------------------------------------------------- */

void *PairSpinExchangeBiquadratic::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut") == 0) return (void *) &cut_spin_exchange_global;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::compute(int eflag, int vflag)
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

  if (nlocal_max < nlocal) {    // grow emag lists if necessary
    nlocal_max = nlocal;
    memory->grow(emag,nlocal_max,"pair/spin:emag");
  }

  // computation of the exchange interaction
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

      local_cut2 = cut_spin_exchange[itype][jtype]*cut_spin_exchange[itype][jtype];

      // compute exchange interaction

      if (rsq <= local_cut2) {
        compute_exchange(i,j,rsq,fmi,spi,spj);

        if (lattice_flag)
          compute_exchange_mech(i,j,rsq,eij,fi,spi,spj);

        if (eflag) {
          evdwl -= compute_energy(i,j,rsq,spi,spj);
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

void PairSpinExchangeBiquadratic::compute_single_pair(int ii, double fmi[3])
{
  int *type = atom->type;
  double **x = atom->x;
  double **sp = atom->sp;
  double local_cut2;
  double xi[3];
  double delx,dely,delz;
  double spi[3],spj[3];

  int j,jnum,itype,jtype,ntypes;
  int k,locflag;
  int *jlist,*numneigh,**firstneigh;

  double rsq;

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
    spi[0] = sp[ii][0];
    spi[1] = sp[ii][1];
    spi[2] = sp[ii][2];

    jlist = firstneigh[ii];
    jnum = numneigh[ii];

    for (int jj = 0; jj < jnum; jj++) {

      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      local_cut2 = cut_spin_exchange[itype][jtype]*cut_spin_exchange[itype][jtype];

      spj[0] = sp[j][0];
      spj[1] = sp[j][1];
      spj[2] = sp[j][2];

      delx = xi[0] - x[j][0];
      dely = xi[1] - x[j][1];
      delz = xi[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= local_cut2) {
        compute_exchange(ii,j,rsq,fmi,spi,spj);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute exchange interaction between spins i and j
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::compute_exchange(int i, int j, double rsq,
    double fmi[3], double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype,jtype;
  double Jex,Kex,r2j,r2k,sdots;
  itype = type[i];
  jtype = type[j];

  r2j = rsq/J3[itype][jtype]/J3[itype][jtype];
  r2k = rsq/K3[itype][jtype]/K3[itype][jtype];

  Jex = 4.0*J1_mag[itype][jtype]*r2j;
  Jex *= (1.0-J2[itype][jtype]*r2j);
  Jex *= exp(-r2j);

  Kex = 4.0*K1_mag[itype][jtype]*r2k;
  Kex *= (1.0-K2[itype][jtype]*r2k);
  Kex *= exp(-r2k);

  sdots = (spi[0]*spj[0]+spi[1]*spj[1]+spi[2]*spj[2]);

  fmi[0] += (Jex*spj[0] + 2.0*Kex*spj[0]*sdots);
  fmi[1] += (Jex*spj[1] + 2.0*Kex*spj[1]*sdots);
  fmi[2] += (Jex*spj[2] + 2.0*Kex*spj[2]*sdots);
}

/* ----------------------------------------------------------------------
   compute the mechanical force due to the exchange interaction between atom i and atom j
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::compute_exchange_mech(int i, int j,
    double rsq, double eij[3], double fi[3],  double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype,jtype;
  double Jex,Jex_mech,Kex,Kex_mech,sdots;
  double rja,rka,rjr,rkr,iJ3,iK3;
  double fx, fy, fz;
  itype = type[i];
  jtype = type[j];

  Jex = J1_mech[itype][jtype];
  iJ3 = 1.0/(J3[itype][jtype]*J3[itype][jtype]);
  Kex = K1_mech[itype][jtype];
  iK3 = 1.0/(K3[itype][jtype]*K3[itype][jtype]);

  rja = rsq*iJ3;
  rjr = sqrt(rsq)*iJ3;
  rka = rsq*iK3;
  rkr = sqrt(rsq)*iK3;

  Jex_mech = 1.0-rja-J2[itype][jtype]*rja*(2.0-rja);
  Jex_mech *= 8.0*Jex*rjr*exp(-rja);

  Kex_mech = 1.0-rka-K2[itype][jtype]*rka*(2.0-rka);
  Kex_mech *= 8.0*Kex*rkr*exp(-rka);

  sdots = (spi[0]*spj[0]+spi[1]*spj[1]+spi[2]*spj[2]);

  // apply or not energy and force offset

  fx = fy = fz = 0.0;
  if (e_offset == 1) { // set offset
    fx = (Jex_mech*(sdots-1.0) + Kex_mech*(sdots*sdots-1.0))*eij[0];
    fy = (Jex_mech*(sdots-1.0) + Kex_mech*(sdots*sdots-1.0))*eij[1];
    fz = (Jex_mech*(sdots-1.0) + Kex_mech*(sdots*sdots-1.0))*eij[2];
  } else if (e_offset == 0) { // no offset ("normal" calculation)
    fx =  (Jex_mech*sdots + Kex_mech*sdots*sdots)*eij[0];
    fy =  (Jex_mech*sdots + Kex_mech*sdots*sdots)*eij[1];
    fz =  (Jex_mech*sdots + Kex_mech*sdots*sdots)*eij[2];
  } else error->all(FLERR,"Illegal option in pair exchange/biquadratic command");

  fi[0] -= 0.5*fx;
  fi[1] -= 0.5*fy;
  fi[2] -= 0.5*fz;
  // fi[0] -= fx;
  // fi[1] -= fy;
  // fi[2] -= fz;
}

/* ----------------------------------------------------------------------
   compute energy of spin pair i and j
------------------------------------------------------------------------- */

double PairSpinExchangeBiquadratic::compute_energy(int i, int j, double rsq,
    double spi[3], double spj[3])
{
  int *type = atom->type;
  int itype,jtype;
  double Jex,Kex,sdots;
  double r2j,r2k;
  double energy = 0.0;
  itype = type[i];
  jtype = type[j];

  r2j = rsq/J3[itype][jtype]/J3[itype][jtype];
  r2k = rsq/K3[itype][jtype]/K3[itype][jtype];

  Jex = 4.0*J1_mech[itype][jtype]*r2j;
  Jex *= (1.0-J2[itype][jtype]*r2j);
  Jex *= exp(-r2j);

  Kex = 4.0*K1_mech[itype][jtype]*r2k;
  Kex *= (1.0-K2[itype][jtype]*r2k);
  Kex *= exp(-r2k);

  sdots = (spi[0]*spj[0]+spi[1]*spj[1]+spi[2]*spj[2]);

  // apply or not energy and force offset

  if (e_offset == 1) { // set offset
    energy = 0.5*(Jex*(sdots-1.0) + Kex*(sdots*sdots-1.0));
  } else if (e_offset == 0) { // no offset ("normal" calculation)
    energy = 0.5*(Jex*sdots + Kex*sdots*sdots);
  } else error->all(FLERR,"Illegal option in pair exchange/biquadratic command");

  return energy;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cut_spin_exchange,n+1,n+1,"pair/spin/exchange:cut_spin_exchange");
  memory->create(J1_mag,n+1,n+1,"pair/spin/exchange:J1_mag");
  memory->create(J1_mech,n+1,n+1,"pair/spin/exchange:J1_mech");
  memory->create(J2,n+1,n+1,"pair/spin/exchange:J2");
  memory->create(J3,n+1,n+1,"pair/spin/exchange:J3");
  memory->create(K1_mag,n+1,n+1,"pair/spin/exchange:J1_mag");
  memory->create(K1_mech,n+1,n+1,"pair/spin/exchange:J1_mech");
  memory->create(K2,n+1,n+1,"pair/spin/exchange:J2");
  memory->create(K3,n+1,n+1,"pair/spin/exchange:J3");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&J1_mag[i][j],sizeof(double),1,fp);
        fwrite(&J1_mech[i][j],sizeof(double),1,fp);
        fwrite(&J2[i][j],sizeof(double),1,fp);
        fwrite(&J3[i][j],sizeof(double),1,fp);
        fwrite(&K1_mag[i][j],sizeof(double),1,fp);
        fwrite(&K1_mech[i][j],sizeof(double),1,fp);
        fwrite(&K2[i][j],sizeof(double),1,fp);
        fwrite(&K3[i][j],sizeof(double),1,fp);
        fwrite(&cut_spin_exchange[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&J1_mag[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&J1_mech[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&J2[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&J3[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&K1_mag[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&K1_mech[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&K2[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&K3[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_spin_exchange[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&J1_mag[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&J1_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&J2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&J3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K1_mag[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K1_mech[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K2[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&K3[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_spin_exchange[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }
}


/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::write_restart_settings(FILE *fp)
{
  fwrite(&cut_spin_exchange_global,sizeof(double),1,fp);
  fwrite(&e_offset,sizeof(int),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairSpinExchangeBiquadratic::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_spin_exchange_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&e_offset,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_spin_exchange_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&e_offset,1,MPI_INT,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}
