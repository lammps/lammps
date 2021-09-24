/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "bond_fene_nm_split.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include <cmath>
#include <cstring>



using namespace LAMMPS_NS;
using MathConst::MY_CUBEROOT2;


/* ---------------------------------------------------------------------- */

BondFENEnm::~BondFENEnm()
{
if (allocated && !copymode) {
       memory->destroy(setflag);
       memory->destroy(k);
       memory->destroy(r0);
       memory->destroy(epsilon);
       memory->destroy(sigma);
       memory->destroy(nn);
       memory->destroy(mm);
    }
}

/* ---------------------------------------------------------------------- */

void BondFENEnm::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r0sq,rlogarg,sr2,sr6;
  double r;

  ebond = sr6 = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    // force from log term

    rsq = delx*delx + dely*dely + delz*delz;
    r0sq = r0[type] * r0[type];
    rlogarg = 1.0 - rsq/r0sq;

      // if r -> r0, then rlogarg < 0.0 which is an error
      // issue a warning and reset rlogarg = epsilon
      // if r > 2*r0 something serious is wrong, abort

      // change cutuff from .1 to .02 so only bond lengths > 1.485 give the warning
      // and crash the run if rlogarg < -.21 rather than < 3
      // Don't print out warnings, only errors
    if (rlogarg < .02) {
      char str[128];
      sprintf(str,"FENE bond too long: " BIGINT_FORMAT " "
              TAGINT_FORMAT " " TAGINT_FORMAT " %g",
              update->ntimestep,atom->tag[i1],atom->tag[i2],sqrt(rsq));
      if (rlogarg <= -.21) error->one(FLERR,"Bad FENE bond");
      rlogarg = 0.02;
    }

    fbond = -k[type]/rlogarg;

      // force from n-m term
      // MY_CUBEROOT2 cutoff assumes sigma = 2^{1/6}
    if (rsq < MY_CUBEROOT2) {
        r = sqrt(rsq);
        fbond += epsilon[type]*(nn[type]*mm[type]/(nn[type]-mm[type]))*(pow(sigma[type]/r,nn[type])-pow(sigma[type]/r,mm[type]))/rsq;
    }

    // energy

    if (eflag) {
      ebond = -0.5 * k[type]*r0sq*log(rlogarg);
      if (rsq < MY_CUBEROOT2) ebond += (epsilon[type]/(nn[type]-mm[type]))*(mm[type]*pow(sigma[type]/r,nn[type])-nn[type]*pow(sigma[type]/r,mm[type]));
    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond;
      f[i1][1] += dely*fbond;
      f[i1][2] += delz*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond;
      f[i2][1] -= dely*fbond;
      f[i2][2] -= delz*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondFENEnm::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");
  memory->create(epsilon,n+1,"bond:epsilon");
  memory->create(sigma,n+1,"bond:sigma");
  memory->create(nn,n+1,"bond:nn");
  memory->create(mm,n+1,"bond:mm");
    memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondFENEnm::coeff(int narg, char **arg)
{
  if (narg != 7) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double k_one = utils::numeric(FLERR,arg[1],false,lmp);
  double r0_one = utils::numeric(FLERR,arg[2],false,lmp);
  double epsilon_one = utils::numeric(FLERR,arg[3],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[4],false,lmp);
    double nn_one = utils::numeric(FLERR,arg[5],false,lmp);
    double mm_one = utils::numeric(FLERR,arg[6],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    epsilon[i] = epsilon_one;
    sigma[i] = sigma_one;
      nn[i] = nn_one;
      mm[i] = mm_one;

    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if special_bond settings are valid
------------------------------------------------------------------------- */

void BondFENEnm::init_style()
{
  // special bonds should be 0 1 1

  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0) {
    if (comm->me == 0)
      error->warning(FLERR,"Use special bonds = 0,1,1 with bond style fene");
  }
}

/* ---------------------------------------------------------------------- */

double BondFENEnm::equilibrium_distance(int i)
{
  return 0.97*sigma[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondFENEnm::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&epsilon[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&sigma[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&nn[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&mm[1],sizeof(double),atom->nbondtypes,fp);

}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondFENEnm::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&r0[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&epsilon[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&sigma[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&nn[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&mm[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);

  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsilon[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma[1],atom->nbondtypes,MPI_DOUBLE,0,world);
    MPI_Bcast(&nn[1],atom->nbondtypes,MPI_DOUBLE,0,world);
    MPI_Bcast(&mm[1],atom->nbondtypes,MPI_DOUBLE,0,world);


  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondFENEnm::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g\n",i,k[i],r0[i],epsilon[i],sigma[i],nn[i],mm[i]);
}

/* ---------------------------------------------------------------------- */

double BondFENEnm::single(int type, double rsq, int /*i*/, int /*j*/,
                        double &fforce)
{
  double r0sq = r0[type] * r0[type];
  double rlogarg = 1.0 - rsq/r0sq;
  double r;
    // if r -> r0, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*r0 something serious is wrong, abort

    // change cutuff from .1 to .02 so only bond lengths > 1.485 give the warning
    // and crash the run if rlogarg < -.21 rather than < 3
    // Don't print out warnings, only errors
    if (rlogarg < 0.02) {
      char str[128];
      sprintf(str,"FENE bond too long: " BIGINT_FORMAT " %g",
              update->ntimestep,sqrt(rsq));
        
      if (rlogarg <= -.21) error->one(FLERR,"Bad FENE bond");
      rlogarg = 0.02;;
  }

  double eng = -0.5 * k[type]*r0sq*log(rlogarg);
  fforce = -k[type]/rlogarg;
    
    // MY_CUBEROOT2 cutoff assumes sigma = 2^1/6
    if (rsq < sigma[type]*sigma[type]) {
       r = sqrt(rsq);
       fforce += epsilon[type]*(nn[type]*mm[type]/(nn[type]-mm[type]))*(pow(sigma[type]/r,nn[type])-pow(sigma[type]/r,mm[type]))/rsq;
       eng += (epsilon[type]/(nn[type]-mm[type]))*(mm[type]*pow(sigma[type]/r,nn[type])-nn[type]*pow(sigma[type]/r,mm[type]));
    }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *BondFENEnm::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str,"kappa")==0) return (void*) k;
  if (strcmp(str,"r0")==0) return (void*) r0;
  return nullptr;
}
