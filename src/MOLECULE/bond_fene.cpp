/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_fene.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondFENE::BondFENE(LAMMPS *lmp) : Bond(lmp)
{
  TWO_1_3 = pow(2.0,(1.0/3.0));
}

/* ---------------------------------------------------------------------- */

BondFENE::~BondFENE()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(k);
    memory->sfree(r0);
    memory->sfree(epsilon);
    memory->sfree(sigma);
  }
}

/* ---------------------------------------------------------------------- */

void BondFENE::compute(int eflag, int vflag)
{
  int i1,i2,n,type,factor;
  double delx,dely,delz,rsq,r0sq,rlogarg,fforce,sr2,sr6,rfactor;

  energy = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

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

    if (newton_bond) factor = 2;
    else {
      factor = 0;
      if (i1 < nlocal) factor++;
      if (i2 < nlocal) factor++;
    }
    rfactor = 0.5*factor;

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);

    // force from log term

    rsq = delx*delx + dely*dely + delz*delz;
    r0sq = r0[type] * r0[type];
    rlogarg = 1.0 - rsq/r0sq;

    // if r -> r0, then rlogarg < 0.0 which is an error
    // issue a warning and reset rlogarg = epsilon
    // if r > 2*r0 something serious is wrong, abort

    if (rlogarg < 0.1) {
      char str[128];
      sprintf(str,"FENE bond too long: %d %d %d %g",
              update->ntimestep,atom->tag[i1],atom->tag[i2],sqrt(rsq));
      error->warning(str);
      if (rlogarg <= -3.0) error->one("Bad FENE bond");
      rlogarg = 0.1;
    }

    fforce = -k[type]/rlogarg;

    // force from LJ term

    if (rsq < TWO_1_3*sigma[type]*sigma[type]) {
      sr2 = sigma[type]*sigma[type]/rsq;
      sr6 = sr2*sr2*sr2;
      fforce += 48.0*epsilon[type]*sr6*(sr6-0.5)/rsq;
    }

    // energy

    if (eflag) {
      energy -= 0.5*rfactor * k[type]*r0sq*log(rlogarg);
      if (rsq < TWO_1_3*sigma[type]*sigma[type])
	energy += rfactor * (4.0*epsilon[type]*sr6*(sr6-1.0) + epsilon[type]);
    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fforce;
      f[i1][1] += dely*fforce;
      f[i1][2] += delz*fforce;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fforce;
      f[i2][1] -= dely*fforce;
      f[i2][2] -= delz*fforce;
    }

    // virial contribution

    if (vflag) {
      virial[0] += rfactor * delx*delx*fforce;
      virial[1] += rfactor * dely*dely*fforce;
      virial[2] += rfactor * delz*delz*fforce;
      virial[3] += rfactor * delx*dely*fforce;
      virial[4] += rfactor * delx*delz*fforce;
      virial[5] += rfactor * dely*delz*fforce;
    }
  }
}

/* ---------------------------------------------------------------------- */

void BondFENE::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  k = (double *) memory->smalloc((n+1)*sizeof(double),"bond:k");
  r0 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:r0");
  epsilon = (double *) memory->smalloc((n+1)*sizeof(double),"bond:epsilon");
  sigma = (double *) memory->smalloc((n+1)*sizeof(double),"bond:sigma");
  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondFENE::coeff(int narg, char **arg)
{
  if (narg != 5) error->all("Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = atof(arg[1]);
  double r0_one = atof(arg[2]);
  double epsilon_one = atof(arg[3]);
  double sigma_one = atof(arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    epsilon[i] = epsilon_one;
    sigma[i] = sigma_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args for bond coefficients");
}

/* ---------------------------------------------------------------------- */

double BondFENE::equilibrium_distance(int i)
{
  return 0.97*sigma[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondFENE::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&epsilon[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&sigma[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondFENE::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&epsilon[1],sizeof(double),atom->nbondtypes,fp);
    fread(&sigma[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsilon[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

void BondFENE::single(int type, double rsq, int i, int j,
		      int eflag, double &fforce, double &eng)
{
  double r0sq = r0[type] * r0[type];
  double rlogarg = 1.0 - rsq/r0sq;

  // if r -> r0, then rlogarg < 0.0 which is an error
  // issue a warning and reset rlogarg = epsilon
  // if r > 2*r0 something serious is wrong, abort

  if (rlogarg < 0.1) {
    char str[128];
    sprintf(str,"FENE bond too long: %d %g",update->ntimestep,sqrt(rsq));
    error->warning(str);
    if (rlogarg <= -3.0) error->one("Bad FENE bond");
    rlogarg = 0.1;
  }

  fforce = -k[type]/rlogarg;

  // force from LJ term

  double sr2,sr6;
  if (rsq < TWO_1_3*sigma[type]*sigma[type]) {
    sr2 = sigma[type]*sigma[type]/rsq;
    sr6 = sr2*sr2*sr2;
    fforce += 48.0*epsilon[type]*sr6*(sr6-0.5)/rsq;
  }

  // energy

  if (eflag) {
    eng = -0.5 * k[type]*r0sq*log(rlogarg);
    if (rsq < TWO_1_3*sigma[type]*sigma[type])
      eng += 4.0*epsilon[type]*sr6*(sr6-1.0) + epsilon[type];
  }
}
