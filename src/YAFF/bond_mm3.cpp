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
   Contributing author: Steven Vandenbrande
------------------------------------------------------------------------- */

#include "bond_mm3.h"

#include <cmath>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondMM3::BondMM3(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondMM3::~BondMM3()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(k2);
  }
}

/* ---------------------------------------------------------------------- */

void BondMM3::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,dr2,de_bond,K3,K4;

  ebond = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  /*
  E = K(r-r0)^2 [1-2.55*(r-r0)+(7/12)*2.55^(2)*(r-r0)^2]
  with -2.55 in angstrom^(-1) and (7/12)*2.55^(2) in angstrom^(-2)
  These prefactors are converted here to the correct units
  */
  K3 = -2.55/force->angstrom;
  K4 = 7.0/12.0*2.55*2.55/force->angstrom/force->angstrom;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    dr2 = dr*dr;

    // force & energy

    de_bond = 2.0*k2[type]*dr*(1.0 + 1.5*K3*dr + 2.0*K4*dr2);
    if (r > 0.0) fbond = -de_bond/r;
    else fbond = 0.0;

    if (eflag) ebond = k2[type]*dr2*(1.0+K3*dr+K4*dr2);

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

void BondMM3::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(r0,n+1,"bond:r0");
  memory->create(k2,n+1,"bond:k2");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs from one line in input script or data file
------------------------------------------------------------------------- */

void BondMM3::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double k2_one = utils::numeric(FLERR,arg[1],false,lmp);
  double r0_one = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k2[i] = k2_one;
    r0[i] = r0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondMM3::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondMM3::write_restart(FILE *fp)
{
  fwrite(&k2[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondMM3::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k2[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&r0[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
  }
  MPI_Bcast(&k2[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondMM3::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k2[i],r0[i]);
}

/* ---------------------------------------------------------------------- */

double BondMM3::single(int type, double rsq,
                       int /* i */, int /* j */, double &fforce)
{
  /*
    E = K(r-r0)^2 [1-2.55*(r-r0)+(7/12)*2.55^(2)*(r-r0)^2]
    with -2.55 in angstrom^(-1) and (7/12)*2.55^(2) in angstrom^(-2)
    These prefactors are converted here to the correct units
  */
  double K3 = -2.55/force->angstrom;
  double K4 = 7.0/12.0*2.55*2.55/force->angstrom/force->angstrom;
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double dr2 = dr*dr;
  double de_bond = 2.0*k2[type]*dr*(1.0 + 1.5*K3*dr + 2.0*K4*dr2);
  if (r > 0.0) fforce = -de_bond/r;
  else fforce = 0.0;
  return k2[type]*dr2*(1.0+K3*dr+K4*dr2);
}
