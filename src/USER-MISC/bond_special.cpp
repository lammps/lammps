/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include "bond_special.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondSpecial::BondSpecial(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondSpecial::~BondSpecial()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(factor_lj);
    memory->destroy(factor_coul);
  }
}

/* ---------------------------------------------------------------------- */

void BondSpecial::init_style()
{
  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support bond style special");

  if (force->special_lj[1] != 0.0 || force->special_coul[1] != 0)
    error->all(FLERR,"Invalid 1-2 setting for bond style special.");

  if (force->special_angle != 1 && (force->special_lj[2] != 0.0 ||
                                    force->special_coul[3] != 0.0))
    error->all(FLERR,"Invalid 1-3 setting for bond style special.");

  if (force->special_dihedral != 1 && (force->special_lj[3] != 0.0 ||
                                       force->special_coul[3] != 0.0))
    error->all(FLERR,"Invalid 1-4 setting for bond style special.");

}

/* ---------------------------------------------------------------------- */

void BondSpecial::compute(int eflag, int vflag)
{
  int i1,i2,i1type,i2type,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq;

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

    i1type = atom->type[i1];
    i2type = atom->type[i2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;

    ebond = force->pair->single(i1,i2,i1type,i2type,rsq,factor_lj[type],
        factor_coul[type],fbond);

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

void BondSpecial::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(factor_lj,n+1,"bond:factor_lj");
  memory->create(factor_coul,n+1,"bond:factor_coul");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondSpecial::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nbondtypes,ilo,ihi);

  double factor_lj_one = force->numeric(FLERR,arg[1]);
  double factor_coul_one = force->numeric(FLERR,arg[2]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    factor_lj[i] = factor_lj_one;
    factor_coul[i] = factor_coul_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondSpecial::equilibrium_distance(int i)
{
  return 0.;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondSpecial::write_restart(FILE *fp)
{
  fwrite(&factor_lj[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&factor_coul[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondSpecial::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&factor_lj[1],sizeof(double),atom->nbondtypes,fp);
    fread(&factor_coul[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&factor_lj[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&factor_coul[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondSpecial::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,factor_lj[i],factor_coul[i]);
}

/* ---------------------------------------------------------------------- */

double BondSpecial::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  int itype = atom->type[i];
  int jtype = atom->type[j];
  double ebond;
  ebond = force->pair->single(i,j,itype,jtype,rsq,factor_lj[type],
      factor_coul[type],fforce);
  return ebond;
}


