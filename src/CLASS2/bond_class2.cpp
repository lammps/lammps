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
   Contributing author: Eric Simon (Cray)
------------------------------------------------------------------------- */

#include "bond_class2.h"

#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondClass2::BondClass2(LAMMPS *lmp) : Bond(lmp)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

BondClass2::~BondClass2()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(r0);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);
  }
}

/* ---------------------------------------------------------------------- */

void BondClass2::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,dr2,dr3,dr4,de_bond;

  ebond = 0.0;
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

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    dr2 = dr*dr;
    dr3 = dr2*dr;
    dr4 = dr3*dr;

    // force & energy

    de_bond = 2.0*k2[type]*dr + 3.0*k3[type]*dr2 + 4.0*k4[type]*dr3;
    if (r > 0.0) fbond = -de_bond/r;
    else fbond = 0.0;

    if (eflag) ebond = k2[type]*dr2 + k3[type]*dr3 + k4[type]*dr4;

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

void BondClass2::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(r0,n+1,"bond:r0");
  memory->create(k2,n+1,"bond:k2");
  memory->create(k3,n+1,"bond:k3");
  memory->create(k4,n+1,"bond:k4");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs from one line in input script or data file
------------------------------------------------------------------------- */

void BondClass2::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nbondtypes,ilo,ihi,error);

  double r0_one = utils::numeric(FLERR,arg[1],false,lmp);
  double k2_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k3_one = utils::numeric(FLERR,arg[3],false,lmp);
  double k4_one = utils::numeric(FLERR,arg[4],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    r0[i] = r0_one;
    k2[i] = k2_one;
    k3[i] = k3_one;
    k4[i] = k4_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondClass2::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondClass2::write_restart(FILE *fp)
{
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k3[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&k4[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondClass2::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&r0[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&k2[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&k3[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
    utils::sfread(FLERR,&k4[1],sizeof(double),atom->nbondtypes,fp,nullptr,error);
  }
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k4[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondClass2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,r0[i],k2[i],k3[i],k4[i]);
}

/* ---------------------------------------------------------------------- */

double BondClass2::single(int type, double rsq, int /*i*/, int /*j*/, double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double dr2 = dr*dr;
  double dr3 = dr2*dr;
  double dr4 = dr3*dr;
  double de_bond = 2.0*k2[type]*dr + 3.0*k3[type]*dr2 + 4.0*k4[type]*dr3;
  if (r > 0.0) fforce = -de_bond/r;
  else fforce = 0.0;
  return (k2[type]*dr2 + k3[type]*dr3 + k4[type]*dr4);
}

/* ---------------------------------------------------------------------- */

void BondClass2::born_matrix(int type, double rsq, int /*i*/, int /*j*/, double &du, double &du2)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  du = 0.0;
  du2 = 2*k2[type] + 6*k3[type]*dr + 12*k4[type]*dr*dr;
  if (r > 0.0) du = 2*k2[type]*dr + 3*k3[type]*dr*dr + 4*k4[type]*dr*dr*dr;
}

/* ---------------------------------------------------------------------- */

void *BondClass2::extract(const char *str, int &dim)
{
  dim = 1;
  if (strcmp(str,"r0")==0) return (void*) r0;
  return nullptr;
}
