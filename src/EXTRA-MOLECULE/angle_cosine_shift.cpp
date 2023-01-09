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
   Contributing author: Carsten Svaneborg, science@zqex.dk
------------------------------------------------------------------------- */

#include "angle_cosine_shift.h"

#include <cmath>
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCosineShift::AngleCosineShift(LAMMPS *lmp) : Angle(lmp)
{
  kcost = nullptr;
}

/* ---------------------------------------------------------------------- */

AngleCosineShift::~AngleCosineShift()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(kcost);
    memory->destroy(ksint);
    memory->destroy(theta);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineShift::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,c,s,cps,kcos,ksin,a11,a12,a22;

  eangle = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];

    // 1st bond

    delx1 = x[i1][0] - x[i2][0];
    dely1 = x[i1][1] - x[i2][1];
    delz1 = x[i1][2] - x[i2][2];

    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    r1 = sqrt(rsq1);

    // 2nd bond

    delx2 = x[i3][0] - x[i2][0];
    dely2 = x[i3][1] - x[i2][1];
    delz2 = x[i3][2] - x[i2][2];

    rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    r2 = sqrt(rsq2);

    // c = cosine of angle
    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // C= sine of angle
    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;

    // force & energy
    kcos=kcost[type];
    ksin=ksint[type];
    if (eflag) eangle = -k[type]-kcos*c-ksin*s;

    cps = c/s;          // NOTE absorbed one c

    a11 = (-kcos +ksin*cps )*c/ rsq1;
    a12 = ( kcos -ksin*cps )  / (r1*r2);
    a22 = (-kcos +ksin*cps )*c/ rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineShift::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(ksint,n+1,"angle:ksint");
  memory->create(kcost,n+1,"angle:kcost");
  memory->create(theta,n+1,"angle:theta");
  memory->create(setflag,n+1,"angle:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void AngleCosineShift::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  double umin   = utils::numeric(FLERR,arg[1],false,lmp);
  double theta0 = utils::numeric(FLERR,arg[2],false,lmp);

// k=Umin/2

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = umin/2;
    kcost[i] = umin/2*cos(theta0*MY_PI / 180.0);
    ksint[i] = umin/2*sin(theta0*MY_PI / 180.0);
    theta[i] = theta0*MY_PI / 180.0;

    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleCosineShift::equilibrium_angle(int i)
{
  return theta[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCosineShift::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kcost[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ksint[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleCosineShift::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0)
    {
      utils::sfread(FLERR,&k[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
      utils::sfread(FLERR,&kcost[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
      utils::sfread(FLERR,&ksint[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
      utils::sfread(FLERR,&theta[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    }
  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kcost[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ksint[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleCosineShift::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g\n",i,2.0*k[i],theta[i]/MY_PI*180.0);
}

/* ---------------------------------------------------------------------- */

double AngleCosineShift::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  double s=sqrt(1.0-c*c);

  return -k[type]-kcost[type]*c-ksint[type]*s;
}
