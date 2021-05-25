// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Steven Vandenbrande
------------------------------------------------------------------------- */

#include "angle_cross.h"

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

AngleCross::AngleCross(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCross::~AngleCross()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(kss);
    memory->destroy(kbs0);
    memory->destroy(kbs1);
    memory->destroy(r00);
    memory->destroy(r01);
    memory->destroy(theta0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCross::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta;
  double dr1,dr2,tk1,tk2,aa1,aa2,aa11,aa12,aa21,aa22;
  double rsq1,rsq2,r1,r2,c,s,b1,b2;
  double vx11,vx12,vy11,vy12,vz11,vz12,vx21,vx22,vy21,vy22,vz21,vz22;

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

    // angle (cos and sin)

    c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // force & energy for bond-bond term
    dr1 = r1 - r00[type];
    dr2 = r2 - r01[type];
    tk1 = kss[type] * dr1;
    tk2 = kss[type] * dr2;

    f1[0] = -delx1*tk2/r1;
    f1[1] = -dely1*tk2/r1;
    f1[2] = -delz1*tk2/r1;

    f3[0] = -delx2*tk1/r2;
    f3[1] = -dely2*tk1/r2;
    f3[2] = -delz2*tk1/r2;

    if (eflag) eangle = kss[type]*dr1*dr2;

    // force & energy for bond-angle term
    dtheta = acos(c) - theta0[type];

    aa1 = s * dr1 * kbs0[type];
    aa2 = s * dr2 * kbs1[type];

    aa11 = aa1 * c / rsq1;
    aa12 = -aa1 / (r1 * r2);
    aa21 = aa2 * c / rsq1;
    aa22 = -aa2 / (r1 * r2);

    vx11 = (aa11 * delx1) + (aa12 * delx2);
    vx12 = (aa21 * delx1) + (aa22 * delx2);
    vy11 = (aa11 * dely1) + (aa12 * dely2);
    vy12 = (aa21 * dely1) + (aa22 * dely2);
    vz11 = (aa11 * delz1) + (aa12 * delz2);
    vz12 = (aa21 * delz1) + (aa22 * delz2);

    aa11 = aa1 * c / rsq2;
    aa21 = aa2 * c / rsq2;

    vx21 = (aa11 * delx2) + (aa12 * delx1);
    vx22 = (aa21 * delx2) + (aa22 * delx1);
    vy21 = (aa11 * dely2) + (aa12 * dely1);
    vy22 = (aa21 * dely2) + (aa22 * dely1);
    vz21 = (aa11 * delz2) + (aa12 * delz1);
    vz22 = (aa21 * delz2) + (aa22 * delz1);

    b1 = kbs0[type] * dtheta / r1;
    b2 = kbs1[type] * dtheta / r2;

    f1[0] -= vx11 + b1*delx1 + vx12;
    f1[1] -= vy11 + b1*dely1 + vy12;
    f1[2] -= vz11 + b1*delz1 + vz12;

    f3[0] -= vx21 + b2*delx2 + vx22;
    f3[1] -= vy21 + b2*dely2 + vy22;
    f3[2] -= vz21 + b2*delz2 + vz22;

    if (eflag) eangle += kbs0[type]*dr1*dtheta + kbs1[type]*dr2*dtheta;

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

void AngleCross::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(kss,n+1,"angle:kss");
  memory->create(kbs0,n+1,"angle:kbs0");
  memory->create(kbs1,n+1,"angle:kbs1");
  memory->create(r00,n+1,"angle:r00");
  memory->create(r01,n+1,"angle:r01");
  memory->create(theta0,n+1,"angle:theta0");
  memory->create(setflag,n+1,"angle:setflag");

  for (int i = 1; i <= n; i++)
    setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs
------------------------------------------------------------------------- */

void AngleCross::coeff(int narg, char **arg)
{
  if (narg != 7) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  int count = 0;

    double kss_one = utils::numeric(FLERR,arg[1],false,lmp);
    double kbs0_one = utils::numeric(FLERR,arg[2],false,lmp);
    double kbs1_one = utils::numeric(FLERR,arg[3],false,lmp);
    double r0_one = utils::numeric(FLERR,arg[4],false,lmp);
    double r1_one = utils::numeric(FLERR,arg[5],false,lmp);
    double theta0_one = utils::numeric(FLERR,arg[6],false,lmp);

    for (int i = ilo; i <= ihi; i++) {
      kss[i] = kss_one;
      kbs0[i] = kbs0_one;
      kbs1[i] = kbs1_one;
      r00[i] = r0_one;
      r01[i] = r1_one;
      // Convert theta0 from degrees to radians
      theta0[i] = theta0_one*MY_PI/180.0;
      setflag[i] = 1;
      count++;
    }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleCross::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCross::write_restart(FILE *fp)
{
  fwrite(&kss[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kbs0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kbs1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&r00[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&r01[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleCross::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&kss[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&kbs0[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&kbs1[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&r00[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&r01[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&theta0[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
  }

  MPI_Bcast(&kss[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kbs0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kbs1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r00[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r01[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleCross::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g\n",
            i,kss[i],kbs0[i],kbs1[i],r00[i],r01[i],theta0[i]/MY_PI*180.0);
}

/* ---------------------------------------------------------------------- */

double AngleCross::single(int type, int i1, int i2, int i3)
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

  double s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  double dtheta = acos(c) - theta0[type];
  double dr1 = r1 - r00[type];
  double dr2 = r2 - r01[type];
  double energy = kss[type]*dr1*dr2+kbs0[type]*dr1*dtheta + kbs1[type]*dr2*dtheta;
  return energy;
}
