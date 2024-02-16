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
   Contributing author: Loukas D. Peristeras (Scienomics SARL)
   [ based on angle_harmonic.cpp]
------------------------------------------------------------------------- */

#include "angle_quartic.h"

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

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

AngleQuartic::AngleQuartic(LAMMPS *lmp) : Angle(lmp)
{
  born_matrix_enable = 1;
}

/* ---------------------------------------------------------------------- */

AngleQuartic::~AngleQuartic()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);
    memory->destroy(theta0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleQuartic::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta,dtheta2,dtheta3,dtheta4,tk;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;

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

    // force & energy

    dtheta = acos(c) - theta0[type];
    dtheta2 = dtheta * dtheta;
    dtheta3 = dtheta2 * dtheta;
    tk =  2.0 * k2[type] * dtheta + 3.0 * k3[type] * dtheta2 +
      4.0 * k4[type] * dtheta3;

    if (eflag) {
      dtheta4 = dtheta3 * dtheta;
      eangle = k2[type] * dtheta2 + k3[type] * dtheta3 + k4[type] * dtheta4;
    }

    a = -tk * s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

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

void AngleQuartic::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k2,n+1,"angle:k2");
  memory->create(k3,n+1,"angle:k3");
  memory->create(k4,n+1,"angle:k4");
  memory->create(theta0,n+1,"angle:theta0");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleQuartic::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  double theta0_one = utils::numeric(FLERR,arg[1],false,lmp);
  double k2_one = utils::numeric(FLERR,arg[2],false,lmp);
  double k3_one = utils::numeric(FLERR,arg[3],false,lmp);
  double k4_one = utils::numeric(FLERR,arg[4],false,lmp);

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k2[i] = k2_one;
    k3[i] = k3_one;
    k4[i] = k4_one;
    theta0[i] = theta0_one/180.0 * MY_PI;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleQuartic::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleQuartic::write_restart(FILE *fp)
{
  fwrite(&k2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k3[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k4[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleQuartic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k2[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k3[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k4[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&theta0[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
  }
  MPI_Bcast(&k2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k4[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleQuartic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,theta0[i]/MY_PI*180.0,k2[i],k3[i],k4[i]);
}

/* ---------------------------------------------------------------------- */

double AngleQuartic::single(int type, int i1, int i2, int i3)
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

  double dtheta = acos(c) - theta0[type];
  double dtheta2 = dtheta * dtheta;
  double dtheta3 = dtheta2 * dtheta;
  double dtheta4 = dtheta3 * dtheta;
  return k2[type] * dtheta2 + k3[type] * dtheta3 + k4[type] * dtheta4;
}

/* ---------------------------------------------------------------------- */

void AngleQuartic::born_matrix(int type, int i1, int i2, int i3, double &du, double &du2)
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
  double theta = acos(c);

  double s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;

  double dtheta = theta - theta0[type];
  double dtheta2 = dtheta * dtheta;
  double dtheta3 = dtheta2 * dtheta;

  du = -(2.0 * k2[type] * dtheta + 3.0 * k3[type] * dtheta2 + 4.0 * k4[type] * dtheta3) / s;
  du2 = (2.0 * k2[type] + 6.0 * k3[type] * dtheta + 12.0 * k4[type] * dtheta2) / (s*s) -
          (2.0 * k2[type] * dtheta + 3.0 * k3[type] * dtheta2 + 4.0 * k4[type] * dtheta3) * c / (s*s*s);
}
