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

#include "angle_gaussian.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMAL 0.001
#define SMALL 1.0e-8

/* ---------------------------------------------------------------------- */

AngleGaussian::AngleGaussian(LAMMPS *lmp)
  : Angle(lmp), nterms(nullptr), angle_temperature(nullptr),
    alpha(nullptr), width(nullptr), theta0(nullptr)
{
}

/* ---------------------------------------------------------------------- */

AngleGaussian::~AngleGaussian()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(nterms);
    memory->destroy(angle_temperature);
    for (int i = 1; i <= atom->nangletypes; i++) {
      delete [] alpha[i];
      delete [] width[i];
      delete [] theta0[i];
    }
    delete [] alpha;
    delete [] width;
    delete [] theta0;
  }
}

/* ---------------------------------------------------------------------- */

void AngleGaussian::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
  double prefactor, exponent, g_i, sum_g_i, sum_numerator;

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
    if (s < SMAL) s = SMAL;
    s = 1.0/s;

    // force & energy
    double theta = acos(c);

    sum_g_i = 0.0;
    sum_numerator = 0.0;
    for (int i = 0; i < nterms[type]; i++) {
      dtheta = theta - theta0[type][i];
      prefactor = (alpha[type][i]/(width[type][i]*sqrt(MY_PI2)));
      exponent = -2*dtheta*dtheta/(width[type][i]*width[type][i]);
      g_i = prefactor*exp(exponent);
      sum_g_i += g_i;
      sum_numerator += g_i*dtheta/(width[type][i]*width[type][i]);
    }

    if (sum_g_i < SMALL) sum_g_i = SMALL;
    if (eflag) eangle = -(force->boltz*angle_temperature[type])*log(sum_g_i);

    // I should check about the sign of this expression
    a = -4.0*(force->boltz*angle_temperature[type])*(sum_numerator/sum_g_i)*s;
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

void AngleGaussian::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(nterms,n+1,"angle:nterms");
  memory->create(angle_temperature,n+1,"angle:angle_temperature");

  alpha = new double *[n+1];
  width = new double *[n+1];
  theta0 = new double *[n+1];
  memset(alpha,0,sizeof(double)*(n+1));
  memset(width,0,sizeof(double)*(n+1));
  memset(theta0,0,sizeof(double)*(n+1));

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleGaussian::coeff(int narg, char **arg)
{
  if (narg < 6) error->all(FLERR,"Incorrect args for angle coefficients");

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  double angle_temperature_one = utils::numeric(FLERR,arg[1],false,lmp);
  int n = utils::inumeric(FLERR,arg[2],false,lmp);
  if (narg != 3*n + 3)
    error->all(FLERR,"Incorrect args for angle coefficients");

  if (!allocated) allocate();

  // convert theta0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    angle_temperature[i] = angle_temperature_one;
    nterms[i] = n;
    delete[] alpha[i];
    alpha[i] = new double [n];
    delete[] width[i];
    width[i] = new double [n];
    delete[] theta0[i];
    theta0[i] = new double [n];
    for (int j = 0; j < n; j++) {
      alpha[i][j] = utils::numeric(FLERR,arg[3+3*j],false,lmp);
      width[i][j] = utils::numeric(FLERR,arg[4+3*j],false,lmp);
      theta0[i][j] = utils::numeric(FLERR,arg[5+3*j],false,lmp)* MY_PI / 180.0;
      setflag[i] = 1;
    }
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleGaussian::equilibrium_angle(int i)
{
  return theta0[i][0];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleGaussian::write_restart(FILE *fp)
{
  fwrite(&angle_temperature[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&nterms[1],sizeof(int),atom->nangletypes,fp);
  for (int i = 1; i <= atom->nangletypes; i++) {
    fwrite(alpha[i],sizeof(double),nterms[i],fp);
    fwrite(width[i],sizeof(double),nterms[i],fp);
    fwrite(theta0[i],sizeof(double),nterms[i],fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleGaussian::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&angle_temperature[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&nterms[1],sizeof(int),atom->nangletypes,fp,nullptr,error);
  }
  MPI_Bcast(&angle_temperature[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&nterms[1],atom->nangletypes,MPI_INT,0,world);

  // allocate
  for (int i = 1; i <= atom->nangletypes; i++) {
    alpha[i] = new double [nterms[i]];
    width[i] = new double [nterms[i]];
    theta0[i] = new double [nterms[i]];
  }

  if (comm->me == 0) {
    for (int i = 1; i <= atom->nangletypes; i++) {
      utils::sfread(FLERR,alpha[i],sizeof(double),nterms[i],fp,nullptr,error);
      utils::sfread(FLERR,width[i],sizeof(double),nterms[i],fp,nullptr,error);
      utils::sfread(FLERR,theta0[i],sizeof(double),nterms[i],fp,nullptr,error);
    }
  }

  for (int i = 1; i <= atom->nangletypes; i++) {
    MPI_Bcast(alpha[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(width[i],nterms[i],MPI_DOUBLE,0,world);
    MPI_Bcast(theta0[i],nterms[i],MPI_DOUBLE,0,world);
  }

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleGaussian::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++) {
    fprintf(fp,"%d %g %d",i,angle_temperature[i],nterms[i]);
    for (int j = 0; j < nterms[i]; j++) {
      fprintf(fp," %g %g %g",alpha[i][j],width[i][j],(theta0[i][j]/MY_PI)*180.0);
    }
    fprintf(fp, "\n");
  }

}

/* ---------------------------------------------------------------------- */

double AngleGaussian::single(int type, int i1, int i2, int i3)
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
  double theta = acos(c) ;

  double sum_g_i = 0.0;
  for (int i = 0; i < nterms[type]; i++) {
    double dtheta = theta - theta0[type][i];
    double prefactor = (alpha[type][i]/(width[type][i]*sqrt(MY_PI2)));
    double exponent = -2*dtheta*dtheta/(width[type][i]*width[type][i]);
    double g_i = prefactor*exp(exponent);
    sum_g_i += g_i;
  }

  if (sum_g_i < SMALL) sum_g_i = SMALL;
  return -(force->boltz*angle_temperature[type])*log(sum_g_i);
}
