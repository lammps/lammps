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

/* ----------------------------------------------------------------------
   Contributing author: Hendrik Heenen (Technical University of Munich)
                        and Rochus Schmid (Ruhr-Universitaet Bochum)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include "angle_class2_p6.h"
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

AngleClass2P6::AngleClass2P6(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleClass2P6::~AngleClass2P6()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(setflag_a);
    memory->destroy(setflag_bb);
    memory->destroy(setflag_ba);

    memory->destroy(theta0);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);
    memory->destroy(k5);
    memory->destroy(k6);

    memory->destroy(bb_k);
    memory->destroy(bb_r1);
    memory->destroy(bb_r2);

    memory->destroy(ba_k1);
    memory->destroy(ba_k2);
    memory->destroy(ba_r1);
    memory->destroy(ba_r2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleClass2P6::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta,dtheta2,dtheta3,dtheta4,dtheta5,dtheta6,de_angle;
  double dr1,dr2,tk1,tk2,aa1,aa2,aa11,aa12,aa21,aa22;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22,b1,b2;
  double vx11,vx12,vy11,vy12,vz11,vz12,vx21,vx22,vy21,vy22,vz21,vz22;

  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

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

    // force & energy for angle term

    dtheta = acos(c) - theta0[type];
    dtheta2 = dtheta*dtheta;
    dtheta3 = dtheta2*dtheta;
    dtheta4 = dtheta3*dtheta;
    dtheta5 = dtheta4*dtheta;
    dtheta6 = dtheta5*dtheta;

    de_angle = 2.0*k2[type]*dtheta + 3.0*k3[type]*dtheta2 +
      4.0*k4[type]*dtheta3 + 5.0*k5[type]*dtheta4 +
      6.0*k6[type]*dtheta5;

    a = -de_angle*s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;

    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    if (eflag) eangle = k2[type]*dtheta2 + k3[type]*dtheta3 + k4[type]*dtheta4
                 + k5[type]*dtheta5 + k6[type]*dtheta6;

    // force & energy for bond-bond term

    dr1 = r1 - bb_r1[type];
    dr2 = r2 - bb_r2[type];
    tk1 = bb_k[type] * dr1;
    tk2 = bb_k[type] * dr2;

    f1[0] -= delx1*tk2/r1;
    f1[1] -= dely1*tk2/r1;
    f1[2] -= delz1*tk2/r1;

    f3[0] -= delx2*tk1/r2;
    f3[1] -= dely2*tk1/r2;
    f3[2] -= delz2*tk1/r2;

    if (eflag) eangle += bb_k[type]*dr1*dr2;

    // force & energy for bond-angle term

    aa1 = s * dr1 * ba_k1[type];
    aa2 = s * dr2 * ba_k2[type];

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

    b1 = ba_k1[type] * dtheta / r1;
    b2 = ba_k2[type] * dtheta / r2;

    f1[0] -= vx11 + b1*delx1 + vx12;
    f1[1] -= vy11 + b1*dely1 + vy12;
    f1[2] -= vz11 + b1*delz1 + vz12;

    f3[0] -= vx21 + b2*delx2 + vx22;
    f3[1] -= vy21 + b2*dely2 + vy22;
    f3[2] -= vz21 + b2*delz2 + vz22;

    if (eflag) eangle += ba_k1[type]*dr1*dtheta + ba_k2[type]*dr2*dtheta;

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

void AngleClass2P6::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(theta0,n+1,"angle:theta0");
  memory->create(k2,n+1,"angle:k2");
  memory->create(k3,n+1,"angle:k3");
  memory->create(k4,n+1,"angle:k4");
  memory->create(k5,n+1,"angle:k5");
  memory->create(k6,n+1,"angle:k6");

  memory->create(bb_k,n+1,"angle:bb_k");
  memory->create(bb_r1,n+1,"angle:bb_r1");
  memory->create(bb_r2,n+1,"angle:bb_r2");

  memory->create(ba_k1,n+1,"angle:ba_k1");
  memory->create(ba_k2,n+1,"angle:ba_k2");
  memory->create(ba_r1,n+1,"angle:ba_r1");
  memory->create(ba_r2,n+1,"angle:ba_r2");

  memory->create(setflag,n+1,"angle:setflag");
  memory->create(setflag_a,n+1,"angle:setflag_a");
  memory->create(setflag_bb,n+1,"angle:setflag_bb");
  memory->create(setflag_ba,n+1,"angle:setflag_ba");
  for (int i = 1; i <= n; i++)
    setflag[i] = setflag_a[i] = setflag_bb[i] = setflag_ba[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
   arg1 = "bb" -> BondBond coeffs
   arg1 = "ba" -> BondAngle coeffs
   else -> Angle coeffs
------------------------------------------------------------------------- */

void AngleClass2P6::coeff(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  int count = 0;

  if (strcmp(arg[1],"bb") == 0) {
    if (narg != 5) error->all(FLERR,"Incorrect args for angle coefficients");

    double bb_k_one = force->numeric(FLERR,arg[2]);
    double bb_r1_one = force->numeric(FLERR,arg[3]);
    double bb_r2_one = force->numeric(FLERR,arg[4]);

    for (int i = ilo; i <= ihi; i++) {
      bb_k[i] = bb_k_one;
      bb_r1[i] = bb_r1_one;
      bb_r2[i] = bb_r2_one;
      setflag_bb[i] = 1;
      count++;
    }

  } else if (strcmp(arg[1],"ba") == 0) {
    if (narg != 6) error->all(FLERR,"Incorrect args for angle coefficients");

    double ba_k1_one = force->numeric(FLERR,arg[2]);
    double ba_k2_one = force->numeric(FLERR,arg[3]);
    double ba_r1_one = force->numeric(FLERR,arg[4]);
    double ba_r2_one = force->numeric(FLERR,arg[5]);

    for (int i = ilo; i <= ihi; i++) {
      ba_k1[i] = ba_k1_one;
      ba_k2[i] = ba_k2_one;
      ba_r1[i] = ba_r1_one;
      ba_r2[i] = ba_r2_one;
      setflag_ba[i] = 1;
      count++;
    }

  } else {
    if (narg != 7) error->all(FLERR,"Incorrect args for angle coefficients");

    double theta0_one = force->numeric(FLERR,arg[1]);
    double k2_one = force->numeric(FLERR,arg[2]);
    double k3_one = force->numeric(FLERR,arg[3]);
    double k4_one = force->numeric(FLERR,arg[4]);
    double k5_one = force->numeric(FLERR,arg[5]);
    double k6_one = force->numeric(FLERR,arg[6]);

    // convert theta0 from degrees to radians

    for (int i = ilo; i <= ihi; i++) {
      theta0[i] = theta0_one/180.0 * MY_PI;
      k2[i] = k2_one;
      k3[i] = k3_one;
      k4[i] = k4_one;
      k5[i] = k5_one;
      k6[i] = k6_one;
      setflag_a[i] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");

  for (int i = ilo; i <= ihi; i++)
    if (setflag_a[i] == 1 && setflag_bb[i] == 1 && setflag_ba[i] == 1)
      setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double AngleClass2P6::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleClass2P6::write_restart(FILE *fp)
{
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k3[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k4[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k5[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k6[1],sizeof(double),atom->nangletypes,fp);

  fwrite(&bb_k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&bb_r1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&bb_r2[1],sizeof(double),atom->nangletypes,fp);

  fwrite(&ba_k1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_k2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_r1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_r2[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleClass2P6::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&theta0[1],sizeof(double),atom->nangletypes,fp);
    fread(&k2[1],sizeof(double),atom->nangletypes,fp);
    fread(&k3[1],sizeof(double),atom->nangletypes,fp);
    fread(&k4[1],sizeof(double),atom->nangletypes,fp);
    fread(&k5[1],sizeof(double),atom->nangletypes,fp);
    fread(&k6[1],sizeof(double),atom->nangletypes,fp);

    fread(&bb_k[1],sizeof(double),atom->nangletypes,fp);
    fread(&bb_r1[1],sizeof(double),atom->nangletypes,fp);
    fread(&bb_r2[1],sizeof(double),atom->nangletypes,fp);

    fread(&ba_k1[1],sizeof(double),atom->nangletypes,fp);
    fread(&ba_k2[1],sizeof(double),atom->nangletypes,fp);
    fread(&ba_r1[1],sizeof(double),atom->nangletypes,fp);
    fread(&ba_r2[1],sizeof(double),atom->nangletypes,fp);
  }

  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k4[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k5[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k6[1],atom->nangletypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&bb_k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb_r1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb_r2[1],atom->nangletypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&ba_k1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_k2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_r1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_r2[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleClass2P6::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g\n",
            i,theta0[i]/MY_PI*180.0,k2[i],k3[i],k4[i],k5[i],k6[i]);

  fprintf(fp,"\nBondBond Coeffs\n\n");
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,bb_k[i],bb_r1[i],bb_r2[i]);

  fprintf(fp,"\nBondAngle Coeffs\n\n");
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,ba_k1[i],ba_k2[i],ba_r1[i],ba_r2[i]);
}

/* ---------------------------------------------------------------------- */

double AngleClass2P6::single(int type, int i1, int i2, int i3)
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
  double dtheta2 = dtheta*dtheta;
  double dtheta3 = dtheta2*dtheta;
  double dtheta4 = dtheta3*dtheta;
  double dtheta5 = dtheta4*dtheta;
  double dtheta6 = dtheta5*dtheta;

  double energy = k2[type]*dtheta2 + k3[type]*dtheta3 + k4[type]*dtheta4
           + k5[type]*dtheta5 + k6[type]*dtheta6;

  double dr1 = r1 - bb_r1[type];
  double dr2 = r2 - bb_r2[type];
  energy += bb_k[type]*dr1*dr2;

  energy += ba_k1[type]*dr1*dtheta + ba_k2[type]*dr2*dtheta;
  return energy;
}
