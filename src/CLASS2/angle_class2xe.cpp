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
   Contributing author: Josh Kemppainen (Michigan Technological University)
                        jdkemppa@mtu.edu
------------------------------------------------------------------------- */

#include "angle_class2xe.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "citeme.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 0.001;

static const char cite_style_classII_xe[] =
  "ClassII-xe Style: https://doi.org/10.1038/s41524-025-01838-5\n\n"
  "@Article{Kemppainen 2025,\n"
  " author = {J. Kemppainen and H. Heinz and G. M. Odegard},\n"
  " title = {Integrating complete bond dissociation in Class II force fields},\n"
  " journal = {npj Computational Materials},\n"
  " year =    2025,\n"
  " volume =  11,\n"
  " pages =   {341}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

AngleClass2xe::AngleClass2xe(LAMMPS *lmp) : Angle(lmp)
{
   if (lmp->citeme) lmp->citeme->add(cite_style_classII_xe);
}

/* ---------------------------------------------------------------------- */

AngleClass2xe::~AngleClass2xe()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(setflag_a);
    memory->destroy(setflag_bb);
    memory->destroy(setflag_ba);

    memory->destroy(theta0);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);

    memory->destroy(bb_d0);
    memory->destroy(bb_alpha);
    memory->destroy(bb_r1);
    memory->destroy(bb_r2);

    memory->destroy(ba_d1);
    memory->destroy(ba_d2);
    memory->destroy(ba_alpha1);
    memory->destroy(ba_alpha2);
    memory->destroy(ba_r1);
    memory->destroy(ba_r2);
  }
}

/* ---------------------------------------------------------------------- */

void AngleClass2xe::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double dtheta,dtheta2,dtheta3,dtheta4,de_angle;
  double dr1,dr2,tk1,tk2,aa1,aa2,aa11,aa12,aa21,aa22;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22,b1,b2;
  double vx11,vx12,vy11,vy12,vz11,vz12,vx21,vx22,vy21,vy22,vz21,vz22;
  double bb_ralpha1,bb_ralpha2;
  double ba_ralpha1,ba_ralpha2;

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

    // force & energy for angle term

    dtheta = acos(c) - theta0[type];
    dtheta2 = dtheta*dtheta;
    dtheta3 = dtheta2*dtheta;
    dtheta4 = dtheta3*dtheta;

    de_angle = 2.0*k2[type]*dtheta + 3.0*k3[type]*dtheta2 +
      4.0*k4[type]*dtheta3;

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

    if (eflag) eangle = k2[type]*dtheta2 + k3[type]*dtheta3 + k4[type]*dtheta4;

    // force & energy for bond-bond term

    dr1 = r1 - bb_r1[type];
    dr2 = r2 - bb_r2[type];
    bb_ralpha1 = exp(-bb_alpha[type]*dr1);
    bb_ralpha2 = exp(-bb_alpha[type]*dr2);
    tk1 = bb_d0[type]*bb_alpha[type]*bb_ralpha1*(1 - bb_ralpha2);
    tk2 = bb_d0[type]*bb_alpha[type]*bb_ralpha2*(1 - bb_ralpha1);

    f1[0] -= delx1*tk2/r1;
    f1[1] -= dely1*tk2/r1;
    f1[2] -= delz1*tk2/r1;

    f3[0] -= delx2*tk1/r2;
    f3[1] -= dely2*tk1/r2;
    f3[2] -= delz2*tk1/r2;

    if (eflag) eangle += bb_d0[type]*(1 - bb_ralpha1)*(1 - bb_ralpha2);

    // force & energy for bond-angle term

    dr1 = r1 - ba_r1[type];
    dr2 = r2 - ba_r2[type];
    ba_ralpha1 = exp(-ba_alpha1[type]*dr1);
    ba_ralpha2 = exp(-ba_alpha2[type]*dr2);
    aa1 = ba_d1[type]*(1 - ba_ralpha1)*s;
    aa2 = ba_d2[type]*(1 - ba_ralpha2)*s;

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

    b1 = ba_d1[type]*ba_alpha1[type]*ba_ralpha1*dtheta/r1;
    b2 = ba_d2[type]*ba_alpha2[type]*ba_ralpha2*dtheta/r2;

    f1[0] -= vx11 + b1*delx1 + vx12;
    f1[1] -= vy11 + b1*dely1 + vy12;
    f1[2] -= vz11 + b1*delz1 + vz12;

    f3[0] -= vx21 + b2*delx2 + vx22;
    f3[1] -= vy21 + b2*dely2 + vy22;
    f3[2] -= vz21 + b2*delz2 + vz22;

    if (eflag) eangle += ba_d1[type]*(1 - ba_ralpha1)*dtheta + ba_d2[type]*(1 - ba_ralpha2)*dtheta;

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

void AngleClass2xe::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(theta0,n+1,"angle:theta0");
  memory->create(k2,n+1,"angle:k2");
  memory->create(k3,n+1,"angle:k3");
  memory->create(k4,n+1,"angle:k4");

  memory->create(bb_d0,n+1,"angle:bb_d0");
  memory->create(bb_alpha,n+1,"angle:bb_alpha");
  memory->create(bb_r1,n+1,"angle:bb_r1");
  memory->create(bb_r2,n+1,"angle:bb_r2");

  memory->create(ba_d1,n+1,"angle:ba_d1");
  memory->create(ba_d2,n+1,"angle:ba_d2");
  memory->create(ba_alpha1,n+1,"angle:ba_alpha1");
  memory->create(ba_alpha2,n+1,"angle:ba_alpha2");
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

void AngleClass2xe::coeff(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  int count = 0;

  if (strcmp(arg[1],"bb") == 0) {
    if (narg != 6) error->all(FLERR,"Incorrect args for angle coefficients");

    double bb_d0_one = utils::numeric(FLERR,arg[2],false,lmp);
    double bb_alpha_one = utils::numeric(FLERR,arg[3],false,lmp);
    double bb_r1_one = utils::numeric(FLERR,arg[4],false,lmp);
    double bb_r2_one = utils::numeric(FLERR,arg[5],false,lmp);

    for (int i = ilo; i <= ihi; i++) {
      bb_d0[i] = bb_d0_one;
      bb_alpha[i] = bb_alpha_one;
      bb_r1[i] = bb_r1_one;
      bb_r2[i] = bb_r2_one;
      setflag_bb[i] = 1;
      count++;
    }

  } else if (strcmp(arg[1],"ba") == 0) {
    if (narg != 8) error->all(FLERR,"Incorrect args for angle coefficients");

    double ba_d1_one = utils::numeric(FLERR,arg[2],false,lmp);
    double ba_d2_one = utils::numeric(FLERR,arg[3],false,lmp);
    double ba_alpha1_one = utils::numeric(FLERR,arg[4],false,lmp);
    double ba_alpha2_one = utils::numeric(FLERR,arg[5],false,lmp);
    double ba_r1_one = utils::numeric(FLERR,arg[6],false,lmp);
    double ba_r2_one = utils::numeric(FLERR,arg[7],false,lmp);

    for (int i = ilo; i <= ihi; i++) {
      ba_d1[i] = ba_d1_one;
      ba_d2[i] = ba_d2_one;
      ba_alpha1[i] = ba_alpha1_one;
      ba_alpha2[i] = ba_alpha2_one;
      ba_r1[i] = ba_r1_one;
      ba_r2[i] = ba_r2_one;
      setflag_ba[i] = 1;
      count++;
    }

  } else {
    if (narg != 5) error->all(FLERR,"Incorrect args for angle coefficients");

    double theta0_one = utils::numeric(FLERR,arg[1],false,lmp);
    double k2_one = utils::numeric(FLERR,arg[2],false,lmp);
    double k3_one = utils::numeric(FLERR,arg[3],false,lmp);
    double k4_one = utils::numeric(FLERR,arg[4],false,lmp);

    // convert theta0 from degrees to radians

    for (int i = ilo; i <= ihi; i++) {
      theta0[i] = theta0_one/180.0 * MY_PI;
      k2[i] = k2_one;
      k3[i] = k3_one;
      k4[i] = k4_one;
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

double AngleClass2xe::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleClass2xe::write_restart(FILE *fp)
{
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k3[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&k4[1],sizeof(double),atom->nangletypes,fp);

  fwrite(&bb_d0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&bb_alpha[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&bb_r1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&bb_r2[1],sizeof(double),atom->nangletypes,fp);

  fwrite(&ba_d1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_d2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_alpha1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_alpha2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_r1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ba_r2[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleClass2xe::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&theta0[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k2[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k3[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&k4[1],sizeof(double),atom->nangletypes,fp,nullptr,error);

    utils::sfread(FLERR,&bb_d0[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&bb_alpha[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&bb_r1[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&bb_r2[1],sizeof(double),atom->nangletypes,fp,nullptr,error);

    utils::sfread(FLERR,&ba_d1[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&ba_d2[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&ba_alpha1[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&ba_alpha2[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&ba_r1[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
    utils::sfread(FLERR,&ba_r2[1],sizeof(double),atom->nangletypes,fp,nullptr,error);
  }

  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k3[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&k4[1],atom->nangletypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&bb_d0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb_alpha[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb_r1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&bb_r2[1],atom->nangletypes,MPI_DOUBLE,0,world);

  MPI_Bcast(&ba_d1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_d2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_alpha1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_alpha2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_r1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ba_r2[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleClass2xe::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",
            i,theta0[i]/MY_PI*180.0,k2[i],k3[i],k4[i]);

  fprintf(fp,"\nBondBond Coeffs\n\n");
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",i,bb_d0[i],bb_alpha[i],bb_r1[i],bb_r2[i]);

  fprintf(fp,"\nBondAngle Coeffs\n\n");
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g\n",i,ba_d1[i],ba_d2[i],ba_alpha1[i],ba_alpha2[i],ba_r1[i],ba_r2[i]);
}

/* ---------------------------------------------------------------------- */

double AngleClass2xe::single(int type, int i1, int i2, int i3)
{
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(FLERR, delx1,dely1,delz1);
  double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(FLERR, delx2,dely2,delz2);
  double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

  double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double dtheta = acos(c) - theta0[type];
  double dtheta2 = dtheta*dtheta;
  double dtheta3 = dtheta2*dtheta;
  double dtheta4 = dtheta3*dtheta;

  double energy = k2[type]*dtheta2 + k3[type]*dtheta3 + k4[type]*dtheta4;

  double dr1 = r1 - bb_r1[type];
  double dr2 = r2 - bb_r2[type];
  double bb_ralpha1 = exp(-bb_alpha[type] * dr1);
  double bb_ralpha2 = exp(-bb_alpha[type] * dr2);
  energy += bb_d0[type]*(1 - bb_ralpha1)*(1 - bb_ralpha2);

  dr1 = r1 - ba_r1[type];
  dr2 = r2 - ba_r2[type];
  double ba_ralpha1 = exp(-ba_alpha1[type]*dr1);
  double ba_ralpha2 = exp(-ba_alpha2[type]*dr2);
  energy += ba_d1[type]*(1 - ba_ralpha1)*dtheta + ba_d2[type]*(1 - ba_ralpha2)*dtheta;

  return energy;
}
