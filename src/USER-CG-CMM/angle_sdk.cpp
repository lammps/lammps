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
   Variant of the harmonic angle potential for use with the
   lj/sdk potential for coarse grained MD simulations.
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "angle_sdk.h"
#include "atom.h"
#include "neighbor.h"
#include "pair.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "lj_sdk_common.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace LJSDKParms;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleSDK::AngleSDK(LAMMPS *lmp) : Angle(lmp) { repflag = 0;}

/* ---------------------------------------------------------------------- */

AngleSDK::~AngleSDK()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(theta0);
    memory->destroy(repscale);

    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- */

void AngleSDK::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3;
  double eangle,f1[3],f3[3],e13,f13;
  double dtheta,tk;
  double rsq1,rsq2,rsq3,r1,r2,c,s,a,a11,a12,a22;

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

    // 1-3 LJ interaction.
    // we only want to use the repulsive part,
    // and it can be scaled (or off).
    // so this has to be done here and not in the
    // general non-bonded code.

    f13 = e13 = delx3 = dely3 = delz3 = 0.0;

    if (repflag) {

      delx3 = x[i1][0] - x[i3][0];
      dely3 = x[i1][1] - x[i3][1];
      delz3 = x[i1][2] - x[i3][2];
      rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;

      const int type1 = atom->type[i1];
      const int type3 = atom->type[i3];

      f13=0.0;
      e13=0.0;

      if (rsq3 < rminsq[type1][type3]) {
        const int ljt = lj_type[type1][type3];
        const double r2inv = 1.0/rsq3;

        if (ljt == LJ12_4) {
          const double r4inv=r2inv*r2inv;

          f13 = r4inv*(lj1[type1][type3]*r4inv*r4inv - lj2[type1][type3]);
          if (eflag) e13 = r4inv*(lj3[type1][type3]*r4inv*r4inv - lj4[type1][type3]);

        } else if (ljt == LJ9_6) {
          const double r3inv = r2inv*sqrt(r2inv);
          const double r6inv = r3inv*r3inv;

          f13 = r6inv*(lj1[type1][type3]*r3inv - lj2[type1][type3]);
          if (eflag) e13 = r6inv*(lj3[type1][type3]*r3inv - lj4[type1][type3]);

        } else if (ljt == LJ12_6) {
          const double r6inv = r2inv*r2inv*r2inv;

          f13 = r6inv*(lj1[type1][type3]*r6inv - lj2[type1][type3]);
          if (eflag) e13 = r6inv*(lj3[type1][type3]*r6inv - lj4[type1][type3]);
        }

        // make sure energy is 0.0 at the cutoff.
        if (eflag) e13 -= emin[type1][type3];

        f13 *= r2inv;
      }
    }

    // force & energy

    dtheta = acos(c) - theta0[type];
    tk = k[type] * dtheta;

    if (eflag) eangle = tk*dtheta;

    a = -2.0 * tk * s;
    a11 = a*c / rsq1;
    a12 = -a / (r1*r2);
    a22 = a*c / rsq2;

    f1[0] = a11*delx1 + a12*delx2;
    f1[1] = a11*dely1 + a12*dely2;
    f1[2] = a11*delz1 + a12*delz2;
    f3[0] = a22*delx2 + a12*delx1;
    f3[1] = a22*dely2 + a12*dely1;
    f3[2] = a22*delz2 + a12*delz1;

    // apply force to each of the 3 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0] + f13*delx3;
      f[i1][1] += f1[1] + f13*dely3;
      f[i1][2] += f1[2] + f13*delz3;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0] - f13*delx3;
      f[i3][1] += f3[1] - f13*dely3;
      f[i3][2] += f3[2] - f13*delz3;
    }

    if (evflag) {
      ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
      if (repflag)
        ev_tally13(i1,i3,nlocal,newton_bond,e13,f13,delx3,dely3,delz3);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AngleSDK::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(theta0,n+1,"angle:theta0");
  memory->create(repscale,n+1,"angle:repscale");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleSDK::coeff(int narg, char **arg)
{
  if ((narg < 3) || (narg > 6))
    error->all(FLERR,"Incorrect args for angle coefficients");

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double theta0_one = force->numeric(FLERR,arg[2]);
  double repscale_one;

  // backward compatibility with old cg/cmm style input:
  // this had <lj_type> <epsilon> <sigma>
  // if epsilon is set to 0.0 we accept it as repscale 0.0
  // otherwise assume repscale 1.0, since we were using
  // epsilon to turn repulsion on or off.
  if (narg == 6) {
    repscale_one = force->numeric(FLERR,arg[4]);
    if (repscale_one > 0.0) repscale_one = 1.0;
  } else if (narg == 4) repscale_one = force->numeric(FLERR,arg[3]);
  else if (narg == 3) repscale_one = 1.0;
  else error->all(FLERR,"Incorrect args for angle coefficients");

  // convert theta0 from degrees to radians and store coefficients

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    theta0[i] = theta0_one/180.0 * MY_PI;
    repscale[i] = repscale_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ----------------------------------------------------------------------
   error check and initialize all values needed for force computation
------------------------------------------------------------------------- */

void AngleSDK::init_style()
{

  // make sure we use an SDK pair_style and that we need the 1-3 repulsion

  repflag = 0;
  for (int i = 1; i <= atom->nangletypes; i++)
    if (repscale[i] > 0.0) repflag = 1;

  // set up pointers to access SDK LJ parameters for 1-3 interactions

  if (repflag) {
    int itmp;
    if (force->pair == NULL)
      error->all(FLERR,"Angle style SDK requires use of a compatible with Pair style");

    lj1 = (double **) force->pair->extract("lj1",itmp);
    lj2 = (double **) force->pair->extract("lj2",itmp);
    lj3 = (double **) force->pair->extract("lj3",itmp);
    lj4 = (double **) force->pair->extract("lj4",itmp);
    lj_type = (int **) force->pair->extract("lj_type",itmp);
    rminsq = (double **) force->pair->extract("rminsq",itmp);
    emin = (double **) force->pair->extract("emin",itmp);

    if (!lj1 || !lj2 || !lj3 || !lj4 || !lj_type || !rminsq || !emin)
      error->all(FLERR,"Angle style SDK is incompatible with Pair style");
  }
}

/* ---------------------------------------------------------------------- */

double AngleSDK::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleSDK::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&repscale[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleSDK::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nangletypes,fp);
    fread(&theta0[1],sizeof(double),atom->nangletypes,fp);
    fread(&repscale[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&repscale[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleSDK::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],theta0[i]/MY_PI*180.0);
}

/* ---------------------------------------------------------------------- */

void AngleSDK::ev_tally13(int i, int j, int nlocal, int newton_bond,
                          double evdwl, double fpair,
                          double delx, double dely, double delz)
{
  double v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) {
        energy += evdwl;
      } else {
        if (i < nlocal)
          energy += 0.5*evdwl;
        if (j < nlocal)
          energy += 0.5*evdwl;
      }
    }
    if (eflag_atom) {
      if (newton_bond || i < nlocal) eatom[i] += 0.5*evdwl;
      if (newton_bond || j < nlocal) eatom[i] += 0.5*evdwl;
    }
  }

  if (vflag_either) {
    v[0] = delx*delx*fpair;
    v[1] = dely*dely*fpair;
    v[2] = delz*delz*fpair;
    v[3] = delx*dely*fpair;
    v[4] = delx*delz*fpair;
    v[5] = dely*delz*fpair;

    if (vflag_global) {
      if (newton_bond) {
        virial[0] += v[0];
        virial[1] += v[1];
        virial[2] += v[2];
        virial[3] += v[3];
        virial[4] += v[4];
        virial[5] += v[5];
      } else {
        if (i < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          virial[0] += 0.5*v[0];
          virial[1] += 0.5*v[1];
          virial[2] += 0.5*v[2];
          virial[3] += 0.5*v[3];
          virial[4] += 0.5*v[4];
          virial[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        vatom[i][0] += 0.5*v[0];
        vatom[i][1] += 0.5*v[1];
        vatom[i][2] += 0.5*v[2];
        vatom[i][3] += 0.5*v[3];
        vatom[i][4] += 0.5*v[4];
        vatom[i][5] += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        vatom[j][0] += 0.5*v[0];
        vatom[j][1] += 0.5*v[1];
        vatom[j][2] += 0.5*v[2];
        vatom[j][3] += 0.5*v[3];
        vatom[j][4] += 0.5*v[4];
        vatom[j][5] += 0.5*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double AngleSDK::single(int type, int i1, int i2, int i3)
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

  double e13=0.0;
  if (repflag) {

    // 1-3 LJ interaction.
    double delx3 = x[i1][0] - x[i3][0];
    double dely3 = x[i1][1] - x[i3][1];
    double delz3 = x[i1][2] - x[i3][2];
    domain->minimum_image(delx3,dely3,delz3);

    const int type1 = atom->type[i1];
    const int type3 = atom->type[i3];

    const double rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;

    if (rsq3 < rminsq[type1][type3]) {
      const int ljt = lj_type[type1][type3];
      const double r2inv = 1.0/rsq3;

      if (ljt == LJ12_4) {
        const double r4inv=r2inv*r2inv;

        e13 = r4inv*(lj3[type1][type3]*r4inv*r4inv - lj4[type1][type3]);

      } else if (ljt == LJ9_6) {
        const double r3inv = r2inv*sqrt(r2inv);
        const double r6inv = r3inv*r3inv;

        e13 = r6inv*(lj3[type1][type3]*r3inv - lj4[type1][type3]);

      } else if (ljt == LJ12_6) {
        const double r6inv = r2inv*r2inv*r2inv;

        e13 = r6inv*(lj3[type1][type3]*r6inv - lj4[type1][type3]);
      }

      // make sure energy is 0.0 at the cutoff.
      e13 -= emin[type1][type3];
    }
  }

  double dtheta = acos(c) - theta0[type];
  double tk = k[type] * dtheta;
  return tk*dtheta + e13;
}
