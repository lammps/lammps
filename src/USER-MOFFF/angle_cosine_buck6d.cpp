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
#include <cstdlib>
#include "angle_cosine_buck6d.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleCosineBuck6d::AngleCosineBuck6d(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCosineBuck6d::~AngleCosineBuck6d()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(multiplicity);
    memory->destroy(th0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCosineBuck6d::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type,itype,jtype;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double eangle,f1[3],f3[3];
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
  double tk;

  // extra lj variables
  double delx3,dely3,delz3,rsq3,r3;
  double rexp,r32inv,r36inv,r314inv,forcebuck6d,fpair;
  double term1,term2,term3,term4,term5,ebuck6d,evdwl;
  double rcu,rqu,sme,smf;

  eangle = evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  // insure pair->ev_tally() will use 1-3 virial contribution

  if (vflag_global == 2)
    force->pair->vflag_either = force->pair->vflag_global = 1;

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int newton_bond = force->newton_bond;
  int *atomtype = atom->type;

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

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;
    s = 1.0/s;

    // force & energy

    // explicit lj-contribution

    itype = atomtype[i1];
    jtype = atomtype[i3];

    delx3 = x[i1][0] - x[i3][0];
    dely3 = x[i1][1] - x[i3][1];
    delz3 = x[i1][2] - x[i3][2];
    rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;

    if (rsq3 < cut_ljsq[itype][jtype]) {
      r3 = sqrt(rsq3);
      r32inv = 1.0/rsq3;
      r36inv = r32inv*r32inv*r32inv;
      r314inv = r36inv*r36inv*r32inv;
      rexp = exp(-r3*buck6d2[itype][jtype]);
      term1 = buck6d3[itype][jtype]*r36inv;
      term2 = buck6d4[itype][jtype]*r314inv;
      term3 = term2*term2;
      term4 = 1.0/(1.0 + term2);
      term5 = 1.0/(1.0 + 2.0*term2 + term3);
      forcebuck6d = buck6d1[itype][jtype]*buck6d2[itype][jtype]*r3*rexp;
      forcebuck6d -= term1*(6.0*term4 - term5*14.0*term2);
      ebuck6d = buck6d1[itype][jtype]*rexp - term1*term4;

      // smoothing term
      if (rsq3 > rsmooth_sq[itype][jtype]) {
        rcu = r3*rsq3;
        rqu = rsq3*rsq3;
        sme = c5[itype][jtype]*rqu*r3 + c4[itype][jtype]*rqu + c3[itype][jtype]*rcu +
              c2[itype][jtype]*rsq3 + c1[itype][jtype]*r3 + c0[itype][jtype];
        smf = 5.0*c5[itype][jtype]*rqu + 4.0*c4[itype][jtype]*rcu +
              3.0*c3[itype][jtype]*rsq3 + 2.0*c2[itype][jtype]*r3 + c1[itype][jtype];
        forcebuck6d = forcebuck6d*sme + ebuck6d*smf;
        ebuck6d *= sme;
      }
    } else forcebuck6d = 0.0;

    // add forces of additional LJ interaction

    fpair = forcebuck6d * r32inv;
    if (newton_pair || i1 < nlocal) {
      f[i1][0] += delx3*fpair;
      f[i1][1] += dely3*fpair;
      f[i1][2] += delz3*fpair;
    }
    if (newton_pair || i3 < nlocal) {
      f[i3][0] -= delx3*fpair;
      f[i3][1] -= dely3*fpair;
      f[i3][2] -= delz3*fpair;
    }

    evdwl = 0.0;
    if (eflag) {
      if (rsq3 < cut_ljsq[itype][jtype]) {
        evdwl = ebuck6d - offset[itype][jtype];
      }
    }

    //update pair energy and velocities

    if (evflag) force->pair->ev_tally(i1,i3,nlocal,newton_pair,
                                      evdwl,0.0,fpair,delx3,dely3,delz3);

    tk = multiplicity[type]*acos(c)-th0[type];

    if (eflag) eangle = k[type]*(1.0+cos(tk));

    a = k[type]*multiplicity[type]*sin(tk)*s;

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

void AngleCosineBuck6d::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(multiplicity,n+1,"angle:multiplicity");
  memory->create(th0,n+1,"angle:th0");
  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleCosineBuck6d::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  double c_one = force->numeric(FLERR,arg[1]);
  int n_one = force->inumeric(FLERR,arg[2]);
  int th0_one = force->numeric(FLERR,arg[3]);
  if (n_one <= 0) error->all(FLERR,"Incorrect args for angle coefficients");


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = c_one;
    multiplicity[i] = n_one;
    // transform offset angle to radians
    th0[i] = th0_one/180.0 * MY_PI;

    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ----------------------------------------------------------------------
   check special_bond settings are valid and initialize vdwl parameters
------------------------------------------------------------------------- */

void AngleCosineBuck6d::init_style()
{
  // set local ptrs to buck6d 13 arrays setup by Pair
  int itmp;
  if (force->pair == NULL)
    error->all(FLERR,"Angle cosine/buck6d is incompatible with Pair style");
  cut_ljsq = (double **) force->pair->extract("cut_ljsq",itmp);
  buck6d1 = (double **) force->pair->extract("buck6d1",itmp);
  buck6d2 = (double **) force->pair->extract("buck6d2",itmp);
  buck6d3 = (double **) force->pair->extract("buck6d3",itmp);
  buck6d4 = (double **) force->pair->extract("buck6d4",itmp);
  rsmooth_sq = (double **) force->pair->extract("rsmooth_sq",itmp);
  c0 = (double **) force->pair->extract("c0",itmp);
  c1 = (double **) force->pair->extract("c1",itmp);
  c2 = (double **) force->pair->extract("c2",itmp);
  c3 = (double **) force->pair->extract("c3",itmp);
  c4 = (double **) force->pair->extract("c4",itmp);
  c5 = (double **) force->pair->extract("c5",itmp);
  offset = (double **) force->pair->extract("offset",itmp);
  if (!buck6d1 || !buck6d2 || !buck6d3 || !buck6d4 || !c0 || !c1 || !c2)
    error->all(FLERR,"Angle cosine/buck6d is incompatible with Pair style");

  // special bonds must be x 0 x or double counting
  if (force->special_lj[2] != 0.0)
    error->all(FLERR,"Angle style requires special_bonds lj = x,0,x;"
                     " otherwise buck6d 1-3 interaction are counted twice");
}

/* ---------------------------------------------------------------------- */

double AngleCosineBuck6d::equilibrium_angle(int /*i*/)
{
  return MY_PI;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCosineBuck6d::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&multiplicity[1],sizeof(int),atom->nangletypes,fp);
  fwrite(&th0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleCosineBuck6d::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nangletypes,fp);
    fread(&multiplicity[1],sizeof(int),atom->nangletypes,fp);
    fread(&th0[1],sizeof(double),atom->nangletypes,fp);
  }

  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&multiplicity[1],atom->nangletypes,MPI_INT,0,world);
  MPI_Bcast(&th0[1],atom->nangletypes,MPI_INT,0,world);
  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void AngleCosineBuck6d::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nangletypes; i++) {
    fprintf(fp,"%d %g %d %g\n",i,k[i],multiplicity[i],th0[i]);
  }
}

/* ---------------------------------------------------------------------- */

double AngleCosineBuck6d::single(int type, int i1, int i2, int i3)
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

  double tk = multiplicity[type]*acos(c)-th0[type];

  return k[type]*(1.0+cos(tk));
}
