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
   Special Angle Potential for the CMM coarse grained MD potentials.
   Contributing author: Axel Kohlmeyer <akohlmey@gmail.com>
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "angle_cg_cmm.h"
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

AngleCGCMM::AngleCGCMM(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleCGCMM::~AngleCGCMM()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(theta0);
    memory->destroy(cg_type);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(rcut);
  }
}

/* ---------------------------------------------------------------------- */

void AngleCGCMM::ev_tally_lj13(int i, int j, int nlocal, int newton_bond,
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

void AngleCGCMM::compute(int eflag, int vflag)
{
  int i1,i2,i3,n,type;
  double delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3;
  double eangle,f1[3],f3[3],e13,f13;
  double dtheta,tk;
  double rsq1,rsq2,rsq3,r1,r2,r3,c,s,a,a11,a12,a22;

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
    // so this has to be done here and not in the
    // general non-bonded code.
    delx3 = x[i1][0] - x[i3][0];
    dely3 = x[i1][1] - x[i3][1];
    delz3 = x[i1][2] - x[i3][2];
    rsq3 = delx3*delx3 + dely3*dely3 + delz3*delz3;
    r3 = sqrt(rsq3);

    f13=0.0;
    e13=0.0;

    if (r3 < rcut[type]) {
      const int cgt = cg_type[type];
      const double cgpow1 = cg_pow1[cgt];
      const double cgpow2 = cg_pow2[cgt];
      const double cgpref = cg_prefact[cgt];

      const double ratio = sigma[type]/r3;
      const double eps = epsilon[type];

      f13 = cgpref*eps / rsq3 * (cgpow1*pow(ratio,cgpow1)
                                  - cgpow2*pow(ratio,cgpow2));

      if (eflag) e13 = eps + cgpref*eps * (pow(ratio,cgpow1)
                                          - pow(ratio,cgpow2));

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

    // apply force to each of 3 atoms

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

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);

    if (evflag) ev_tally_lj13(i1,i3,nlocal,newton_bond,
                         e13,f13,delx3,dely3,delz3);

  }
}

/* ---------------------------------------------------------------------- */

void AngleCGCMM::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(k,n+1,"angle:k");
  memory->create(theta0,n+1,"angle:theta0");
  memory->create(epsilon,n+1,"angle:epsilon");
  memory->create(sigma,n+1,"angle:sigma");
  memory->create(rcut,n+1,"angle:rcut");

  memory->create(cg_type,n+1,"angle:cg_type");
  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) {
    cg_type[i] = CG_NOT_SET;
    setflag[i] = 0;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleCGCMM::coeff(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double theta0_one = force->numeric(FLERR,arg[2]);

  int cg_type_one=find_cg_type(arg[3]);
  if (cg_type_one == CG_NOT_SET) error->all(FLERR,"Error reading CG type flag.");

  double epsilon_one = force->numeric(FLERR,arg[4]);
  double sigma_one = force->numeric(FLERR,arg[5]);

  // find minimum of LJ potential. we only want to include
  // the repulsive part of the 1-3 LJ.
  double rcut_one = sigma_one*exp(
           1.0/(cg_pow1[cg_type_one]-cg_pow2[cg_type_one])
           *log(cg_pow1[cg_type_one]/cg_pow2[cg_type_one])
         );

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    // convert theta0 from degrees to radians
    theta0[i] = theta0_one/180.0 * MY_PI;
    epsilon[i] = epsilon_one;
    sigma[i] = sigma_one;
    rcut[i] = rcut_one;
    cg_type[i] = cg_type_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ---------------------------------------------------------------------- */

double AngleCGCMM::equilibrium_angle(int i)
{
  return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleCGCMM::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&epsilon[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&sigma[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&rcut[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&cg_type[1],sizeof(int),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleCGCMM::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nangletypes,fp);
    fread(&theta0[1],sizeof(double),atom->nangletypes,fp);
    fread(&epsilon[1],sizeof(double),atom->nangletypes,fp);
    fread(&sigma[1],sizeof(double),atom->nangletypes,fp);
    fread(&rcut[1],sizeof(double),atom->nangletypes,fp);
    fread(&cg_type[1],sizeof(int),atom->nangletypes,fp);
  }
  MPI_Bcast(&k[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsilon[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&rcut[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&cg_type[1],atom->nangletypes,MPI_INT,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double AngleCGCMM::single(int type, int i1, int i2, int i3)
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

  // 1-3 LJ interaction.
  double delx3 = x[i1][0] - x[i3][0];
  double dely3 = x[i1][1] - x[i3][1];
  double delz3 = x[i1][2] - x[i3][2];
  domain->minimum_image(delx3,dely3,delz3);

  const double r3 = sqrt(delx3*delx3 + dely3*dely3 + delz3*delz3);

  double e13=0.0;

  if (r3 < rcut[type]) {
    const int cgt = cg_type[type];
    const double cgpow1 = cg_pow1[cgt];
    const double cgpow2 = cg_pow2[cgt];
    const double cgpref = cg_prefact[cgt];

    const double ratio = sigma[type]/r3;
    const double eps = epsilon[type];

    e13 = eps + cgpref*eps * (pow(ratio,cgpow1)
                              - pow(ratio,cgpow2));
  }

  double dtheta = acos(c) - theta0[type];
  double tk = k[type] * dtheta;
  return tk*dtheta + e13;
}
