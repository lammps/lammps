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

#include "improper_harmonic.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperHarmonic::ImproperHarmonic(LAMMPS *lmp) : Improper(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperHarmonic::~ImproperHarmonic()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(chi);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperHarmonic::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double ss1,ss2,ss3,r1,r2,r3,c0,c1,c2,s1,s2;
  double s12,c,s,domega,a,a11,a22,a33,a12,a13,a23;
  double sx2,sy2,sz2;

  eimproper = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nimproperlist; n++) {
    i1 = improperlist[n][0];
    i2 = improperlist[n][1];
    i3 = improperlist[n][2];
    i4 = improperlist[n][3];
    type = improperlist[n][4];

    // geometry of 4-body

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    ss1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
    ss2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
    ss3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

    r1 = sqrt(ss1);
    r2 = sqrt(ss2);
    r3 = sqrt(ss3);

    // sin and cos of angle

    c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
    c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
    c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

    s1 = 1.0 - c1*c1;
    if (s1 < SMALL) s1 = SMALL;
    s1 = 1.0 / s1;

    s2 = 1.0 - c2*c2;
    if (s2 < SMALL) s2 = SMALL;
    s2 = 1.0 / s2;

    s12 = sqrt(s1*s2);
    c = (c1*c2 + c0) * s12;

    // error check

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
        char str[128];
        sprintf(str,"Improper problem: %d " BIGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT,
                me,update->ntimestep,
                atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
        error->warning(FLERR,str,0);
        fprintf(screen,"  1st atom: %d %g %g %g\n",
                me,x[i1][0],x[i1][1],x[i1][2]);
        fprintf(screen,"  2nd atom: %d %g %g %g\n",
                me,x[i2][0],x[i2][1],x[i2][2]);
        fprintf(screen,"  3rd atom: %d %g %g %g\n",
                me,x[i3][0],x[i3][1],x[i3][2]);
        fprintf(screen,"  4th atom: %d %g %g %g\n",
                me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    s = sqrt(1.0 - c*c);
    if (s < SMALL) s = SMALL;

    // force & energy

    domega = acos(c) - chi[type];
    a = k[type] * domega;

    if (eflag) eimproper = a*domega;

    a = -a * 2.0/s;
    c = c * a;
    s12 = s12 * a;
    a11 = c*ss1*s1;
    a22 = -ss2 * (2.0*c0*s12 - c*(s1+s2));
    a33 = c*ss3*s2;
    a12 = -r1*r2*(c1*c*s1 + c2*s12);
    a13 = -r1*r3*s12;
    a23 = r2*r3*(c2*c*s2 + c1*s12);

    sx2  = a22*vb2x + a23*vb3x + a12*vb1x;
    sy2  = a22*vb2y + a23*vb3y + a12*vb1y;
    sz2  = a22*vb2z + a23*vb3z + a12*vb1z;

    f1[0] = a12*vb2x + a13*vb3x + a11*vb1x;
    f1[1] = a12*vb2y + a13*vb3y + a11*vb1y;
    f1[2] = a12*vb2z + a13*vb3z + a11*vb1z;

    f2[0] = -sx2 - f1[0];
    f2[1] = -sy2 - f1[1];
    f2[2] = -sz2 - f1[2];

    f4[0] = a23*vb2x + a33*vb3x + a13*vb1x;
    f4[1] = a23*vb2y + a33*vb3y + a13*vb1y;
    f4[2] = a23*vb2z + a33*vb3z + a13*vb1z;

    f3[0] = sx2 - f4[0];
    f3[1] = sy2 - f4[1];
    f3[2] = sz2 - f4[2];

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f1,f3,f4,
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperHarmonic::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(k,n+1,"improper:k");
  memory->create(chi,n+1,"improper:chi");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperHarmonic::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double chi_one = force->numeric(FLERR,arg[2]);

  // convert chi from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    chi[i] = chi_one/180.0 * MY_PI;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperHarmonic::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&chi[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperHarmonic::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&chi[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&k[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&chi[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperHarmonic::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],chi[i]/MY_PI*180.0);
}
