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
   Contributing author: Tod A Pascal (Caltech)
------------------------------------------------------------------------- */

#include "improper_umbrella.h"
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
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

ImproperUmbrella::ImproperUmbrella(LAMMPS *lmp) : Improper(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

ImproperUmbrella::~ImproperUmbrella()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(kw);
    memory->destroy(w0);
    memory->destroy(C);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperUmbrella::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double domega,c,a,s,projhfg,dhax,dhay,dhaz,dahx,dahy,dahz,cotphi;
  double ax,ay,az,ra2,rh2,ra,rh,rar,rhr,arx,ary,arz,hrx,hry,hrz;

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

    // 1st bond

    vb1x = x[i2][0] - x[i1][0];
    vb1y = x[i2][1] - x[i1][1];
    vb1z = x[i2][2] - x[i1][2];

    // 2nd bond

    vb2x = x[i3][0] - x[i1][0];
    vb2y = x[i3][1] - x[i1][1];
    vb2z = x[i3][2] - x[i1][2];

    // 3rd bond

    vb3x = x[i4][0] - x[i1][0];
    vb3y = x[i4][1] - x[i1][1];
    vb3z = x[i4][2] - x[i1][2];

    // c0 calculation
    // A = vb1 X vb2 is perpendicular to IJK plane

    ax = vb1y*vb2z-vb1z*vb2y;
    ay = vb1z*vb2x-vb1x*vb2z;
    az = vb1x*vb2y-vb1y*vb2x;
    ra2 = ax*ax+ay*ay+az*az;
    rh2 = vb3x*vb3x+vb3y*vb3y+vb3z*vb3z;
    ra = sqrt(ra2);
    rh = sqrt(rh2);
    if (ra < SMALL) ra = SMALL;
    if (rh < SMALL) rh = SMALL;

    rar = 1/ra;
    rhr = 1/rh;
    arx = ax*rar;
    ary = ay*rar;
    arz = az*rar;
    hrx = vb3x*rhr;
    hry = vb3y*rhr;
    hrz = vb3z*rhr;

    c = arx*hrx+ary*hry+arz*hrz;

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
    cotphi = c/s;

    projhfg = (vb3x*vb1x+vb3y*vb1y+vb3z*vb1z) /
      sqrt(vb1x*vb1x+vb1y*vb1y+vb1z*vb1z);
    projhfg += (vb3x*vb2x+vb3y*vb2y+vb3z*vb2z) /
      sqrt(vb2x*vb2x+vb2y*vb2y+vb2z*vb2z);
    if (projhfg > 0.0) {
      s *= -1.0;
      cotphi *= -1.0;
    }

    //  force and energy
    // if w0 = 0: E = k * (1 - cos w)
    // if w0 != 0: E = 0.5 * C (cos w - cos w0)^2, C = k/(sin(w0)^2

    if (w0[type] == 0.0) {
      if (eflag) eimproper = kw[type] * (1.0-s);
      a = -kw[type];
    } else {
      domega = s - cos(w0[type]);
      a = 0.5 * C[type] * domega;
      if (eflag) eimproper = a * domega;
      a *= 2.0;
    }

    // dhax = diffrence between H and A in X direction, etc

    a = a*cotphi;
    dhax = hrx-c*arx;
    dhay = hry-c*ary;
    dhaz = hrz-c*arz;

    dahx = arx-c*hrx;
    dahy = ary-c*hry;
    dahz = arz-c*hrz;

    f2[0] = (dhay*vb1z - dhaz*vb1y)*rar*a;
    f2[1] = (dhaz*vb1x - dhax*vb1z)*rar*a;
    f2[2] = (dhax*vb1y - dhay*vb1x)*rar*a;

    f3[0] = (-dhay*vb2z + dhaz*vb2y)*rar*a;
    f3[1] = (-dhaz*vb2x + dhax*vb2z)*rar*a;
    f3[2] = (-dhax*vb2y + dhay*vb2x)*rar*a;

    f4[0] = dahx*rhr*a;
    f4[1] = dahy*rhr*a;
    f4[2] = dahz*rhr*a;

    f1[0] = -(f2[0] + f3[0] + f4[0]);
    f1[1] = -(f2[1] + f3[1] + f4[1]);
    f1[2] = -(f2[2] + f3[2] + f4[2]);

    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f3[0];
      f[i2][1] += f3[1];
      f[i2][2] += f3[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f2[0];
      f[i3][1] += f2[1];
      f[i3][2] += f2[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag) {

      // get correct 4-body geometry for virial tally

      vb1x = x[i1][0] - x[i2][0];
      vb1y = x[i1][1] - x[i2][1];
      vb1z = x[i1][2] - x[i2][2];

      vb2x = x[i3][0] - x[i2][0];
      vb2y = x[i3][1] - x[i2][1];
      vb2z = x[i3][2] - x[i2][2];

      vb3x = x[i4][0] - x[i3][0];
      vb3y = x[i4][1] - x[i3][1];
      vb3z = x[i4][2] - x[i3][2];

      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f1,f2,f4,
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
    }
  }
}

/* ---------------------------------------------------------------------- */

void ImproperUmbrella::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(kw,n+1,"improper:kw");
  memory->create(w0,n+1,"improper:w0");
  memory->create(C,n+1,"improper:C");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperUmbrella::coeff(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Incorrect args for improper coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double w_one = force->numeric(FLERR,arg[2]);

  // convert w0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kw[i] = k_one;
    w0[i] = w_one/180.0 * MY_PI;
    if (w_one == 0) C[i] = 1.0;
    else C[i] = kw[i]/(pow(sin(w0[i]),2.0));
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperUmbrella::write_restart(FILE *fp)
{
  fwrite(&kw[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&w0[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&C[1],sizeof(double),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperUmbrella::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&kw[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&w0[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&C[1],sizeof(double),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&kw[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&w0[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&C[1],atom->nimpropertypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperUmbrella::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g %g\n",i,kw[i],w0[i]/MY_PI*180.0);
}
