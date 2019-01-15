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
   Contributing author: Loukas D. Peristeras (Scienomics SARL)
   [ based on improper_umbrella.cpp Tod A Pascal (Caltech) ]
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "improper_fourier.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

ImproperFourier::ImproperFourier(LAMMPS *lmp) : Improper(lmp) {}

/* ---------------------------------------------------------------------- */

ImproperFourier::~ImproperFourier()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(C0);
    memory->destroy(C1);
    memory->destroy(C2);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperFourier::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;

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

    addone(i1,i2,i3,i4, type,evflag,eflag,
           vb1x, vb1y, vb1z,
           vb2x, vb2y, vb2z,
           vb3x, vb3y, vb3z);
    if ( all[type] ) {
      addone(i1,i4,i2,i3, type,evflag,eflag,
             vb3x, vb3y, vb3z,
             vb1x, vb1y, vb1z,
             vb2x, vb2y, vb2z);
      addone(i1,i3,i4,i2, type,evflag,eflag,
             vb2x, vb2y, vb2z,
             vb3x, vb3y, vb3z,
             vb1x, vb1y, vb1z);
    }
  }
}

void ImproperFourier::addone(const int &i1,const int &i2,const int &i3,const int &i4,
            const int &type,const int &evflag,const int &eflag,
            const double &vb1x, const double &vb1y, const double &vb1z,
            const double &vb2x, const double &vb2y, const double &vb2z,
            const double &vb3x, const double &vb3y, const double &vb3z)
{
  double eimproper,f1[3],f2[3],f3[3],f4[3];
  double c,c2,a,s,projhfg,dhax,dhay,dhaz,dahx,dahy,dahz,cotphi;
  double ax,ay,az,ra2,rh2,ra,rh,rar,rhr,arx,ary,arz,hrx,hry,hrz;

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  eimproper = 0.0;

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
  //  E = k ( C0 + C1 cos(w) + C2 cos(w)

  c2 = 2.0*s*s-1.0;
  if (eflag) eimproper =  k[type]*(C0[type]+C1[type]*s+C2[type]*c2);

  // dhax = diffrence between H and A in X direction, etc

  a = k[type]*(C1[type]+4.0*C2[type]*s)*cotphi;
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

  if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,eimproper,f1,f2,f4,
               -vb1x,-vb1y,-vb1z,vb2x-vb1x,vb2y-vb1y,vb2z-vb1z,vb3x-vb2x,vb3y-vb2y,vb3z-vb2z);
}

/* ---------------------------------------------------------------------- */

void ImproperFourier::allocate()
{
  allocated = 1;
  int n = atom->nimpropertypes;

  memory->create(k,n+1,"improper:k");
  memory->create(C0,n+1,"improper:C0");
  memory->create(C1,n+1,"improper:C1");
  memory->create(C2,n+1,"improper:C2");
  memory->create(all,n+1,"improper:C2");

  memory->create(setflag,n+1,"improper:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void ImproperFourier::coeff(int narg, char **arg)
{

  if ( narg != 5 && narg != 6 ) error->all(FLERR,"Incorrect args for improper coefficients");

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nimpropertypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double C0_one = force->numeric(FLERR,arg[2]);
  double C1_one = force->numeric(FLERR,arg[3]);
  double C2_one = force->numeric(FLERR,arg[4]);
  int all_one = 1;
  if ( narg == 6 ) all_one = force->inumeric(FLERR,arg[5]);

  // convert w0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    C0[i] = C0_one;
    C1[i] = C1_one;
    C2[i] = C2_one;
    all[i] = all_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for improper coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void ImproperFourier::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&C0[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&C1[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&C2[1],sizeof(double),atom->nimpropertypes,fp);
  fwrite(&all[1],sizeof(int),atom->nimpropertypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void ImproperFourier::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&C0[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&C1[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&C2[1],sizeof(double),atom->nimpropertypes,fp);
    fread(&all[1],sizeof(int),atom->nimpropertypes,fp);
  }
  MPI_Bcast(&k[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&C0[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&C1[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&C2[1],atom->nimpropertypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&all[1],atom->nimpropertypes,MPI_INT,0,world);

  for (int i = 1; i <= atom->nimpropertypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void ImproperFourier::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp,"%d %g %g %g %g %d\n",i,k[i],C0[i],C1[i],C2[i],all[i]);
}
